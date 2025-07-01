# core.py - core part of allele-specific feature counting.


import gc
import math
import numpy as np
import os
import pysam

from logging import debug, error, info
from .mcount import MCount, SNPCount
from ..utils.gfeature import load_feature_objects, save_feature_objects
from ..utils.hapidx import hap2idx, idx2hap
from ..utils.sam import check_read
from ..xlib.xfile import zopen, ZF_F_GZIP
from ..xlib.xsam import sam_fetch


# NOTE: 
# 1. bgzf errors when using pysam.AlignmentFile.fetch in parallel (with
#    multiprocessing): https://github.com/pysam-developers/pysam/issues/397


# solved.
# - issue of double counting SNPs in one UMI.
# - post-QC check of orphan reads in PE.
# - UMI collapsing when determine allele of one SNP.



def fc_features(
    reg_obj_fn,
    sam_fn,
    out_mtx_fns,
    samples,
    batch_idx,
    conf
):
    """Feature counting for a list of features.

    This function does feature counting for a list of features. When iterating
    one feature, it 
    (1) calculates UMI/read counts of this feature in single cells and output 
        them into each allele-specific count matrix file.
    (2) outputs post-filtering reads of this feature into each allele-specific
        BAM file.
    (3) outputs the corresponding pileuped cell and UMI IDs (CUMIs) into each
        allele-specific CUMI file.
    
    Parameters
    ----------
    reg_obj_fn : str
        File storing a list of :class:`~..utils.gfeature.Feautre` objects.
    sam_fn : str
        Input BAM file.
    out_mtx_fns : dict of {str : str}
        Output allele-specific sparse matrix files.
        Keys are alleles and values are sparse matrix files.
    samples : list of str
        A list of sample/cell IDs.
    batch_idx : int
        The index of this batch.
    conf : .config.Config
        The :class:`~.config.Config` object.

    Returns
    -------
    dict
        Results of this batch.
    """
    info("[Batch-%d] start ..." % batch_idx)

    sam = pysam.AlignmentFile(sam_fn, "r", require_index = True)
    reg_list = load_feature_objects(reg_obj_fn)
    alleles = list(out_mtx_fns.keys())
    fp_mtx = {ale: zopen(fn, "wt", ZF_F_GZIP, is_bytes = False) \
                for ale, fn in out_mtx_fns.items()}


    # core part.
    m = float(len(reg_list))
    l = 0                    # fraction of processed genes, used for verbose.
    nr_mtx = {ale:0 for ale in alleles}      # number of records.
    for idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            debug("[Batch-%d] processing feature '%s' ..." % \
                  (batch_idx, reg.name))

        r, cnt = fc_fet1(reg, sam, samples, alleles, conf)
        if r < 0 or cnt is None:
            error("errcode -9 (%s)." % reg.name)
            raise RuntimeError

        sr = {ale:"" for ale in alleles}              # string of one record.
        for i, cell in enumerate(samples):
            nu = {ale:cnt[ale][cell] for ale in alleles}     # number of UMIs.
            if np.sum([v for v in nu.values()]) <= 0:
                continue
            for ale in alleles:
                if nu[ale] > 0:
                    sr[ale] += "%d\t%d\t%d\n" % (idx + 1, i + 1, nu[ale])
                    nr_mtx[ale] += 1

        if np.any([len(s) > 0 for s in sr.values()]):
            for ale in alleles:
                fp_mtx[ale].write(sr[ale])

        n = idx + 1
        frac = n / m
        if frac - l >= 0.1 or n == m:
            if conf.debug > 0:
                debug("[Batch-%d] %d%% genes processed" % 
                    (batch_idx, math.floor(frac * 100)))
            l = frac

    nr_reg = len(reg_list)

    
    # clean files.
    if conf.debug > 0:
        debug("[Batch-%d] clean files ..." % batch_idx)

    sam.close()
    for ale in alleles:
        fp_mtx[ale].close()

    
    # reg objects, each containing post-filtering SNPs.
    save_feature_objects(reg_list, reg_obj_fn)
    
    info("[Batch-%d] done!" % batch_idx)


    del reg_list
    gc.collect()

    res = dict(
        # nr_reg : int
        #   Number of unique features in this batch that are outputted.
        nr_reg = nr_reg,
        
        # nr_mtx : dict of {str : int}
        #   Number of records in each allele-specific count matrix file.
        #   Keys are allele names, values are number of records.
        nr_mtx = nr_mtx,
        
        # out_mtx_fns : dict of {str : str}
        #   Output allele-specific sparse matrix files.
        #   Keys are alleles and values are sparse matrix files.
        out_mtx_fns = out_mtx_fns
    )
    
    return(res)



def fc_fet1(reg, sam, samples, alleles, conf):
    """Feature counting for one feature.
    
    Parameters
    ----------
    reg : :class:`~..utils.gfeature.Feature`
        The feature to be counted.
    sam : :class:`pysam.AlignmentFile`
        File object for input SAM/BAM file.
    samples: list of str
        A list of cell barcodes.
    alleles : list of str
        A list of allele names.
    conf : :class:`.config.Config`
        Global configuration object.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    dict or None
        The *allele x cell* counts of this feature.
        It is a two-layer dict, with "allele name (str)" and "cell ID (str)"
        as keys, respectively, and "counts (int)" as values.
        None if error happens or no any UMIs fetched.
    """
    mcnt = MCount(samples, conf)
    out_sam_list = {
        ale: pysam.AlignmentFile(dat.seed_sam_fn, "wb", template = sam) \
            for ale, dat in reg.allele_data.items()
    }

    
    # add UMIs fetched by SNPs into counting machine.
    snp_cnt_list = []
    for snp in reg.snp_list:
        itr = sam_fetch(sam, snp.chrom, snp.pos, snp.pos)
        if not itr:    
            continue
        snp_cnt = SNPCount(snp, mcnt, mcnt.cells, conf)
        for read in itr:
            if check_read(read, reg, conf) < 0:
                continue
            ret = snp_cnt.add_read(read)
            if ret < 0:      # read filtered if ret > 0
                return(-3, None)
        snp_cnt_list.append(snp_cnt)
        
    
    # add all UMIs into counting machine.
    itr = sam_fetch(sam, reg.chrom, reg.start, reg.end - 1)
    if not itr:    
        return(0, None)
    for read in itr:
        if check_read(read, reg, conf) < 0:
            continue
        if conf.use_barcodes():
            ret, ucnt = mcnt.add_read(read)
        else:
            raise ValueError
        if ret < 0:         # read filtered if ret > 0
            return(-5, None)
    
    
    # stat SNP and QC.
    snp_list = []
    for snp_cnt in snp_cnt_list:
        qc_fail = False
        stat_allele = snp_cnt.stat()
        n_all = sum(stat_allele.values())
        if n_all < conf.min_count:
            qc_fail = True
        n_ref = stat_allele[snp.ref]
        n_alt = stat_allele[snp.alt]
        n_minor = min(n_ref, n_alt)
        if n_minor < n_all * conf.min_maf:
            qc_fail = True
        if qc_fail:
            snp_id = snp_cnt.get_id()
            mcnt.del_snp(snp_id)
        else:
            snp_list.append(snp_cnt.snp)
    reg.snp_list = snp_list
        
            
    # stat counting machine.
    stat_hap, stat_umi = mcnt.stat()

    
    # output allele-specific reads.
    hap_idx = None
    itr = sam_fetch(sam, reg.chrom, reg.start, reg.end - 1)
    if not itr:    
        return(-7, None)
    for read in itr:
        if check_read(read, reg, conf) < 0:
            continue
        if conf.use_barcodes():
            hap_idx = mcnt.get_hap_idx(read)
        else:
            raise ValueError
        if hap_idx is None:
            continue   
        read.set_tag(conf.hap_idx_tag, hap_idx)
            
        # output reads to feature-allele-specific SAM/BAM file.
        # Note that these reads are superset of the reads used for read
        # sampling in `rs` module, because after the read iteration loop,
        # there could be a few UMI filtering steps.
        ale = idx2hap(hap_idx)
        if ale not in out_sam_list:
            continue
        out_sam = out_sam_list[ale]
        out_sam.write(read)
                
    for s in out_sam_list.values():
        s.close()
    for ale, dat in reg.allele_data.items():
        pysam.index(dat.seed_sam_fn)


    # output allele-specific CUMI list.
    cumi_fps = {
        ale: open(dat.seed_cumi_fn, "w")   \
            for ale, dat in reg.allele_data.items()
    }
    for ale, fp in cumi_fps.items():
        assert ale in conf.cumi_alleles
        for cell in samples:
            umis = set()
            for i in hap2idx(ale):
                if i not in stat_umi:
                    continue
                if cell not in stat_umi[i]:
                    continue
                umis.update(stat_umi[i][cell])
            for umi in sorted(list(umis)):
                fp.write("%s\t%s\n" % (cell, umi))
        fp.close()
        
        
    # output allele-specific CUMI counts.
    ale_cnt = {ale: {cell:0 for cell in samples} for ale in alleles}
    for ale in alleles:      #  {"A", "B", "D", "O", "U"}
        for cell in samples:
            for i in hap2idx(ale):
                if i not in stat_hap:
                    continue
                if cell not in stat_hap[i]:
                    continue
                ale_cnt[ale][cell] += stat_hap[i][cell]

    return((0, ale_cnt))
