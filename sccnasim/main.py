# main.py - cmdline interface.



import logging
import numpy as np
import os
import sys
import time

from logging import info, error
from .afc.main import afc_wrapper
from .app import APP, VERSION
from .config import Config
from .cs.main import cs_wrapper
from .pp.main import pp_wrapper
from .rs.main import rs_wrapper
from .utils.cellanno import load_cells
from .utils.gfeature import assign_feature_batch, load_features
from .xlib.xbase import assert_e
from .xlib.xdata import load_h5ad, save_h5ad
from .xlib.xlog import init_logging



def main():
    pass



def main_wrapper(
    sam_fn,
    cell_anno_fn, feature_fn, phased_snp_fn,
    clone_anno_fn, cna_profile_fn, 
    refseq_fn,
    out_dir,
    overlap_features_how = "raw",
    size_factor = "libsize",
    marginal = "auto",
    libsize_ratio = 1.0,
    loss_allele_freq = 0.01,
    kwargs_fit_sf = None,
    kwargs_fit_rd = None,
    chroms = "human_autosome",
    cell_tag = "CB", umi_tag = "UB", umi_len = 10,
    barcode_whitelist_fn = None,
    ncores = 1, seed = 123, verbose = False,
    min_count = 1, min_maf = 0,
    strandness = "forward", min_include = 0.5, multi_mapper_how = "discard",
    xf_tag = "xf", gene_tag = "GN",
    min_mapq = 20, min_len = 30,
    incl_flag = 0, excl_flag = -1,
    no_orphan = True,
    debug_level = 0
):
    """Wrapper for running the main pipeline.

    Parameters
    ----------
    sam_fn : str
        Indexed BAM file.
    cell_anno_fn : str
        The cell annotation file. 
        It is a header-free TSV file and its first two columns are:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
    feature_fn : str
        A TSV file listing target features. 
        It is header-free and its first 5 columns shoud be: 
        - "chrom" (str): chromosome name of the feature.
        - "start" (int): start genomic position of the feature, 1-based
          and inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
        - "strand" (str): feature strand, either "+" (positive) or 
          "-" (negative).
    phased_snp_fn : str
        A TSV or VCF file listing phased SNPs.
        If TSV, it is a header-free file containing SNP annotations, whose
        first six columns should be:
        - "chrom" (str): chromosome name of the SNP.
        - "pos" (int): genomic position of the SNP, 1-based.
        - "ref" (str): the reference allele of the SNP.
        - "alt" (str): the alternative allele of the SNP.
        - "ref_hap" (int): the haplotype index of `ref`, one of {0, 1}.
        - "alt_hap" (int): the haplotype index of `alt`, one of {1, 0}.
        If VCF, it should contain "GT" in its "FORMAT" field.
    clone_anno_fn : str
        A TSV file listing clonal anno information.
        It is header-free and its first 3 columns are:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the source cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, 
          then it will be set as the number of cells in `source_cell_type`.
    cna_profile_fn : str
        A TSV file listing clonal CNA profiles. 
        It is header-free and its first 6 columns are:
        - "chrom" (str): chromosome name of the CNA region.
        - "start" (int): start genomic position of the CNA region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNA region, 1-based and
          inclusive.
        - "clone" (str): clone ID.
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
    refseq_fn : str
        A FASTA file storing reference genome sequence.
    out_dir : str
        The output folder.
    overlap_features_how : str, default "raw"
        How to process overlapping features.
        - "raw": Leave all input gene annotations unchanged.
        - "quantile": remove highly overlapping genes.
           Remove genes with number of overlapping genes larger than a given
           value (default is the 0.99 quantile among all genes that have 
           overlaps).
        - "union": keep the union range of gene overlaps.
           Replace consecutive overlapping genes with their union genomic 
           range, i.e., aggregate overlapping genes into non-overlapping
           super-genes.
    size_factor : str or None, default "libsize"
        The type of size factor.
        Currently, only support "libsize" (library size).
        Set to `None` if do not use size factors for model fitting.
    marginal : {"auto", "poi", "nb", "zinb"}
        Type of marginal distribution.
        One of
        - "auto" (auto select).
        - "poi" (Poisson).
        - "nb" (Negative Binomial).
        - "zinb" (Zero-Inflated Negative Binomial).
    libsize_ratio : float, default 1.0
        Ratio of library size of simulated cells compared to seed cells.
    loss_allele_freq : float, default 0.01
        The frequency of the lost allele, to mimic real error rate, i.e.,
        sometimes we observe reads from the lost allele.
    kwargs_fit_sf : dict or None, default None
        The additional kwargs passed to function 
        :func:`~.marginal.fit_libsize_wrapper` for fitting size factors.
        The available arguments are:
        - dist : {"lognormal", "swr", "normal", "t"}
            Type of distribution.
        If None, set to `{}`.
    kwargs_fit_rd : dict or None, default None
        The additional kwargs passed to function 
        :func:`~.marginal.fit_RD_wrapper` for fitting read depth.
        The available arguments are:
        - min_nonzero_num : tuple of int, default (1, 1, 3)
            The minimum number of cells that have non-zeros in one feature,
            for alleles 'A', 'B', and 'U', respectively.
            If smaller than the cutoff, then the feature will not be fitted
            (i.e., its mean will be directly treated as 0).
        - max_iter : int, default 1000
            Number of maximum iterations in model fitting.
        - pval_cutoff : float, default 0.05
            The p-value cutoff for model selection with GLR test.
        If None, set to `{}`.
    chroms : str, default "human_autosome"
        Comma separated chromosome names.
        Reads in other chromosomes will not be used for sampling and hence
        will not be present in the output BAM file(s).
        If "human_autosome", set to `"1,2,...22"`.
    cell_tag : str or None, default "CB"
        Tag for cell barcodes, set to None when using sample IDs.
    umi_tag : str or None, default "UB"
        Tag for UMI, set to None when reads only.
    umi_len : int, default 10
        Length of output UMI barcode.
    barcode_whitelist_fn : str or None, default None
        File containing whitelist cell barcodes to be sampled for simulated
        data.
        If None, use randomly generated cell barcodes.
    ncores : int, default 1
        Number of cores.
    seed : int or None, default 123
        Seed for random numbers.
        None means not using a fixed seed.
    verbose : bool, default False
        Whether to show detailed logging information.
    min_count : int, default 1
        Minimum aggragated count for SNP.
    min_maf : float, default 0
        Minimum minor allele fraction for SNP.
    strandness : {"forward", "reverse", "unstranded"}
        Strandness of the sequencing protocol.
        - "forward": SE sense; PE R1 antisense and R2 sense;
            e.g., 10x 3' data.
        - "reverse": SE antisense; PE R1 sense and R2 antisense;
            e.g., 10x 5' data.
        - "unstranded": no strand information.
    min_include : int or float, default 0.5
        Minimum length of included part within specific feature.
        If float between (0, 1), it is the minimum fraction of included length.
    multi_mapper_how : {"discard", "duplicate"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "duplicate": count the UMI for every mapped gene.
    xf_tag : str or None, default "xf"
        The extra alignment flags set by tools like CellRanger or SpaceRanger.
        If set, only reads with tag's value 17 or 25 will count.
        If `None`, turn this tag off.
    gene_tag : str or None, default "GN"
        The tag for gene name set by tools like CellRanger or SpaceRanger.
        If `None`, turn this tag off.
    min_mapq : int, default 20
        Minimum MAPQ for read filtering.
    min_len : int, default 30
        Minimum mapped length for read filtering.
    incl_flag : int, default 0
        Required flags: skip reads with all mask bits unset.
    excl_flag : int, default -1
        Filter flags: skip reads with any mask bits set.
        Value -1 means setting it to 772 when using UMI, or 1796 otherwise.
    no_orphan : bool, default True
        If `False`, do not skip anomalous read pairs.

    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    conf = Config()


    # input and output files.
    conf.sam_fn = sam_fn
    conf.cell_anno_fn = cell_anno_fn
    conf.feature_fn = feature_fn
    conf.snp_fn = phased_snp_fn
    conf.clone_anno_fn = clone_anno_fn
    conf.cna_profile_fn = cna_profile_fn
    conf.refseq_fn = refseq_fn
    conf.out_dir = out_dir
    
    
    # preprocessing.
    conf.overlap_features_how = overlap_features_how


    # count simulation.
    conf.size_factor = size_factor
    conf.marginal = marginal
    conf.libsize_ratio = libsize_ratio
    conf.loss_allele_freq = loss_allele_freq
    if kwargs_fit_sf is None:
        conf.kwargs_fit_sf = dict()
    else:
        conf.kwargs_fit_sf = kwargs_fit_sf
    if kwargs_fit_rd is None:
        conf.kwargs_fit_rd = dict()
    else:
        conf.kwargs_fit_rd = kwargs_fit_rd


    # optional arguments.
    if chroms == "human_autosome":
        conf.chroms = ",".join([str(c) for c in range(1, 23)])
    else:
        conf.chroms = chroms
    conf.cell_tag = cell_tag
    conf.umi_tag = umi_tag
    conf.umi_len = umi_len
    conf.barcode_whitelist_fn = barcode_whitelist_fn
    conf.ncores = ncores
    conf.seed = seed
    conf.verbose = verbose
    
    
    # SNP filtering.
    conf.min_count = min_count
    conf.min_maf = min_maf
    
    
    # read assignment.
    conf.strandness = strandness
    conf.min_include = min_include
    conf.multi_mapper_how = multi_mapper_how
    conf.xf_tag = xf_tag
    conf.gene_tag = gene_tag


    # read filtering.
    conf.min_mapq = min_mapq
    conf.min_len = min_len
    conf.incl_flag = incl_flag
    conf.excl_flag = excl_flag
    conf.no_orphan = no_orphan
    
    conf.debug_level = debug_level


    ret, res = main_run(conf)
    return((ret, res))



def main_core(conf):
    ret = main_init(conf)
    if ret < 0:
        raise ValueError
    conf.show(fp = sys.stdout, prefix = "\t")
    os.makedirs(conf.out_dir, exist_ok = True)

    step = 1

    # Note:
    # Use `xx_wrapper()`` function in each step instead of directly accessing
    # or modifying the internal `config` object, to keep codes independent.

    # preprocessing.
    info("start preprocessing ...")
    pp_ret, pp_res = pp_wrapper(
        cell_anno_fn = conf.cell_anno_fn,
        feature_fn = conf.feature_fn,
        snp_fn = conf.snp_fn,
        clone_anno_fn = conf.clone_anno_fn,
        cna_profile_fn = conf.cna_profile_fn,
        out_dir = os.path.join(conf.out_dir, "%d_pp" % step),
        chroms = conf.chroms,
        strandness = conf.strandness,
        overlap_features_how = conf.overlap_features_how
    )
    if pp_ret < 0:
        error("preprocessing failed (%d)." % pp_ret)
        raise ValueError
    info("pp results:")
    info(str(pp_res))
    step += 1
    
    
    # create feature-specific result folders.
    feature_dir = os.path.join(conf.out_dir, "features")
    os.makedirs(feature_dir, exist_ok = True)
    features = load_features(pp_res["feature_fn_new"])
    out_feature_dirs = assign_feature_batch(
        feature_names = features["feature"].values,
        root_dir = feature_dir,
        batch_size = 1000
    )


    # allele-specific feature counting.
    info("start allele-specific feature counting ...")
    afc_ret, afc_res = afc_wrapper(
        sam_fn = conf.sam_fn,
        barcode_fn = pp_res["barcode_fn_new"] if conf.use_barcodes() else None,
        feature_fn = pp_res["feature_fn_new"],
        phased_snp_fn = pp_res["snp_fn_new"],
        out_dir = os.path.join(conf.out_dir, "%d_afc" % step),
        debug_level = conf.debug_level,
        ncores = conf.ncores,
        cell_tag = conf.cell_tag,
        umi_tag = conf.umi_tag,
        min_count = conf.min_count,
        min_maf = conf.min_maf,
        strandness = conf.strandness,
        min_include = conf.min_include,
        multi_mapper_how = conf.multi_mapper_how,
        xf_tag = conf.xf_tag,
        gene_tag = conf.gene_tag,
        min_mapq = conf.min_mapq,
        min_len = conf.min_len,
        incl_flag = conf.incl_flag,
        excl_flag = conf.excl_flag,
        no_orphan = conf.no_orphan,
        out_feature_dirs = out_feature_dirs
    )
    if afc_ret < 0:
        error("allele-specific feature counting failed (%d)." % afc_ret)
        raise ValueError
    info("afc results:")
    info(str(afc_res))
    step += 1


    # count simulation.
    info("start count simulation ...")
    count_fn_new = afc_res["count_fn"].replace(".h5ad", ".cell_anno.h5ad")
    add_cell_anno(
        count_fn = afc_res["count_fn"],
        cell_anno_fn = pp_res["cell_anno_fn_new"],
        out_count_fn = count_fn_new
    )
    info("new input (annotated) count adata file is saved to '%s'." % \
        count_fn_new)

    cs_ret, cs_res = cs_wrapper(
        count_fn = count_fn_new,
        clone_anno_fn = pp_res["clone_anno_fn_new"],
        cna_profile_fn = pp_res["cna_profile_fn_new"],
        out_dir = os.path.join(conf.out_dir, "%d_cs" % step),
        size_factor = conf.size_factor,
        marginal = conf.marginal,
        libsize_ratio = conf.libsize_ratio,
        loss_allele_freq = conf.loss_allele_freq,
        cna_mode = "hap-aware",
        barcode_whitelist_fn = conf.barcode_whitelist_fn,
        ncores = conf.ncores,
        verbose = conf.verbose,
        kwargs_fit_sf = conf.kwargs_fit_sf,
        kwargs_fit_rd = conf.kwargs_fit_rd
    )
    if cs_ret < 0:
        error("count simulation failed (%d)." % cs_ret)
        raise ValueError
    info("cs results:")
    info(str(cs_res))
    step += 1


    # read simulation.
    info("start read simulation ...")
    rs_ret, rs_res = rs_wrapper(
        count_fn = cs_res["count_fn"],
        feature_fn = afc_res["fet_obj_fn"],
        refseq_fn = conf.refseq_fn,
        out_dir = os.path.join(conf.out_dir, "%d_rs" % step),
        debug_level = conf.debug_level,
        ncores = conf.ncores,
        cell_tag = conf.cell_tag,
        umi_tag = conf.umi_tag,
        umi_len = conf.umi_len
    )
    if rs_ret < 0:
        error("read simulation failed (%d)." % rs_ret)
        raise ValueError
    info("rs results:")
    info(str(rs_res))
    step += 1


    # construct returned values.
    res = rs_res
    return(res)



def main_run(conf):
    if conf.debug_level > 0:
        init_logging(stream = sys.stdout, ch_level = logging.DEBUG)
    else:
        init_logging(stream = sys.stdout)
        
    # currently the whole simulation results are not reproducible with a seed,
    # possibly due to the parallel computing.
    # TODO: make it reproducible
    # Ref:
    # - https://albertcthomas.github.io/good-practices-random-number-generators/
    # - https://numpy.org/doc/stable/reference/random/parallel.html
    if conf.seed is not None:
        np.random.seed(conf.seed)
    
    
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)
    info("%s (VERSION %s)." % (APP, VERSION))

    try:
        res = main_core(conf)
    except ValueError as e:
        error(str(e))
        error("Running program failed.")
        error("Quiting ...")
        ret = -1
    else:
        info("All Done!")
        ret = 0
    finally:
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return((ret, res))



def main_init(conf):
    if conf.sam_fn is not None:
        assert_e(conf.sam_fn)
    assert_e(conf.cell_anno_fn)
    assert_e(conf.feature_fn)
    assert_e(conf.snp_fn)
    assert_e(conf.clone_anno_fn)
    assert_e(conf.cna_profile_fn)
    assert_e(conf.refseq_fn)
    
    assert conf.strandness in ("forward", "reverse", "unstranded")
    
    assert conf.multi_mapper_how in ("discard", "duplicate")

    if conf.barcode_whitelist_fn is not None:
        if conf.barcode_whitelist_fn.lower() == "none":
            conf.barcode_whitelist_fn = None
        else:
            assert_e(conf.barcode_whitelist_fn)
    
    return(0)



def add_cell_anno(count_fn, cell_anno_fn, out_count_fn):
    cell_anno = load_cells(cell_anno_fn)
    assert "cell" in cell_anno.columns
    assert "cell_type" in cell_anno.columns

    adata = load_h5ad(count_fn)
    assert "cell" in adata.obs.columns

    assert np.all(adata.obs["cell"].isin(cell_anno["cell"]))
    adata.obs = adata.obs.merge(cell_anno, how = "left", on = "cell")

    assert "cell_type" in adata.obs.columns

    save_h5ad(adata, out_count_fn)
