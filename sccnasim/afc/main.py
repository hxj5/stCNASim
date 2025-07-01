# main.py - allele-specific feature counting.


import gc
import multiprocessing
import os
import shutil
import sys
import time

from logging import debug, error, info
from logging import warning as warn
from .config import Config, COMMAND
from .core import fc_features
from .io import load_feature_from_txt,  \
    load_snp_from_vcf, load_snp_from_tsv,  \
    merge_mtx
from .mfu import mfu_main
from .sam import detect_pe_mode
from ..app import APP, VERSION
from ..utils.gfeature import assign_feature_batch,  \
    load_feature_objects, save_feature_objects
from ..utils.io import load_bams, load_barcodes, load_samples
from ..xlib.xbase import assert_e
from ..xlib.xdata import load_adata, save_h5ad
from ..xlib.xfile import ZF_F_GZIP, ZF_F_PLAIN
from ..xlib.xio import load_pickle, save_pickle
from ..xlib.xlog import init_logging
from ..xlib.xthread import split_n2batch, mp_error_handler



def afc_wrapper(
    sam_fn, barcode_fn,
    feature_fn, phased_snp_fn, 
    out_dir,
    debug_level = 0,
    ncores = 1,
    cell_tag = "CB", umi_tag = "UB",
    min_count = 20, min_maf = 0.1,
    strandness = "forward",
    min_include = 0.5,
    multi_mapper_how = "discard",
    xf_tag = "xf", gene_tag = "GN",
    min_mapq = 20, min_len = 30,
    incl_flag = 0, excl_flag = -1,
    no_orphan = True,
    out_feature_dirs = None
):
    """Wrapper for running the afc (allele-specific counting) module.

    Parameters
    ----------
    sam_fn : str
        Input indexed BAM file.
    barcode_fn : str
        A plain file listing all effective cell barcodes.
    feature_fn : str
        A TSV file listing target features.
        It is header-free and its first 4 columns shoud be:
        - "chrom" (str): chromosome name. 
        - "start" (int): start genomic position of the feature, 1-based and
          inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
    phased_snp_fn : str
        A TSV or VCF file listing phased SNPs.
        If TSV, it should be header-free, containing six columns:
        - "chrom" (str): chromosome name. 
        - "pos" (int): genomic position of the SNP, 1-based.
        - "ref" (str): the reference allele (REF) of the SNP.
        - "alt" (str): the alternative allele (ALT) of the SNP.
        - "ref_hap" (int): the haplotype index of the "ref", 0 or 1.
        - "alt_hap" (int): the haplotype index of the "alt", 0 or 1.
        If VCF, it should store the phased genotype in the "GT" within the
        "FORMAT" field.
    out_dir : str
        Output directory.
    debug_level : {0, 1, 2}
        The debugging level, the larger the number is, more detailed debugging
        information will be outputted.
    ncores : int, default 1
        Number of cores.
    cell_tag : str or None, default "CB"
        Tag for cell barcodes, set to None when using sample IDs.
    umi_tag : str or None, default "UB"
        Tag for UMI, set to None when reads only.
    min_count : int, default 20
        Minimum aggragated count for SNP.
    min_maf : float, default 0.1
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
    out_feature_dirs : list of str or None, default None
        A list of output folders for feature-specific results.
        If None, subfolders will be created under the `out_dir/alignments`.

    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    conf = Config()
    #init_logging(stream = sys.stdout)

    conf.sam_fn = sam_fn
    conf.barcode_fn = barcode_fn
    conf.feature_fn = feature_fn
    conf.snp_fn = phased_snp_fn
    conf.out_dir = out_dir
    conf.debug = debug_level

    conf.cell_tag = cell_tag
    conf.umi_tag = umi_tag
    conf.ncores = ncores

    conf.min_count = min_count
    conf.min_maf = min_maf

    conf.strandness = strandness
    conf.min_include = min_include
    conf.multi_mapper_how = multi_mapper_how
    conf.xf_tag = xf_tag
    conf.gene_tag = gene_tag
    
    conf.min_mapq = min_mapq
    conf.min_len = min_len
    conf.incl_flag = incl_flag
    conf.excl_flag = excl_flag
    conf.no_orphan = no_orphan
    
    conf.out_feature_dirs = out_feature_dirs
    conf.no_orphan_post_qc = no_orphan

    ret, res = afc_run(conf)
    return((ret, res))



def afc_core(conf):
    info("preprocessing ...")
    data = afc_pp(conf)
    
    reg_list = data["reg_list"]
    snp_set = data["snp_set"]
    samples = data["samples"]
    
    n_samples = len(samples)
    
    count_dir = os.path.join(conf.out_dir, "matrix")
    os.makedirs(count_dir, exist_ok = True)
    
    info("save feature annotations ...")
    out_feature_fn = os.path.join(
        count_dir, conf.out_prefix + ".features.tsv")
    with open(out_feature_fn, "w") as fp:
        for reg in reg_list:
            fp.write("%s\t%d\t%d\t%s\t%s\n" % \
                    (reg.chrom, reg.start, reg.end - 1, reg.name, reg.strand))
            
    info("save cell IDs ...")
    out_sample_fn = os.path.join(
        count_dir, conf.out_prefix + ".samples.tsv")
    with open(out_sample_fn, "w") as fp:
        fp.write("".join([smp + "\n" for smp in samples]))
    

    # extract SNPs for each feature.
    n = 0
    for reg in reg_list:
        snp_list = snp_set.fetch(reg.chrom, reg.start, reg.end)
        if snp_list and len(snp_list) > 0:
            reg.snp_list = snp_list
            n += 1
        else:
            reg.snp_list = []
            if conf.debug > 0:
                debug("no SNP fetched for feature '%s'." % reg.name)
    info("%d features contain SNPs." % n)


    # assign features to several batches of result folders, to avoid exceeding
    # the maximum number of files/sub-folders in one folder.
    assert len(reg_list) == len(conf.out_feature_dirs)
    for reg, feature_dir in zip(reg_list, conf.out_feature_dirs):
        reg.res_dir = feature_dir
        reg.init_allele_data(alleles = conf.cumi_alleles)
    conf.out_feature_dirs.clear()
    conf.out_feature_dirs = None
        

    tmp_dir = os.path.join(conf.out_dir, "tmp_afc")
    os.makedirs(tmp_dir, exist_ok = True)
    

    # split feature list and save to file.
    info("split feature list and save to file ...")

    fet_obj_fn = os.path.join(
        conf.out_dir, conf.out_prefix + ".features.pickle")
    save_feature_objects(reg_list, fet_obj_fn)
    
    # Note, here
    # - max_n_batch: to account for the max allowed files and subfolders in
    #   one folder.
    #   Currently, 6 files output in each batch.
    # - batch_per_core: it seems the overall running time becomes longer when
    #   the number of batches increases too much, possibly due to the overhead
    #   of loading BAM and features in each batch.
    m_reg = len(reg_list)
    bd_m, bd_n, bd_reg_indices = split_n2batch(
            m_reg, conf.ncores, batch_per_core = 3, max_n_batch = 5000)
    info("features are split into %d batches." % bd_m)
    
    bd_dir_list = []
    for idx in range(bd_m):
        d = os.path.join(tmp_dir, "%d" % idx)
        os.makedirs(d, exist_ok = True)
        bd_dir_list.append(d)


    reg_fn_list = []
    for idx, (b, e) in enumerate(bd_reg_indices):
        fn = os.path.join(bd_dir_list[idx], "fet.b%d.pickle" % idx)
        save_feature_objects(reg_list[b:e], fn)
        reg_fn_list.append(fn)
    
        
    # prepare args for multiprocessing.
    out_mtx_fns = {}
    for ale in conf.alleles:
        out_mtx_fns[ale] = os.path.join(
            count_dir, conf.out_prefix + ".%s.mtx" % ale)

    args_fn_list = []
    for idx in range(bd_m):
        mtx_fns = {ale: os.path.join(bd_dir_list[idx], "%s.b%d.mtx" % \
                    (ale, idx)) for ale in conf.alleles}
        args = dict(
            reg_obj_fn = reg_fn_list[idx],
            sam_fn = conf.sam_fn,
            out_mtx_fns = mtx_fns,
            samples = samples,
            batch_idx = idx,
            conf = conf
        )
        fn = os.path.join(bd_dir_list[idx], "args.b%d.pickle" % idx)
        save_pickle(args, fn)
        args_fn_list.append(fn)
        del args
        
    for reg in reg_list:  # save memory
        del reg
    snp_set.destroy()
    del reg_list
    del snp_set
    del samples
    del data
    gc.collect()
    reg_list = snp_set = samples = data = None


    # allele-specific counting with multi-processing.
    info("allele-specific counting with %d cores ..." % min(conf.ncores, bd_m))

    pool = multiprocessing.Pool(processes = min(conf.ncores, bd_m))
    mp_result = []
    for i in range(bd_m):
        args = load_pickle(args_fn_list[i])
        mp_result.append(pool.apply_async(
            func = fc_features, 
            kwds = args,
            callback = show_progress,
            error_callback = mp_error_handler
        ))
        del args
        gc.collect()
    pool.close()
    pool.join()

    info("multiprocessing done!")

    mp_result = [res.get() for res in mp_result]
            

    # merge feature objects containing post-filtering SNPs.
    info("merge feature objects ...")

    fet_obj_snp_filter_fn = fet_obj_fn.replace(".pickle", ".snp_filter.pickle")
    reg_list = []
    for fn in reg_fn_list:
        lst = load_feature_objects(fn)
        reg_list.extend(lst)
    save_feature_objects(reg_list, fet_obj_snp_filter_fn)


    # merge count matrices.
    info("merge output count matrices ...")

    for ale in out_mtx_fns.keys():
        if merge_mtx(
            [res["out_mtx_fns"][ale] for res in mp_result], ZF_F_GZIP,
            out_mtx_fns[ale], "w", ZF_F_PLAIN,
            [res["nr_reg"] for res in mp_result], n_samples,
            sum([res["nr_mtx"][ale] for res in mp_result]),
            remove = True
        ) < 0:
            error("errcode -17")
            raise ValueError
    
    out_fet_obj_fn = fet_obj_snp_filter_fn
    shutil.rmtree(tmp_dir)
    
    
    # process multi-feature UMIs.
    if conf.multi_mapper_how == "discard":
        info("multi_mapper_how = '%s'; processing multi-feature UMIs ..." % \
            conf.multi_mapper_how)
        
        count_dir = os.path.join(conf.out_dir, "matrix_uniq")
        os.makedirs(count_dir, exist_ok = True)
        
        tmp_dir = os.path.join(conf.out_dir, "tmp_mfu")
        os.makedirs(tmp_dir, exist_ok = True)
        
        fet_obj_mfu_fn = fet_obj_snp_filter_fn.replace(".pickle", ".mfu.pickle")
        
        res = mfu_main(
            alleles = conf.cumi_alleles,
            multi_mapper_how = conf.multi_mapper_how,
            fet_obj_fn = fet_obj_snp_filter_fn,
            out_fet_obj_fn = fet_obj_mfu_fn,
            sample_fn = out_sample_fn,
            feature_fn = out_feature_fn,
            count_dir = count_dir,
            tmp_dir = tmp_dir,
            out_prefix = "afc",
            ncores = conf.ncores
        )
        
        out_sample_fn = res["out_sample_fn"]
        out_feature_fn = res["out_feature_fn"]
        out_fet_obj_fn = res["out_fet_obj_fn"]
        out_mtx_fns = res["out_mtx_fns"]
        
        shutil.rmtree(tmp_dir)
        
    elif conf.multi_mapper_how == "duplicate":
        info("multi_mapper_how = '%s'; skip processing multi-feature UMIs ..." % \
             conf.multi_mapper_how)
        
    else:
        raise ValueError("invalid multi_mapper_how = '%s'." % conf.multi_mapper_how)


    # construct adata and save into h5ad file.
    info("construct adata and save into h5ad file ...")
    
    out_count_fn = os.path.join(
        conf.out_dir, conf.out_prefix + ".counts.h5ad")

    adata = None
    for idx, ale in enumerate(out_mtx_fns.keys()):
        dat = load_adata(
            mtx_fn = out_mtx_fns[ale],
            cell_fn = out_sample_fn,
            feature_fn = out_feature_fn,
            cell_columns = ["cell"],
            feature_columns = ["chrom", "start", "end", "feature", "strand"],
            row_is_cell = False,
            sparse_type = "csr"
        )
        if idx == 0:
            adata = dat
            adata.layers[ale] = dat.X
            adata.X = None
        else:
            adata.layers[ale] = dat.X

    # TODO: when adata .X is None, saving the adata into “.h5py” file can raise error 
    # `while writing key 'x' of <class 'h5py._hl.group.group'> to / `
    # for specific version of anndata and h5py.
    save_h5ad(adata.transpose(), out_count_fn)


    # clean
    info("clean ...")


    res = dict(
        # fet_obj_fn : str
        #   Path to a python pickle file storing the `reg_list`.
        #   It will be re-loaded for read sampling.
        fet_obj_fn = out_fet_obj_fn,

        # count_fn : str
        #   Path to a ".adata" file storing a :class:`~anndata.Anndata`
        #   object, which contains all allele-specific *cell x feature* count
        #   matrices.
        count_fn = out_count_fn
    )
    return(res)



def afc_run(conf):
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    try:
        res = afc_core(conf)
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



def afc_pp(conf):
    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")
    
    assert_e(conf.sam_fn)
    
    assert_e(conf.barcode_fn)
    barcodes = sorted(load_barcodes(conf.barcode_fn))
    assert len(set(barcodes)) == len(barcodes)
        
    samples = barcodes
    info("load %d cells." % len(samples))

    
    if not conf.out_dir:
        raise ValueError("out dir needed!")
    os.makedirs(conf.out_dir, exist_ok = True)

    
    assert_e(conf.feature_fn)
    reg_list = load_feature_from_txt(conf.feature_fn)
    if not reg_list:
        error("failed to load feature file.")
        raise ValueError
    info("load %d features." % len(reg_list))


    assert_e(conf.snp_fn)
    snp_set = None
    fn = conf.snp_fn.lower()
    if fn.endswith(".vcf") or fn.endswith(".vcf.gz") or \
                fn.endswith(".vcf.bgz"):
        snp_set = load_snp_from_vcf(conf.snp_fn)
    else:
        snp_set = load_snp_from_tsv(conf.snp_fn)
    if not snp_set or snp_set.get_n() <= 0:
        raise ValueError
    info("load %d SNPs." % snp_set.get_n())


    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None
    if conf.cell_tag and barcodes:
        pass       
    elif (not conf.cell_tag) ^ (not barcodes):
        raise ValueError("should not specify cell_tag or barcodes alone.")
    else:
        pass    

    if conf.umi_tag:
        if conf.umi_tag.upper() == "AUTO":
            if barcodes is None:
                conf.umi_tag = None
            else:
                conf.umi_tag = conf.defaults.UMI_TAG_BC
        elif conf.umi_tag.upper() == "NONE":
            conf.umi_tag = None
    else:
        pass

    
    assert conf.strandness in ("forward", "reverse", "unstranded")


    if conf.excl_flag < 0:
        if conf.use_umi():
            conf.excl_flag = conf.defaults.EXCL_FLAG_UMI
        else:
            conf.excl_flag = conf.defaults.EXCL_FLAG_XUMI
            
    
    if conf.out_feature_dirs is None:
        feature_dir = os.path.join(conf.out_dir, "features")
        os.makedirs(feature_dir, exist_ok = True)
        conf.out_feature_dirs = assign_feature_batch(
            feature_names = [reg.name for reg in reg_list],
            root_dir = feature_dir,
            batch_size = 1000
        )
    else:
        for fet_dir in conf.out_feature_dirs:
            assert_e(fet_dir)
            
            
    conf.pe_mode = detect_pe_mode(conf.sam_fn)
            
            
    info("updated configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")


    data = dict(
        # barcodes : list of str or None
        #   A list of cell barcodes.
        #   None if sample IDs are used.
        barcodes = barcodes,

        # samples : list of str
        #   A list of cell barcodes (droplet-based data) or sample IDs (
        #   well-based data).
        #   It will be used as output IDs of each cell.
        samples = samples,

        # reg_list : list of utils.gfeature.Feature
        #   A list of features.
        reg_list = reg_list,
        
        # snp_set : utils.gfeature.SNPSet
        #   The object storing a set of SNPs.
        snp_set = snp_set
    )
        
    return(data)



def show_progress(rv = None):
    return(rv)
