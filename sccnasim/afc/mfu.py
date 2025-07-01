# mfu.py - multi-feature UMIs.

# To process the multi-feature UMIs given the feature-specific CUMIs.
# While it is easy to load all CUMIs into an DataFrame and then group by
# cells and UMIs, it may be memory consuming when the total number of CUMIs
# is huge.
# Therefore, we use the strategy of divide-and-conquer to split the CUMIs
# into batches and merge them after processing CUMIs in each batch.



import gc
import multiprocessing
import numpy as np
import os
import pandas as pd
import shutil
from logging import info
from ..utils.cumi import load_cumi
from ..utils.gfeature import load_feature_objects, save_feature_objects, load_features
from ..utils.io import load_samples
from ..xlib.xbase import is_file_empty
from ..xlib.xfile import zopen, zconcat, ZF_F_GZIP, ZF_F_PLAIN
from ..xlib.xio import load_pickle, save_pickle
from ..xlib.xthread import split_n2batch, mp_error_handler



def mfu_main(
    alleles,
    multi_mapper_how,
    fet_obj_fn,
    out_fet_obj_fn,
    sample_fn,
    feature_fn,
    count_dir,
    tmp_dir,
    out_prefix,
    ncores
):
    """Main function for processing multi-feature UMIs.
    
    This function aims to process multi-feature UMIs, specifically,
    (1) find the multi-feature UMIs from combined allele-specific UMIs of all
        features and process them according to `multi_mapper_how`.
    (2) re-calculate the allele-specific cell x gene count matrices, and
        update the corresponding file paths in each feature object.
    
    Parameters
    ----------
    alleles : list of str
        A list of alleles.
    multi_mapper_how : {"discard", "duplicate"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "duplicate": count the UMI for every mapped gene.
    fet_obj_fn : str
        The python pickle file storing `~..utils.gfeature.Feature` objects.
    out_fet_obj_fn : str
        The python pickle file storing updated `~..utils.gfeature.Feature` 
        objects.
    sample_fn : str
        File storing cell IDs or barcodes.
        Its order should match that in count matrices.
    feature_fn : str
        File storing feature annotations.
    count_dir : str
        The output folder to store the count matrices.
    tmp_dir : str
        Folder to store temporary data.
    out_prefix : str
        The prefix to output files.
    ncores : int
        Number of cores.

    Returns
    -------
    dict
        Results.
    """
    # check args.
    reg_list = load_feature_objects(fet_obj_fn)
    cells = load_samples(sample_fn)         # ndarray
    features = load_features(feature_fn)      # DataFrame
    features = features["feature"].to_numpy()
    
    assert len(reg_list) == features.shape[0]
    for i, reg in enumerate(reg_list):
        assert reg.name == features[i]
    
    os.makedirs(count_dir, exist_ok = True)
    os.makedirs(tmp_dir, exist_ok = True)
    
    del reg_list
    gc.collect()
    
    step = 1
    
    
    # copy files.
    out_sample_fn = os.path.join(count_dir, "%s.samples.tsv" % out_prefix)
    shutil.copy(sample_fn, out_sample_fn)
    out_feature_fn = os.path.join(count_dir, "%s.features.tsv" % out_prefix)
    shutil.copy(feature_fn, out_feature_fn)
    
    
    # merge allele-specific CUMI files of all features.
    info("merge allele-specific CUMI files of all features ...")
    
    cumi_merge_fn = os.path.join(count_dir, "%s.cumi.merged.tsv" % out_prefix)
    res_dir = os.path.join(tmp_dir, "%d_cumi_merge" % step)
    os.makedirs(res_dir, exist_ok = True)
    merge_cumis(
        alleles = alleles,
        fet_obj_fn = fet_obj_fn,
        out_fn = cumi_merge_fn,
        tmp_dir = res_dir,
        ncores = ncores
    )
    step += 1
    
    
    # process cell-specific multi-feature CUMIs.
    info("process cell-specific multi-feature CUMIs ...")
    
    cumi_cs_fn = os.path.join(count_dir, "%s.cumi.merged.cs.tsv" % out_prefix)
    res_dir = os.path.join(tmp_dir, "%d_cumi_cs" % step)
    os.makedirs(res_dir, exist_ok = True)
    mfu_cs_main(
        in_fn = cumi_merge_fn,
        out_fn = cumi_cs_fn,
        cells = cells,
        multi_mapper_how = multi_mapper_how,
        tmp_dir = res_dir,
        ncores = ncores
    )
    step += 1
    
    
    # extract allele-specific CUMIs.
    info("extract allele-specific CUMIs ...")
    
    res_dir = os.path.join(tmp_dir, "%d_cumi_as" % step)
    os.makedirs(res_dir, exist_ok = True)
    cumi_as_fns = [os.path.join(res_dir, "%s.cumi.tsv" % ale) \
                     for ale in alleles]
    mfu_as_main(
        in_fn = cumi_cs_fn,
        out_fns = cumi_as_fns,
        alleles = alleles
    )
    step += 1
    
    
    # extract feature-specific CUMIs and do counting.
    info("extract feature-specific CUMIs and do counting ...")
    
    res_dir = os.path.join(tmp_dir, "%d_cumi_fs" % step)
    os.makedirs(res_dir, exist_ok = True)
    
    cumi_obj_fns = []
    out_mtx_fns = {}
    
    for ale, cumi_fn in zip(alleles, cumi_as_fns):
        reg_list = load_feature_objects(fet_obj_fn)
        out_fns = [reg.allele_data[ale].seed_cumi_fn.replace(
            ".cumi.tsv", ".mfu.cumi.tsv") for reg in reg_list]
        del reg_list
        gc.collect()
        
        mtx_fn = os.path.join(count_dir, "%s.%s.mtx" % (out_prefix, ale))
        out_mtx_fns[ale] = mtx_fn
        
        n_rec_mtx = mfu_fs_main(
            in_fn = cumi_fn,
            out_files = out_fns,
            matrix_fn = mtx_fn,
            allele = ale,
            cells = cells,
            features = features,
            tmp_dir = os.path.join(res_dir, ale),
            ncores = ncores
        )
        
        fn = os.path.join(res_dir, "%s.out.cumi_fn.list.pickle" % ale)
        save_pickle(out_fns, fn)
        cumi_obj_fns.append(fn)
        
        info("%d matrix records extracted for allele '%s'." % (n_rec_mtx, ale))
    step += 1
    
    
    # update CUMI files in the feature objects.
    info("update CUMI files in the feature objects ...")
    
    reg_list = load_feature_objects(fet_obj_fn)
    for ale, obj_fn in zip(alleles, cumi_obj_fns):
        fn_list = load_pickle(obj_fn)
        for reg, fn in zip(reg_list, fn_list):
            reg.allele_data[ale].seed_cumi_fn = fn
    save_feature_objects(reg_list, out_fet_obj_fn)
    
    
    res = dict(
        # out_sample_fn : str
        #   File storing cells. It is in `count_dir`.
        out_sample_fn = out_sample_fn,
        
        # out_feature_fn : str
        #   File storing features. It is in `count_dir`.        
        out_feature_fn = out_feature_fn,
        
        # out_fet_obj_fn : str
        #   The python pickle file storing updated feature objects.
        out_fet_obj_fn = out_fet_obj_fn,
        
        # out_mtx_fns : dict of {str : str}
        #   Key is allele, value is the allele-specific sparse count 
        #   matrix file.
        out_mtx_fns = out_mtx_fns
    )
    return(res)

    
    
def mfu_cs_main(
    in_fn,
    out_fn,
    cells,
    multi_mapper_how,
    tmp_dir,
    ncores
):
    """Main function of processing cell-specific multi-feature CUMIs.
    
    Parameters
    ----------
    in_fn : str
        Path to the input 4-column CUMI file.
    out_fn : str
        Path to the output 4-column CUMI file.
    cells : list of str
        A list of cell barcodes.
        Its order should match that in count matrices.
    multi_mapper_how : {"discard", "duplicate"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "duplicate": count the UMI for every mapped gene.
    tmp_dir : str
        Folder to store temporary data.
    ncores : int
        Number of cores.

    Returns
    -------
    Void.
    """
    # check args.
    assert multi_mapper_how in ("discard", "duplicate")
    os.makedirs(tmp_dir, exist_ok = True)

    __mfu_cs_batch(
        in_fn = in_fn,
        out_fn = out_fn,
        cells = cells,
        multi_mapper_how = multi_mapper_how,
        tmp_dir = tmp_dir,
        ncores = ncores,
        max_per_batch = 500,
        depth = 0
    )
    

    
def __mfu_cs_batch(
    in_fn,
    out_fn,
    cells,
    multi_mapper_how,
    tmp_dir,
    ncores,
    max_per_batch,
    depth
):
    """Recursive function for `mfu_cs_main()`.
    
    To avoid too many cells (and hence CUMIs) in one batch, this function
    recursively splits large combined file into smaller batches, until the 
    batch size is small than given `max_per_batch`.
    
    Parameters
    ----------
    in_fn : str
        Path to the input 4-column CUMI file.
    out_fn : str
        Path to the output 4-column CUMI file.
    cells : list of str
        A list of cell barcodes.
    multi_mapper_how : {"discard", "duplicate"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "duplicate": count the UMI for every mapped gene.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int, default 1
        Number of cores.
    max_per_batch : int
        Maximum number of `cells` allowed to be processed simultaneously.
    depth : int
        Depth index, 0-based.

    Returns
    -------
    Void.
    """
    n = len(cells)   
    
    if n <= max_per_batch:
        mfu_cs(
            in_fn = in_fn,
            out_fn = out_fn,
            multi_mapper_how = multi_mapper_how
        )
        return
    
    os.makedirs(tmp_dir, exist_ok = True)


    # split the input CUMI file into smaller batches.
    # Note, here
    # - max_n_batch: to account for the issue of "max open files" when
    #   splitting the large combined file into smaller batches.
    #   It will open every batch-specific splitted file simultaneously, 
    #   in total `n_batch` files.
    bd_m, bd_n, bd_indices = split_n2batch(
        n, ncores, min_n_batch = 10, max_n_batch = 300)
    
    fp_list = []
    idx_map = {}
    batches = []
    for idx, (b, e) in enumerate(bd_indices):
        bd_in_fn = os.path.join(tmp_dir, "%d_%d.in.cumi.tsv" % (depth, idx))
        bd_out_fn = os.path.join(tmp_dir, "%d_%d.out.cumi.tsv" % (depth, idx))
        fp = zopen(bd_in_fn, "w", ZF_F_PLAIN)
        fp_list.append(fp)
        for i in range(b, e):
            assert cells[i] not in idx_map
            idx_map[cells[i]] = fp
        batches.append((b, e, bd_in_fn, bd_out_fn))
    
    in_fp = open(in_fn, "r")
    for line in in_fp:
        cell, _, s = line.partition("\t")
        assert cell in idx_map
        fp = idx_map[cell]
        fp.write(line)
    in_fp.close()
    
    for fp in fp_list:
        fp.close()
        

    # next round of extracting and splitting.
    if ncores <= 1:
        for idx, (b, e, bd_in_fn, bd_out_fn) in enumerate(batches):
            res_dir = os.path.join(tmp_dir, "%d_%d" % (depth, idx))
            os.makedirs(res_dir, exist_ok = True)
            __mfu_cs_batch(
                in_fn = bd_in_fn,
                out_fn = bd_out_fn,
                cells = cells[b:e],
                multi_mapper_how = multi_mapper_how,
                tmp_dir = res_dir,
                ncores = 1,
                max_per_batch = max_per_batch,
                depth = depth + 1
            )
    else:
        mp_res = []
        pool = multiprocessing.Pool(processes = min(ncores, bd_m))
        for idx, (b, e, bd_in_fn, bd_out_fn) in enumerate(batches):
            res_dir = os.path.join(tmp_dir, "%d_%d" % (depth, idx))
            os.makedirs(res_dir, exist_ok = True)
            mp_res.append(pool.apply_async(
                func = __mfu_cs_batch,
                kwds = dict(
                    in_fn = bd_in_fn,
                    out_fn = bd_out_fn,
                    cells = cells[b:e],
                    multi_mapper_how = multi_mapper_how,
                    tmp_dir = res_dir,
                    ncores = 1,
                    max_per_batch = max_per_batch,
                    depth = depth + 1
                ),
                callback = None,
                error_callback = mp_error_handler
            ))
        pool.close()
        pool.join()
        
        
    # merge batch-specific CUMI files.
    zconcat(
        in_fn_list = [x[3] for x in batches],
        in_format = ZF_F_PLAIN,
        out_fn = out_fn, 
        out_fmode = "w", 
        out_format = ZF_F_PLAIN, 
        remove = False
    )

    del fp_list
    del idx_map
    del batches
    gc.collect()



def mfu_cs(
    in_fn,
    out_fn,
    multi_mapper_how
):
    """Process cell-specific multi-feature CUMIs for one batch of cells.
    
    Parameters
    ----------
    in_fn : str
        Path to the input 4-column CUMI file.
    out_fn : str
        Path to the output 4-column CUMI file.
    multi_mapper_how : {"discard", "duplicate"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "duplicate": count the UMI for every mapped gene.
        
    Returns
    -------
    Void.
    """
    df = mfu_load_cumi(in_fn)
    if multi_mapper_how == "discard":
        df = df.drop_duplicates(["cell", "umi"], keep = False, 
                                ignore_index = True)
    else:
        pass
    mfu_save_cumi(df, out_fn)



def mfu_as_main(in_fn, out_fns, alleles):
    """Extract allele-specific CUMIs.
    
    Parameters
    ----------
    in_fn : str
        Path to the input 4-column CUMI file.
    out_fns : list of str
        Path to the allele-specific 3-column CUMI files.
    alleles : list of str
        Alleles.
        
    Returns
    -------
    Void.    
    """
    assert len(out_fns) == len(alleles)
    
    idx_map = {}
    fp_list = []
    for ale, fn in zip(alleles, out_fns):
        fp = zopen(fn, "w", ZF_F_PLAIN)
        fp_list.append(fp)
        idx_map[ale] = fp
    
    in_fp = open(in_fn, "r")
    for line in in_fp:
        s, _, ale = line.rstrip().rpartition("\t")
        assert ale in idx_map
        fp = idx_map[ale]
        fp.write(s + "\n")
    in_fp.close()
    
    for fp in fp_list:
        fp.close()
        
        
        
def mfu_fs_main(
    in_fn,
    out_files,
    matrix_fn,
    allele,
    cells,
    features,
    tmp_dir,
    ncores
):
    """Main function of extracting feature-specific CUMIs and do counting for
    specific allele.
    
    Parameters
    ----------
    in_fn : str
        Path to the input 3-column CUMI file.
    out_files : list of str
        A list of output files storing feature-specific CUMIs from all cells.
    matrix_fn : str
        Path to the output sparse matrix file.
    allele : str
        Allele.
    cells : list of str
        All available cell barcodes.
    features : list of str
        All available feature names.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int
        Number of cores.

    Returns
    -------
    int
        Number of records in the `matrix_fn`.
    """
    os.makedirs(tmp_dir, exist_ok = True)
    
    tmp_mtx_fn = os.path.join(tmp_dir, "tmp.matrix.mtx")
    n_rec_mtx = __mfu_fs_batch(
        in_fn = in_fn,
        b0 = 0,
        e0 = len(out_files) - 1,
        out_files = out_files,
        matrix_fn = tmp_mtx_fn,
        allele = allele,
        cells = cells,
        features = features,
        tmp_dir = tmp_dir,
        ncores = ncores,
        matrix_sep = "\t",
        max_per_batch = 300,
        depth = 0
    )
    
    # add header line into the matrix file.
    tmp_header_fn = os.path.join(tmp_dir, "tmp.matrix.header.tsv")
    s  = "%%MatrixMarket matrix coordinate integer general\n"
    s += "%%\n"
    s += "%d %d %d\n" % (len(features), len(cells), n_rec_mtx)
    with open(tmp_header_fn, "w") as fp:
        fp.write(s)
        
    zconcat(
        in_fn_list = [tmp_header_fn, tmp_mtx_fn],
        in_format = ZF_F_PLAIN, 
        out_fn = matrix_fn, 
        out_fmode = "w", 
        out_format = ZF_F_PLAIN, 
        remove = False
    )
    return(n_rec_mtx)



def __mfu_fs_batch(
    in_fn,
    b0,
    e0,
    out_files,
    matrix_fn,
    allele,
    cells,
    features,
    tmp_dir,
    ncores,
    matrix_sep,
    max_per_batch,
    depth
):
    """Recursive function for `mfu_fs_main()`.
    
    To avoid the issue of `max open files` in mfu_fs_main(), this function
    recursively splits large combined file into smaller batches, until the 
    batch size is small than given `max_per_batch`.
    
    Parameters
    ----------
    in_fn   
    out_files
    matrix_fn
    allele
    cells
    features
    tmp_dir
    ncores
        See :func:`mfu_fs_main()`.
    matrix_sep : str
        Delimiter in the `matrix_fn`.        
    b0 : int
        The transcriptomics-scale index of the first feature in this batch.
        0-based, inclusive.
    e0 : int
        The transcriptomics-scale index of the last feature in this batch.
        0-based, inclusive. 
    max_per_batch : int
        Maximum number of features allowed to be processed simultaneously.
    depth : int
        Depth index, 0-based.

    Returns
    -------
    int
        Number of records in the `matrix_fn`.
    """
    p = len(out_files)
    assert p == e0 - b0 + 1
    
    n_rec_mtx = 0
    
    if p <= max_per_batch:
        n_rec_mtx = mfu_fs(
            in_fn = in_fn,
            out_files = out_files,
            matrix_fn = matrix_fn,
            matrix_sep = matrix_sep,
            b = b0,
            e = e0,
            cells = cells,
            features = features
        )
        return(n_rec_mtx)

    os.makedirs(tmp_dir, exist_ok = True)

    # split the input CUMI file into smaller batches.
    # Note, here
    # - max_n_batch: to account for the issue of "max open files" when
    #   splitting the large combined file into smaller batches.
    #   It will open every batch-specific splitted file simultaneously, 
    #   in total `n_batch` files.
    bd_m, bd_n, bd_indices = split_n2batch(
        p, ncores, min_n_batch = 30, max_n_batch = 300)
    
    fp_list = []
    idx_map = {}
    batches = []
    for idx, (b, e) in enumerate(bd_indices):
        cumi_fn = os.path.join(tmp_dir, "%d_%d.in.cumi.tsv" % (depth, idx))
        fp = zopen(cumi_fn, "w", ZF_F_PLAIN)
        fp_list.append(fp)
        
        b += b0
        e += b0
        for i in range(b, e):
            fet = features[i]
            assert fet not in idx_map
            idx_map[fet] = fp
            
        mtx_fn = os.path.join(tmp_dir, "%d_%d.out.matrix.mtx" % (depth, idx))
        batches.append((b, e - 1, cumi_fn, mtx_fn))
    
    in_fp = open(in_fn, "r")
    for line in in_fp:
        s, _, fet = line.rstrip().rpartition("\t")
        assert fet in idx_map
        fp = idx_map[fet]
        fp.write(line)
    in_fp.close()
    
    for fp in fp_list:
        fp.close()
        

    # next round of extracting and splitting.
    if ncores <= 1:
        for idx, (b, e, cumi_fn, mtx_fn) in enumerate(batches):
            res_dir = os.path.join(tmp_dir, "%d_%d" % (depth, idx))
            os.makedirs(res_dir, exist_ok = True)
            n_rec = __mfu_fs_batch(
                in_fn = cumi_fn,
                b0 = b,
                e0 = e,
                out_files = out_files[(b-b0):(e+1-b0)],
                matrix_fn = mtx_fn,
                allele = allele,
                cells = cells,
                features = features,
                tmp_dir = res_dir,
                ncores = 1,
                matrix_sep = matrix_sep,
                max_per_batch = max_per_batch,
                depth = depth + 1
            )
            n_rec_mtx += n_rec
    else:
        mp_res = []
        pool = multiprocessing.Pool(processes = min(ncores, bd_m))
        for idx, (b, e, cumi_fn, mtx_fn) in enumerate(batches):
            res_dir = os.path.join(tmp_dir, "%d_%d" % (depth, idx))
            os.makedirs(res_dir, exist_ok = True)
            mp_res.append(pool.apply_async(
                func = __mfu_fs_batch,
                kwds = dict(
                    in_fn = cumi_fn,
                    b0 = b,
                    e0 = e,
                    out_files = out_files[(b-b0):(e+1-b0)],
                    matrix_fn = mtx_fn,
                    allele = allele,
                    cells = cells,
                    features = features,
                    tmp_dir = res_dir,
                    ncores = 1,
                    matrix_sep = matrix_sep,
                    max_per_batch = max_per_batch,
                    depth = depth + 1
                ),
                callback = None,
                error_callback = mp_error_handler
            ))
        pool.close()
        pool.join()
        
        n_rec_mtx = np.sum([res.get() for res in mp_res])
        
        
    # merge batch-specific matrix files.
    zconcat(
        in_fn_list = [x[3] for x in batches],
        in_format = ZF_F_PLAIN, 
        out_fn = matrix_fn, 
        out_fmode = "w", 
        out_format = ZF_F_PLAIN, 
        remove = False
    )
        
        
    del fp_list
    del idx_map
    del batches
    gc.collect()
    
    return(n_rec_mtx)



def mfu_fs(
    in_fn,
    out_files,
    matrix_fn,
    matrix_sep,
    b,
    e,
    cells,
    features
):
    """Extract feature-specific CUMIs from combined file and do counting.
    
    Parameters
    ----------
    in_fn : str
        Path to the file storing CUMIs of combined features.
    out_files : list of str
        A list of feature-specific CUMI files.
        Its length should be equal to `e-b+1`.
    matrix_fn : str
        Path to the output sparse matrix file.
    matrix_sep : str
        Delimiter in the `matrix_fn`.
    b : int
        The transcriptomics-scale index of the first feature in this batch.
        0-based, inclusive.
    e : int
        The transcriptomics-scale index of the last feature in this batch.
        0-based, inclusive.
    cells : list of str
        All available cell barcodes.
    features : list of str
        All available feature names.
        
    Returns
    -------
    int
        Number of records in the `matrix_fn`.
    """
    def __save_cumi(d, fn_map):
        fet = d["feature"].values[0]
        fn = fn_map[fet]
        d[["cell", "umi"]].to_csv(
            fn, sep = "\t", header = False, index = False)


    # check args.
    assert len(out_files) == e - b + 1
    
    df = None
    if is_file_empty(in_fn):
        df = pd.DataFrame(columns = ["cell", "umi", "feature"],
                        dtype = "string")
    else:
        df = pd.read_table(in_fn, header = None)
        df.columns = ["cell", "umi", "feature"]
    
    
    # create output CUMI file for every features in this batch.
    for fn in out_files:
        with open(fn, "w") as fp:
            pass
    
    # save CUMIs file for features having CUMI records.
    fn_map = {}
    for fet, fn in zip(features[b:(e+1)], out_files):
        fn_map[fet] = fn

    # use group_keys=False to avoid the group key being set as index.
    df.groupby("feature").apply(__save_cumi, fn_map = fn_map)
    
    
    # count CUMIs.
    df["cell_idx"] = df["cell"].map({c:i+1 for i, c in enumerate(cells)})
    df["feature_idx"] = df["feature"].map({f:i+1 for i, f in enumerate(features)})
    stat = df[["feature_idx", "cell_idx", "umi"]].groupby(
        ["feature_idx", "cell_idx"], as_index = False).count()
    stat = stat[["feature_idx", "cell_idx", "umi"]].sort_values(     # here column "umi" contains <n_umi>.
        by = ["feature_idx", "cell_idx"],
        ignore_index = False)
    stat.to_csv(matrix_fn, sep = matrix_sep, header = False, index = False)
    
    return(stat.shape[0])
    


def merge_cumis(
    alleles,
    fet_obj_fn,
    out_fn,
    tmp_dir,
    ncores
):
    """Merge all allele-specific CUMI files of all features.
    
    Parameters
    ----------
    alleles : list of str
        A list of alleles.
    fet_obj_fn : str
        The python pickle file storing `~..utils.gfeature.Feature` objects.
    out_fn : str
        Output file storing all CUMIs of all features.
    tmp_dir : str
        Folder to store temporary data.
    ncores : int
        Number of cores.

    Returns
    -------
    Void.
    """
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # split features into batches.
    reg_list = load_feature_objects(fet_obj_fn)
    p = len(reg_list)
    bd_m, bd_n, bd_indices = split_n2batch(
        p, ncores, min_per_batch = 200, max_per_batch = 500)
    
    batches = []
    for idx, (b, e) in enumerate(bd_indices):
        reg_fn = os.path.join(tmp_dir, "%d.features.pickle" % idx)
        save_pickle(reg_list[b:e], reg_fn)
        cumi_fn = os.path.join(tmp_dir, "%d.cumis.tsv" % idx)
        batches.append((reg_fn, cumi_fn))
    del reg_list
    gc.collect()
    
    
    # merge CUMIs in each batch.
    if ncores <= 1:
        for idx, (reg_fn, cumi_fn) in enumerate(batches):
            merge_cumis_batch(
                fet_obj_fn = reg_fn,
                alleles = alleles,
                out_fn = cumi_fn
            )
    else:
        mp_res = []
        pool = multiprocessing.Pool(processes = min(ncores, bd_m))
        for idx, (reg_fn, cumi_fn) in enumerate(batches):
            mp_res.append(pool.apply_async(
                func = merge_cumis_batch,
                kwds = dict(
                    fet_obj_fn = reg_fn,
                    alleles = alleles,
                    out_fn = cumi_fn
                ),
                callback = None,
                error_callback = mp_error_handler
            ))
        pool.close()
        pool.join()
        
        
    # merge batch-specific CUMI files.
    zconcat(
        in_fn_list = [x[1] for x in batches],
        in_format = ZF_F_PLAIN, 
        out_fn = out_fn, 
        out_fmode = "w", 
        out_format = ZF_F_PLAIN, 
        remove = False
    )
    
    
    
def merge_cumis_batch(
    fet_obj_fn,
    alleles,
    out_fn
):
    """Merge CUMIs for a batch of features.
    
    Parameters
    ----------
    fet_obj_fn : str
        The python pickle file storing `~..utils.gfeature.Feature` objects.
    alleles : list of str
        A list of alleles.
    out_fn : str
        Path to output file.

    Returns
    -------
    Void.
    """
    reg_list = load_feature_objects(fet_obj_fn)
    dat = []
    for idx, reg in enumerate(reg_list):
        df = merge_cumis_feature(
            fn_list = [reg.allele_data[ale].seed_cumi_fn for ale in alleles],
            alleles = alleles,
            name = reg.name
        )
        dat.append(df)
    df = pd.concat(dat, ignore_index = True)
    mfu_save_cumi(df, out_fn)



def merge_cumis_feature(
    fn_list,
    alleles,
    name
):
    """Merge CUMI files for one feature.
    
    Parameters
    ----------
    fn_list : list of str
        A list of allele-specific CUMI files.
    alleles : list of str
        A list of alleles.
    name : str
        Feature name.

    Returns
    -------
    pandas.DataFrame
        The merged CUMI data, containing four columns:
        - "cell" (str): cell barcode;
        - "umi" (str): UMI barcode;
        - "feature" (str): feature name;
        - "allele" (str): allele
    """
    # check args.
    assert len(fn_list) == len(alleles)
    
    dat = []
    for ale, fn in zip(alleles, fn_list):
        df = load_cumi(fn)
        if df.shape[0] == 0:
            df = pd.DataFrame(columns = ["cell", "umi", "feature", "allele"],
                             dtype = "string")
        else:
            df["feature"] = name
            df["allele"] = ale
        dat.append(df)
    df = pd.concat(dat, ignore_index = True)
    return(df)



def mfu_load_cumi(fn):
    """Load 4-column CUMIs from file."""
    df = None
    if is_file_empty(fn):
        df = pd.DataFrame(columns = ["cell", "umi", "feature", "allele"],
                        dtype = "string")
    else:
        df = pd.read_table(fn, header = None)
        df.columns = ["cell", "umi", "feature", "allele"]
    return(df)


def mfu_save_cumi(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)
