# main.py - main() function for count simulation.


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
import random
import sys
import time
import warnings

from anndata import ImplicitModificationWarning
from logging import info, error
from .config import Config
from .core import cs_init, cs_pp, cs_cna_core
from .gcna import calc_cn_ratio
from .io import cs_save_adata2mtx, save_params
from ..xlib.xbarcode import rand_cell_barcodes
from ..xlib.xdata import load_h5ad, save_h5ad
from ..xlib.xio import load_pickle, save_pickle



def cs_core(conf):
    info("init ...")
    data = cs_init(conf)
    
    
    info("preprocessing ...")
    
    pp_dir = os.path.join(conf.out_dir, "pp")
    os.makedirs(pp_dir, exist_ok = True)
    data, pp_res = cs_pp(data, conf, out_dir = pp_dir)
    
    adata = data["adata"]
    clone_anno = data["clone_anno"]
    cna_profile = data["cna_profile"]
    
    params = dict(
        # clones : pandas.Series
        #   The ID of CNA clones.
        clones = clone_anno["clone"],

        # cell_types : pandas.Series
        #   The source cell types used by `clones`.
        cell_types = clone_anno["cell_type"],

        # n_cell_each : list of int
        #   Number of cells in each of `clones`.
        n_cell_each = pp_res["n_cell_each"],

        # cna_profile : pandas.DataFrame
        #   The clonal CNA profile.
        cna_profile = cna_profile,

        # cna_features : dict of {str : numpy.ndarray of int}
        #   The overlapping features of each CNA region.
        #   Keys are ID of CNA region, values are the (transcriptomics scale)
        #   indexes of their overlapping features.
        cna_features = pp_res["cna_fet"],

        # size_factors_type : str or None
        #   The type of size factors, e.g., "libsize".
        #   None means that size factors are not used.
        size_factors_type = conf.size_factor,

        # size_factors_train : numpy.ndarray of float
        #   The cell-wise size factors from trainning data.
        size_factors_train = pp_res["size_factors_train"],

        # size_factors_simu : list of float
        #   The cell-wise simulated size factors.
        size_factors_simu = pp_res["size_factors_simu"],

        # marginal : {"auto", "poi", "nb", "zinb"}
        #   Type of marginal distribution.
        marginal = conf.marginal,

        # kwargs_fit_rd : dict
        #   The additional kwargs passed to function 
        #   :func:`~marginal.fit_RD_wrapper` for fitting read depth.
        kwargs_fit_rd = conf.kwargs_fit_rd,
        
        # libsize_ratio : float
        #   Ratio of library size of simulated cells compared to seed cells.
        libsize_ratio = conf.libsize_ratio
    )
    
    fn = os.path.join(pp_dir, "pp.output.params.pickle")
    save_pickle(params, fn)
    del params

    
    info("process allele separately ...")
    
    # some variables used downstream.
    features = adata.var["feature"]
    n, p = adata.shape
    
    # prepare allele-specific count data
    info("prepare allele-specific count data ...")
    
    count_fn_list = []
    dir_list = []
    for idx, allele in enumerate(conf.alleles):
        dir_ale = os.path.join(conf.out_dir, allele)
        os.makedirs(dir_ale, exist_ok = True)
        dir_list.append(dir_ale)
        
        adata_ale = ad.AnnData(
            X = adata.layers[allele],
            obs = adata.obs,
            var = adata.var
        )
        
        fn = os.path.join(dir_ale, "%s.input.counts.h5ad" % allele)
        save_h5ad(adata_ale, fn)
        count_fn_list.append(fn)
    
    
    # prepare allele-specific args.
    info("prepare allele-specific args ...")
    
    args_fn_list = []
    for idx, allele in enumerate(conf.alleles):
        profile = cna_profile.copy()
        normal = None            # copy number in normal cells.
        if allele == "A":
            profile["cn"] = profile["cn_ale0"]
            normal = 1
        elif allele == "B":
            profile["cn"] = profile["cn_ale1"]
            normal = 1
        elif allele == "U":
            profile["cn"] = profile["cn_ale0"] + profile["cn_ale1"]
            normal = 2
        else:
            raise ValueError

        cn_ratio = calc_cn_ratio(
            cna_profile = profile,
            cna_features = pp_res["cna_fet"],
            p = p,
            normal = normal,
            loss_allele_freq = conf.loss_allele_freq
        )
        
        kwargs_fit_rd = conf.kwargs_fit_rd.copy()
        if "min_nonzero_num" in kwargs_fit_rd:
            kwargs_fit_rd["min_nonzero_num"] = kwargs_fit_rd["min_nonzero_num"][idx]
        
        dir_ale = dir_list[idx]
        args = dict(
            count_fn = count_fn_list[idx],
            out_dir = dir_ale,
            out_prefix = allele,
            clones = clone_anno["clone"],
            cell_types = clone_anno["cell_type"],
            n_cell_each = pp_res["n_cell_each"],
            cn_ratio = cn_ratio,
            size_factor = conf.size_factor,
            size_factors_train = pp_res["size_factors_train"],
            size_factors_simu = pp_res["size_factors_simu"],
            marginal = conf.marginal, 
            kwargs_fit_rd = kwargs_fit_rd,
            libsize_ratio = conf.libsize_ratio,
            ncores = conf.ncores, 
            verbose = conf.verbose
        )
        
        fn = os.path.join(dir_ale, "%s.input.args.pickle" % allele)
        save_pickle(args, fn)
        args_fn_list.append(fn)

    del adata
    del clone_anno
    del cna_profile
    del data
    del pp_res
    del count_fn_list
    gc.collect()
    adata = clone_anno = cna_profile = data = pp_res = count_fn_list = None
        

    # simulate CNAs for each allele.
    info("simulate CNAs for each allele ...")
    
    count_fn_list = []
    params_fn_list = []
    for idx, allele in enumerate(conf.alleles):
        info("start simulating counts for allele '%s'." % allele)
        
        args = load_pickle(args_fn_list[idx])
        adata_ale = cs_cna_core(**args)
        
        assert np.all(adata_ale.var["feature"] == features)
        
        dir_ale = dir_list[idx]
        
        fn = os.path.join(dir_ale, "%s.simu.output.counts.h5ad" % allele)
        save_h5ad(adata_ale, fn)
        count_fn_list.append(fn)
        
        del adata_ale
        gc.collect()
        

    # merge results.
    info("merge results ...")
    
    adata_simu = None
    for idx, allele in enumerate(conf.alleles):
        adata_ale = load_h5ad(count_fn_list[idx])
        if idx == 0:
            adata_simu = ad.AnnData(
                X = None,
                obs = adata_ale.obs,
                var = adata_ale.var
            )
        else:
            assert np.all(
                adata_ale.obs["cell_type"] == adata_simu.obs["cell_type"])
            assert np.all(
                adata_ale.var["feature"] == adata_simu.var["feature"])
        adata_simu.layers[allele] = adata_ale.X

        
    info("generate random cell barcodes for simulated adata ...")
    
    if conf.barcode_whitelist_fn is None:
        adata_simu.obs["cell"] = rand_cell_barcodes(
            m = 16,
            n = adata_simu.shape[0],
            suffix = "-1",
            sort = True
        )
    else:
        adata_simu.obs["cell"] = sample_barcodes(
            fn = conf.barcode_whitelist_fn,
            n = adata_simu.shape[0],
            suffix = "-1",
            sort = True
        )
    
    
    info("update .var of simulated adata ...")
    
    adata = load_h5ad(conf.count_fn)
    assert np.all(adata_simu.var["feature"] == adata.var["feature"])
    adata_simu.var = adata.var
    del adata


    info("save simulated counts ...")
    
    count_fn = os.path.join(conf.out_dir, "%s.counts.h5ad" % \
                                conf.out_prefix)
    save_h5ad(adata_simu, count_fn)
    

    cs_save_adata2mtx(
        adata = adata_simu,
        layers = conf.alleles,
        out_dir = os.path.join(conf.out_dir, "matrix"),
        row_is_cell = True,
        cell_columns = ["cell", "cell_type"],
        barcode_columns = ["cell"]
    )
    
    res = dict(
        count_fn = count_fn
    )
    return(res)



def cs_run(conf):
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            warnings.simplefilter("ignore", UserWarning)
            warnings.simplefilter("ignore", ImplicitModificationWarning)
            res = cs_core(conf)
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



def cs_wrapper(
    count_fn,
    clone_anno_fn, cna_profile_fn,
    out_dir,
    size_factor = "libsize", 
    marginal = "auto",
    libsize_ratio = 1.0,
    loss_allele_freq = 0.01,
    cna_mode = "hap-aware",
    barcode_whitelist_fn = None,
    ncores = 1, verbose = False,
    kwargs_fit_sf = None, kwargs_fit_rd = None
):
    """Wrapper for running the cs (count simulation) module.

    Parameters
    ----------
    count_fn : str
        A h5ad file storing the *cell x feature* count matrices for allele
        A, B, U in three layers "A", "B", "U", respectively.
        Its `.obs` should contain columns:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
        Its `.var` should contain columns:
        - "chrom" (str): chromosome name of the feature.
        - "start" (int): start genomic position of the feature, 1-based
          and inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
    clone_anno_fn : str
        A TSV file listing clonal annotation information.
        It is header-free and its first 3 columns are:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the source cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, 
          then it will be set as the number of cells in `source_cell_type`.
    cna_profile_fn : str
        A TSV file listing clonal CNA profiles.
        It is header-free and its first several columns are:
        - "chrom" (str): chromosome name of the CNA region.
        - "start" (int): start genomic position of the CNA region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNA region, 1-based and
          inclusive.
        - "clone" (str): clone ID.
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
    out_dir : str
        The output folder.
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
    cna_mode : {"hap-aware", "hap-unknown"}
        The mode of CNA profiles.
        - hap-aware: haplotype/allele aware.
        - hap-unknown: haplotype/allele unknown.
    barcode_whitelist_fn : str or None, default None
        File containing whitelist cell barcodes to be sampled for simulated
        data.
        If None, use randomly generated cell barcodes.
    ncores : int, default 1
        The number of cores/sub-processes.
    verbose : bool, default False
        Whether to show detailed logging information.
    kwargs_fit_sf : dict or None, default None
        The additional kwargs passed to function 
        :func:`~marginal.fit_libsize_wrapper` for fitting size factors.
        The available arguments are:
        - dist : {"lognormal", "swr", "normal", "t"}
            Type of distribution.
        If None, set as `{}`.
    kwargs_fit_rd : dcit or None, default None
        The additional kwargs passed to function 
        :func:`~marginal.fit_RD_wrapper` for fitting read depth.
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
        If None, set as `{}`.

    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    conf = Config()
    conf.count_fn = count_fn
    conf.clone_anno_fn = clone_anno_fn
    conf.cna_profile_fn = cna_profile_fn
    conf.out_dir = out_dir

    conf.size_factor = size_factor
    conf.marginal = marginal
    conf.libsize_ratio = libsize_ratio
    conf.loss_allele_freq = loss_allele_freq
    conf.cna_mode = cna_mode
    conf.barcode_whitelist_fn = barcode_whitelist_fn
    conf.ncores = ncores
    conf.verbose = verbose

    conf.kwargs_fit_sf = {} if kwargs_fit_sf is None else kwargs_fit_sf
    conf.kwargs_fit_rd = {} if kwargs_fit_rd is None else kwargs_fit_rd

    ret, res = cs_run(conf)
    return((ret, res))



def sample_barcodes(fn, n, suffix = "-1", sort = True):
    """Sample barcodes from file.
    
    Parameters
    ----------
    fn : str
        File containing cell barcodes to be sampled.
    n : int
        Number of barcodes to generate.
    suffix : str, default "-1"
        Suffix appended to the barcodes.
    sort : bool, default True
        Whether to sort the generated barcodes.
    
    Returns
    -------
    numpy.ndarray
        Sampled barcodes.
    """
    seed = np.loadtxt(fn, dtype = str)
    assert seed.shape[0] >= n
        
    # Note, standard random.sample() is much more efficient than
    # numpy.random.choice() when k << sample-space-size.
    x = np.array(random.sample(seed.tolist(), k = n))
    if sort:
        x = np.sort(x)
    if suffix:
        x = np.array([i + suffix for i in x])
    return(x)
