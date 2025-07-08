# copy.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
import shutil
import sys

from logging import info
from .io import load_params, save_params
from .marginal import fit_RD, simu_RD
from .pp import calc_size_factors, clone_calc_n_cell_each, \
    cna_get_overlap_features, qc_libsize, subset_adata_by_cell_types
from ..utils.cdata import sum_layers
from ..utils.clone import load_clones
from ..utils.gcna import load_cnas
from ..xlib.xdata import array_to_sparse, load_h5ad



def cs_init(conf):
    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")


    assert os.path.exists(conf.count_fn)
    adata = load_h5ad(conf.count_fn)
    for obs_key in ("cell", "cell_type"):
        assert obs_key in adata.obs.columns
    for var_key in ("feature", "chrom", "start", "end"):
        assert var_key in adata.var.columns
    info("load count data, shape = %s." % str(adata.shape))
    
    
    # here use "csr" to make row (cell) slicing efficient.
    adata = array_to_sparse(adata, which = "csr", layers = conf.alleles)


    assert conf.cna_mode in ("hap-aware", "hap-unknown")
    if conf.cna_mode == "hap-aware":
        for ale in conf.alleles:
            assert ale in adata.layers
        adata.X = sum_layers(adata, layers = conf.alleles)
        info("set adata.X as sum of the allele layers.")
    else:
        assert adata.X is not None


    assert os.path.exists(conf.cna_profile_fn)
    cna_profile = load_cnas(
        conf.cna_profile_fn, sep = "\t", cna_mode = conf.cna_mode)
    info("load CNA profiles, shape = %s." % str(cna_profile.shape))


    assert os.path.exists(conf.clone_anno_fn)
    clone_anno = load_clones(conf.clone_anno_fn, sep = "\t")
    info("load clone annotations, shape = %s." % str(clone_anno.shape))


    os.makedirs(conf.out_dir, exist_ok = True)

    if conf.barcode_whitelist_fn is not None:
        if conf.barcode_whitelist_fn.lower() == "none":
            conf.barcode_whitelist_fn = None
        else:
            assert os.path.exists(conf.barcode_whitelist_fn)


    if conf.size_factor is not None:
        assert conf.size_factor in ("libsize", )
    
    assert conf.marginal in ("auto", "poi", "nb", "zinb")


    kwargs_fit_sf = conf.def_kwargs_fit_sf.copy()
    for k, v in kwargs_fit_sf.items():
        if k in conf.kwargs_fit_sf:
            kwargs_fit_sf[k] = conf.kwargs_fit_sf[k]
    conf.kwargs_fit_sf = kwargs_fit_sf
    #info("complete kwargs_fit_sf is: %s." % str(conf.kwargs_fit_sf))

    
    kwargs_fit_rd = conf.def_kwargs_fit_rd.copy()
    for k, v in kwargs_fit_rd.items():
        if k in conf.kwargs_fit_rd:
            kwargs_fit_rd[k] = conf.kwargs_fit_rd[k]
    conf.kwargs_fit_rd = kwargs_fit_rd
    if "min_nonzero_num" in conf.kwargs_fit_rd:
        assert len(conf.kwargs_fit_rd["min_nonzero_num"]) == len(conf.alleles)
    #info("complete kwargs_fit_rd is: %s." % str(conf.kwargs_fit_rd))
    
    
    info("updated configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")

    
    data = dict(
        # adata : anndata.AnnData
        #   The count matrices.
        adata = adata,

        # clone_anno : pandas.DataFrame
        #   The clone annotations.
        clone_anno = clone_anno,

        # cna_profile : pandas.DataFrame
        #   The clonal CNA profile.
        cna_profile = cna_profile
    )
    return(data)



def cs_pp(data, conf, out_dir):
    adata = data["adata"]
    cna_profile = data["cna_profile"]
    clone_anno = data["clone_anno"]
    
    # check args.
    os.makedirs(out_dir, exist_ok = True)

    
    # calc number of cells in each clone before filtering seed cells.
    n_cell_each = clone_calc_n_cell_each(    # list of int
        clone_anno = clone_anno,
        adata = adata
    )
    info("#cells in each clone: %s." % str(n_cell_each))

    
    # subset adata (count matrices) by cell types.
    # only keep cell types listed in clone annotations.
    adata = subset_adata_by_cell_types(adata, clone_anno)
    info("subset adata by cell types. current shape = %s." % str(adata.shape))


    # filter low-quality cells, e.g, with very small library size or small
    # number of expressed features.
    adata, n_cells_filtered = qc_libsize(
        adata, 
        conf, 
        out_dir = out_dir,
        out_prefix = "pp"
    )
    info("QC: %d cells filtered. current shape = %s." %  \
         (n_cells_filtered, str(adata.shape)))
    

    cna_clones = np.unique(cna_profile["clone"])
    all_clones = np.unique(clone_anno["clone"])
    assert np.all(np.isin(cna_clones, all_clones))
    info("there are %d CNA clones in all %d clones." % (
        len(cna_clones), len(all_clones)))


    cna_fet = cna_get_overlap_features(    # dict of {reg_id (str) : feature indexes (list of int)}
        cna_profile = cna_profile,
        adata = adata
    )
    info("overlapping features extracted for %d CNA records." % \
        cna_profile.shape[0])


    # get cell-wise size factors.
    size_factors_train, size_factors_simu = calc_size_factors(
        adata = adata,
        size_factor = conf.size_factor,
        clone_cell_types = clone_anno["cell_type"],
        n_cell_each = n_cell_each,
        kwargs_fit_sf = conf.kwargs_fit_sf,
        verbose = conf.verbose
    )
    info("size factors calculated.")


    data["adata"] = adata.copy()
    res = dict(
        # n_cell_each : list of int
        #   Number of cells in each clone.
        n_cell_each = n_cell_each,
        
        # cna_fet : dict of {str : numpy.ndarray}
        #   The indices of overlapping features of each CNA region.
        cna_fet = cna_fet,
        
        # size_factors_train : numpy.ndarray
        #   The size factors for cells in `adata`.
        size_factors_train = size_factors_train,
        
        # size_factors_test : list of numpy.ndarray
        #   The clone-specific size factors.
        #   Its length and order match `n_cell_each`.
        #   Its elements are size factors of cells in corresponding clone.
        size_factors_simu = size_factors_simu
    )
    return((data, res))



def cs_cna_core(
    count_fn,
    out_dir,
    out_prefix,
    clones,
    cell_types,
    n_cell_each,
    cn_ratio,
    size_factor,
    size_factors_train,
    size_factors_simu,
    marginal,
    kwargs_fit_rd,
    libsize_ratio,
    ncores,
    verbose
):
    """Core function of simulating CNA counts trained on *cell x feature*
    count matrix.
    
    Parameters
    ----------
    count_fn : str
        The :class:`~anndata.AnnData` file storing *cell x feature* matrix in
        its ".X".
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    clones : pandas.Series
        The ID of CNA clones.
    cell_types : pandas.Series
        The source cell types used by `clones`.
    n_cell_each : list of int
        Number of cells in each of `clones`.
    cn_ratio : dict of {str : numpy.ndarray}
        Keys are clones, values are feature-specific copy ratios of 
        corresponding clone.
        The copy ratio, e.g., 1.0 for copy neutral; >1.0 for copy gain;
        and <1.0 for copy loss.
        Note that you can specify clones with copy number alterations only,
        since all features are assumed have ratio 1.0 unless specified.
    size_factor : str or None, default "libsize"
        The type of size factor.
        Currently, only support "libsize" (library size).
        Set to `None` if do not use size factors for model fitting.
    size_factors_train : numpy.ndarray of float or None
        The cell-wise size factors from trainning data.
        Its order should match '.obs["cell"]' in `count_fn`.
        None means do not use it.
    size_factors_simu : list of numpy.ndarray of float or None
        The cell-wise simulated size factors in each of `clones`.
        Its length and order should match `clones`.
        None means do not use it.
    marginal : {"auto", "poi", "nb", "zinb"}
        Type of marginal distribution.
    kwargs_fit_rd : dict
        The additional kwargs passed to function 
        :func:`~.marginal.fit_RD` for fitting read depth.
    libsize_ratio : float, default 1.0
        Ratio of library size of simulated cells compared to seed cells.
    ncores : int
        Number of cores.
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    anndata.AnnData
        Simulated *cell x feature* RD values.
        It has one column "cell_type" in `.obs` and one column "feature" in
        `.var`.
    """
    # check args.
    os.makedirs(out_dir, exist_ok = True)

    clones = clones.values
    assert len(cell_types) == len(clones)
    assert len(n_cell_each) == len(clones)
    for clone in cn_ratio.keys():
        assert clone in clones

    if size_factors_simu is not None:
        assert len(size_factors_simu) == len(clones)
        

    info("start fit RD ...")
    
    fit_dir = os.path.join(out_dir, "tmp_fit_RD")
    os.makedirs(fit_dir, exist_ok = True)
    
    params, features = fit_RD(
        count_fn = count_fn,
        tmp_dir = fit_dir,
        cell_type_fit = np.unique(cell_types),
        s = size_factors_train,
        s_type = size_factor,
        marginal = marginal,
        ncores = ncores,
        verbose = verbose,
        **kwargs_fit_rd
    )
    
    params_fn = os.path.join(out_dir, "%s.fit.output.params.pickle" % out_prefix)
    save_params(params, params_fn)
    del params
    gc.collect()
    
    
    info("start simulate RD ...")
    
    simu_dir = os.path.join(out_dir, "tmp_simu_RD")
    os.makedirs(simu_dir, exist_ok = True)
    
    adata, params = simu_RD(
        params_fn = params_fn,
        features = features,
        tmp_dir = simu_dir,
        cell_type_new = clones,
        cell_type_old = cell_types,
        n_cell_each = n_cell_each,
        s = size_factors_simu,
        cn_ratio = cn_ratio,
        total_count_new = None,
        libsize_ratio = libsize_ratio,
        dtype = np.int32,
        ncores = ncores, 
        verbose = verbose
    )
    
    params_fn = os.path.join(out_dir, "%s.simu.output.params.pickle" % out_prefix)
    save_params(params, params_fn)
    del params
    gc.collect()    
    
    
    info("clean tmp files ...")
    
    shutil.rmtree(fit_dir)
    shutil.rmtree(simu_dir)

    return(adata)
