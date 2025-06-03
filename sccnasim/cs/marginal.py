# marginal.py - fit gene-wise marginal distributions and then simulate.


import anndata as ad
import copy
import gc
import multiprocessing
import numpy as np
import os
import pandas as pd
import scipy as sp

from collections import OrderedDict
from logging import info, error, debug
from .io import load_params, save_params
from ..xlib.xbase import is_scalar_numeric
from ..xlib.xdata import load_h5ad, save_h5ad
from ..xlib.xio import load_pickle, save_pickle
from ..xlib.xmath import   \
    estimate_dist_nb, estimate_dist_poi,  \
    fit_dist_nb, fit_dist_poi, fit_dist_zinb, fit_dist_zip,  \
    rand_zinb
from ..xlib.xmatrix import sparse2array, array2sparse
from ..xlib.xthread import mp_error_handler, split_n2batch



ALL_DIST = ("auto", "zinb", "nb", "poi")



### Marginals of features

def __fit_dist_wrapper(mar, *args, **kwargs):
    func = None
    if mar == "nb":     func = fit_dist_nb
    elif mar == "poi":  func = fit_dist_poi
    elif mar == "zinb": func = fit_dist_zinb
    elif mar == "zip":  func = fit_dist_zip
    else:
        error("invalid marginal '%s'." % mar)
        raise ValueError
    ret, par, stat = func(*args, **kwargs)
    mres = {
        "loglik": stat["loglik"] if stat else None,
        "converged": stat["converged"] if stat else None
    }
    return((ret, par, mres))


def __fit_dist_nb(*args, **kwargs):
    return(__fit_dist_wrapper("nb", *args, **kwargs))

def __fit_dist_poi(*args, **kwargs):
    return(__fit_dist_wrapper("poi", *args, **kwargs))

def __fit_dist_zinb(*args, **kwargs):
    return(__fit_dist_wrapper("zinb", *args, **kwargs))

def __fit_dist_zip(*args, **kwargs):
    return(__fit_dist_wrapper("zip", *args, **kwargs))
    


# TODO: check the significance level of over-dispersion (parameter alpha)?
def fit_RD_feature(
    x,
    s = None,
    marginal = "auto",
    max_iter = 1000,
    pval_cutoff = 0.05
):
    """Fit one feature.

    Parameters
    ----------
    x : numpy.ndarray
        The vector of counts.
    s : float or None, default None
        The size factor, typically library size.
        Set to `None` if do not use it.
    marginal
    max_iter
    pval_cutoff
        See :func:`fit_RD()` for details.

    Returns
    -------
    flag : int
        Bit-wise return code:
        0  - any model is selected?
        1  - selected model is from fitted model (i.e., has loglik)?
             otherwise, is from heuristic model?
        2  - Poisson was fitted?
        3  - NB was fitted?
        4  - ZIP was fitted?
        5  - ZINB was fitted?
        6  - Poisson fit failed?
        7  - NB fit failed?
        8  - ZIP fit failed?
        9  - ZINB fit failed?
        10 - Poisson not converged?
        11 - NB not converged?
        12 - ZIP not converged?
        13 - ZINB not converged?
    model : str
        Model name, should be one of {"poi", "nb", "zip", "zinb"}.
    par : dict
        Fitted parameters.
    mres : dict
        Model result. It includes several keys, including "loglik",
        "converged" etc.
    stat : dict
        Statistics of the feature. It includes several keys, including
        "mean" (mean), "var" (variance), "min" (min), "max" (max), 
        and "median" (median).

    Examples
    --------
    # test poisson fitting
    >>> x = np.random.poisson(3, size = 100)
    >>> res = rdr.mar.fit_RD_feature(x)
    >>> print(res)

    # test NB fitting
    >>> mu, alpha = 5, 2
    >>> var = mu + alpha * mu ** 2
    >>> n = mu ** 2 / (var - mu)
    >>> p = mu / var
    >>> x = np.random.negative_binomial(n, p, size = 1000)
    >>> res = rdr.mar.fit_RD_feature(x, max_iter = 1000, verbose = False)
    >>> print(res)
    """
    if marginal not in ALL_DIST:
        error("invalid marginal '%s'." % marginal)
        raise ValueError

    flag = 0
    model, par = None, None
    mres = {"loglik": None, "converged": None}

    m, v = np.mean(x), np.var(x)
    stat = {"mean": m, "var": v, "min": np.min(x), "max": np.max(x),
            "median": np.median(x)}

    verbose = False
    while True:
        if m <= 0.0:
            model, par = "poi", {"infl": 0, "disp": 0, "mu": 0}
            break

        if marginal == "auto" or marginal == "zinb":
            if m >= v:                     # no over-dispersion.
                ret_poi, par_poi, mres_poi = __fit_dist_poi(
                    x, s, max_iter, verbose)
                flag |= (1 << 2)
                if ret_poi != 0:
                    flag |= (1 << 6)
                    model, par = "poi", estimate_dist_poi(x, s)
                    break
                elif mres_poi["converged"] is False:
                    flag |= (1 << 10)
                    model, par = "poi", estimate_dist_poi(x, s)
                    break                    
                if np.min(x) > 0:          # no zero-inflation.
                    model, par, mres = "poi", par_poi, mres_poi
                    break
                
                ret_zip, par_zip, mres_zip = __fit_dist_zip(
                    x, s, max_iter, verbose)
                flag |= (1 << 4)
                if ret_zip != 0:
                    flag |= (1 << 8)
                    model, par, mres = "poi", par_poi, mres_poi
                    break
                elif mres_zip["converged"] is False:
                    flag |= (1 << 12)
                    model, par = "poi", estimate_dist_poi(x, s)
                    break

                # check significance level of zero-inflation via GLR test.
                chisq_val = 2 * (mres_zip["loglik"] - mres_poi["loglik"])
                p_val = 1 - sp.stats.chi2.cdf(chisq_val, df = 1)
                if p_val < pval_cutoff:
                    model, par, mres = "zip", par_zip, mres_zip
                else:
                    model, par, mres = "poi", par_poi, mres_poi
                break
            
            else:                          # potential over-dispension.
                mres_nb = None
                if marginal == "auto" or np.min(x) > 0:
                    ret_nb, par_nb, mres_nb = __fit_dist_nb(
                        x, s, max_iter, verbose)
                    flag |= (1 << 3)
                    if ret_nb != 0:
                        flag |= (1 << 7)
                        model, par = "poi", estimate_dist_poi(x, s)
                        break
                    elif mres_nb["converged"] is False:
                        flag |= (1 << 11)
                        model, par = "poi", estimate_dist_poi(x, s)
                        break
                    if np.min(x) > 0:      # no zero-inflation.
                        model, par, mres = "nb", par_nb, mres_nb
                        break
                    
                # assert (marginal == "auto" and np.min(x) <= 0) or \
                #       (marginal == "zinb" and np.min(x) <= 0)

                ret_zinb, par_zinb, mres_zinb = __fit_dist_zinb(
                    x, s, max_iter, verbose)
                flag |= (1 << 5)
                if ret_zinb != 0:
                    flag |= (1 << 9)
                    if mres_nb is None:
                        ret_nb, par_nb, mres_nb = __fit_dist_nb(
                            x, s, max_iter, verbose)
                        flag |= (1 << 3)
                        if ret_nb != 0:
                            flag |= (1 << 7)
                            model, par = "poi", estimate_dist_poi(x, s)
                            break
                        elif mres_nb["converged"] is False:
                            flag |= (1 << 11)
                            model, par = "poi", estimate_dist_poi(x, s)
                            break
                    model, par, mres = "nb", par_nb, mres_nb
                    break
                elif mres_zinb["converged"] is False:
                    flag |= (1 << 13)
                    model, par = "poi", estimate_dist_poi(x, s)
                    break
                
                if marginal == "auto":      # assert mres_nb is not None
                    chisq_val = 2 * (mres_zinb["loglik"] - mres_nb["loglik"])
                    p_val = 1 - sp.stats.chi2.cdf(chisq_val, df = 1)
                    if p_val < pval_cutoff:
                        model, par, mres = "zinb", par_zinb, mres_zinb
                    else:
                        model, par, mres = "nb", par_nb, mres_nb
                    break
                else:
                    model, par, mres = "zinb", par_zinb, mres_zinb
                    break
                
        elif marginal == "nb":
            if (m >= v):
                model, par = "poi", estimate_dist_poi(x, s)
                break
            else:
                ret_nb, par_nb, mres_nb = __fit_dist_nb(
                    x, s, max_iter, verbose)
                flag |= (1 << 3)
                if ret_nb != 0:
                    flag |= (1 << 7)
                    model, par = "poi", estimate_dist_poi(x, s)
                    break
                elif mres_nb["converged"] is False:
                    flag |= (1 << 11)
                    model, par = "poi", estimate_dist_poi(x, s)
                    break
                model, par, mres = "nb", par_nb, mres_nb
                break
        
        elif marginal == "poi":
            model, par = "poi", estimate_dist_poi(x, s)
            break
        
        else:
            error("invalid marginal '%s'." % marginal)
            raise ValueError

    if model is not None:
        flag |= (1 << 0)
    if mres and mres["loglik"] is not None:
        flag |= (1 << 1)

    return((flag, model, par, mres, stat))


def __fit_RD_feature(
    x,
    index,
    s = None,
    marginal = "auto",
    max_iter = 100,
    pval_cutoff = 0.05
):
    flag, model, par, mres, stat = fit_RD_feature(
        x = x,
        s = s,
        marginal = marginal,
        max_iter = max_iter,
        pval_cutoff = pval_cutoff
    )
    return((flag, model, par, mres, stat, index))



def __fit_RD_cell_type_batch(
    count_fn,
    idx_b = None,      # 0-based transcriptomics-scale index of the first feature in this batch.
    s = None,
    marginal = "auto",
    min_nonzero_num = 3,
    max_iter = 100,
    pval_cutoff = 0.05    
):
    adata = load_h5ad(count_fn)
    X = adata.X         # should be csc_array or csc_matrix.
    n, p = X.shape
    
    result = []
    fet_idx = {
        "nz": set(),               # Indexes of non-zero features.
        "oth": set()               # Indexes of other features.
    }
    for i in range(p):
        x = X[:, i].toarray()
        idx = idx_b + i
        if (x > 0).sum() >= min_nonzero_num:
            fet_idx["nz"].add(idx)
            res = __fit_RD_feature(
                x = x,
                index = idx,
                s = s,
                marginal = marginal,
                max_iter = max_iter,
                pval_cutoff = pval_cutoff
            )
            result.append(res)
        else:
            fet_idx["oth"].add(idx)
            
    del adata
    gc.collect()
    
    return((result, fet_idx))
    


def fit_RD_cell_type(
    count_fn,
    tmp_dir,
    s,
    s_type,
    marginal = "auto",
    min_nonzero_num = 3,
    ncores = 1,
    max_iter = 100,
    pval_cutoff = 0.05,
    verbose = True
):
    """Fit all features in one cell type.

    Parameters
    ----------
    count_fn : str
        AnnData file storing *cell x feature* count matrix.
    tmp_dir : str
        Folder to store temporary data.
    s
    s_type
    marginal
    min_nonzero_num
    ncores
    max_iter
    pval_cutoff
    verbose
        See :func:`fit_RD()` for details.

    Returns
    -------
    dict
        The fitted parameters, will be used by downstream simulation.
    """
    if verbose:
        info("start ...")
    
    # check args.
    adata = load_h5ad(count_fn)
    os.makedirs(tmp_dir, exist_ok = True)

    if marginal not in ALL_DIST:
        error("invalid marginal '%s'." % marginal)
        raise ValueError
    
    n, p = adata.shape
    if s is None:
        assert s_type is None
    if s_type is None:
        assert s is None
    if s is not None:
        assert len(s) == n
        
    n_read = adata.X.sum()

    
    # feature-specific counts fitting.
    if verbose:
        info("processing %d features in %d cells (ncores = %d) ..." % \
            (p, n, ncores))
        
    # here, use "csc" to make column (feature) slicing efficient.
    adata.X = array2sparse(adata.X, which = "csc")
    
    # split features into batches.
    # here, max_n_batch to account for max number of open files.
    bd_m, bd_n, bd_indices = split_n2batch(p, ncores, max_n_batch = 1000)

    if verbose:
        info("split features into %d batches for multiprocessing." % bd_m)

    count_fn_list = []
    for i, (b, e) in enumerate(bd_indices):
        fn = os.path.join(tmp_dir, "fet.b%d.counts.h5ad" % i)
        save_h5ad(adata[:, b:e], fn)
        count_fn_list.append(fn)
    del adata
    gc.collect()
    adata = None


    # multi-processing features.
    if verbose:
        info("begin multiprocessing with %d cores." % min(ncores, bd_m))
        
    mp_res = []
    pool = multiprocessing.Pool(processes = min(ncores, bd_m))
    for i, (b, e) in enumerate(bd_indices):
        fn = count_fn_list[i]
        mp_res.append(pool.apply_async(
            __fit_RD_cell_type_batch,
            kwds = dict(
                count_fn = fn,
                idx_b = b,
                s = s,
                marginal = marginal,
                min_nonzero_num = min_nonzero_num,
                max_iter = max_iter,
                pval_cutoff = pval_cutoff
            ),
            callback = None,
            error_callback = mp_error_handler
        ))
    pool.close()
    pool.join()
    mp_res = [res.get() for res in mp_res]

    if verbose:
        info("multi-processing finished.")
        info("merge results ...")

    result = []
    fet_idx = {}
    for res, idx in mp_res:
        result.extend(res)
        for k in idx.keys():
            if k not in fet_idx:
                fet_idx[k] = set()
            fet_idx[k].update(idx[k])
    

    # TODO: implement a class for cell type specific params.
    params_nz = []
    params = {
        "params_nz": None,                    # params for non-zero features (pd.DataFrame).
        "fet_idx_nz": fet_idx["nz"],          # index (0-based) of non-zero features (set).
        "fet_idx_oth": fet_idx["oth"],        # index (0-based) of other features (set).
        "size_factor_type": s_type,           # size factor type (str).
        "size_factor_value": s,               # size factor values (np.array or None).
        "min_nonzero_num": min_nonzero_num,   # min number of non-zero entries (int).
        "n_cell": n,                          # number of cells (int).
        "n_read": n_read                      # number of reads (int).
    }
    for flag, model, par, mres, stat, index in result:
        if flag & 0x1:      # any model is selected.
            params_nz.append({
                "flag": flag,
                "model": model,
                "par": par,
                "mres": mres,
                "stat": stat,
                "index": index
            })
        else:
            params["fet_idx_nz"].remove(index)
            params["fet_idx_oth"].add(index)

    assert len(params_nz) == len(params["fet_idx_nz"])
    params_nz = {
        "index": [res["index"] for res in params_nz],
        "mu": [res["par"]["mu"] for res in params_nz],
        "dispersion": [res["par"]["disp"] for res in params_nz],
        "inflation": [res["par"]["infl"] for res in params_nz],
        "mean": [res["stat"]["mean"] for res in params_nz],
        "var": [res["stat"]["var"] for res in params_nz],
        "min": [res["stat"]["min"] for res in params_nz],
        "max": [res["stat"]["max"] for res in params_nz],
        "median": [res["stat"]["median"] for res in params_nz],
        "flag": [res["flag"] for res in params_nz],
        "model": [res["model"] for res in params_nz],
        "loglik": [res["mres"]["loglik"] for res in params_nz],
        "converged": [res["mres"]["converged"] for res in params_nz]
    }
    params_nz = pd.DataFrame(data = params_nz)
    params_nz = params_nz.sort_values(by = ["index"])
    params["params_nz"] = params_nz
    
    return(params)



# TODO: consider Cox–Reid bias adjustment (e.g., in the DESeq2 paper).
def fit_RD(
    count_fn,
    tmp_dir,
    cell_type_fit = None,
    s = None,
    s_type = None,
    marginal = "auto",
    min_nonzero_num = 3,
    ncores = 1,
    max_iter = 100,
    pval_cutoff = 0.05,
    verbose = True
):
    """Fit all features in all cell types.

    Parameters
    ----------
    count_fn : str
        AnnData file containing the *cell x feature* count matrix.
        It should have a column `cell_type` in `.obs`, and a column `feature`
        in `.var`.
    tmp_dir : str
        Folder to store temporary data.
    cell_type_fit : list of str or None, default None
        A list of cell types (str) whose features will be fitted.
        If `None`, use all unique cell types in `count_fn`.
    s : numpy.ndarray of float or None, default None
        The size factors of cells in `count_fn`.
        None means do not use it.
    s_type : str or None, default None
        The type of size factors.
        None means do not use it.
    marginal : {"auto", "poi", "nb", "zinb"}
        Type of marginal distribution.
        One of "auto" (auto select), "poi" (Poisson), 
        "nb" (Negative Binomial),
        and "zinb" (Zero-Inflated Negative Binomial).
    min_nonzero_num : int, default 3
        The minimum number of cells that have non-zeros for one feature.
        If smaller than the cutoff, then the feature will not be fitted
        (i.e., its mean will be directly treated as 0).
    ncores : int, default 1
        The number of cores/sub-processes.
    max_iter : int, default 100
        Number of maximum iterations in model fitting.
    pval_cutoff : float, default 0.05
        The p-value cutoff for model selection with GLR test.
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    OrderedDict
        The fitted parameters, will be used by downstream simulation.
        In each item (pair), the key is the cell type (str) and the value
        is the cell-type-specific parameters returned by 
        :func:`~cs.marginal.fit_RD_cell_type`.
    numpy.ndarray of str
        The feature names from `count_fn`.
    """
    if verbose:
        info("start ...")

    # check args
    adata = load_h5ad(count_fn)
    os.makedirs(tmp_dir, exist_ok = True)

    #X = sparse2array(adata.X)
    X = adata.X            # should be "csr_array" or "csr_matrix".
    n, p = X.shape

    assert "cell_type" in adata.obs.columns
    assert "feature" in adata.var.columns

    cell_types = adata.obs["cell_type"].values
    all_cell_types = list(set(cell_types))
    if cell_type_fit is None:
        cell_type_fit = sorted(all_cell_types)
    else:
        assert len(cell_type_fit) == len(set(cell_type_fit))
        assert np.all(np.isin(cell_type_fit, all_cell_types))

        
    if s is None:
        assert s_type is None
    if s_type is None:
        assert s is None
    if s is not None:
        assert len(s) == n

        
    if marginal not in ALL_DIST:
        error("invalid marginal '%s'." % marginal)
        raise ValueError
        
    features = np.array(adata.var["feature"])

    
    # process by cell type.
    info("fit %d features in %d cell types (ncores = %d) ..." %  \
            (p, len(cell_type_fit), ncores))

    # split cells into batches by cell types.
    count_fn_list = []
    for i, c in enumerate(cell_type_fit):
        fn = os.path.join(tmp_dir, "celltype.b%d.counts.h5ad" % i)
        adata_s = adata[cell_types == c, :]
        save_h5ad(adata_s, fn)
        count_fn_list.append(fn)
    del adata
    gc.collect()
    adata = None


    # model fitting
    params = OrderedDict()
    for c, fn in zip(cell_type_fit, count_fn_list):
        info("fitting RD for cell type '%s'." % c)

        c_dir = os.path.join(tmp_dir, c)
        os.makedirs(c_dir, exist_ok = True)
        
        idx = cell_types == c
        par = fit_RD_cell_type(
            count_fn = fn,
            tmp_dir = c_dir,
            s = s[idx] if s is not None else None,
            s_type = s_type,
            marginal = marginal,
            min_nonzero_num = min_nonzero_num,
            ncores = ncores,
            max_iter = max_iter,
            pval_cutoff = pval_cutoff,
            verbose = verbose
        )
        assert len(par["fet_idx_nz"]) + len(par["fet_idx_oth"]) == p
        params[c] = par


    info("fitting statistics:")
    df = pd.DataFrame(data = {
        "cell_type": cell_type_fit,
        "fet_idx_nz": [len(r["fet_idx_nz"]) for r in params.values()],
        "fet_idx_oth": [len(r["fet_idx_oth"]) for r in params.values()]
    })
    info(str(df))

    return((params, features))



def simu_RD_feature(params, n, s = None, s_type = None):
    """Simulate RD values for one feature.
    
    Parameters
    ----------
    params : dict
        The distribution parameters.
    n : int
        Number of cells.
    s : numpy.ndarray of float or None, default None
        The size factor.
        Its length should be `n`.
        Set to `None` if do not use it.
    s_type : str or None, default None
        The type of size factors.
        Currently only "libsize" is supported.
        Set to `None` if do not use it.
    
    Returns
    -------
    numpy.ndarray of int
        Simulated RD values of length `n`.
    """
    mu, disp, infl = [params[k] for k in ("mu", "dispersion", "inflation")]
    if s is not None:
        if s_type == "libsize":
            mu = s * mu
        else:
            error("invalid size factor '%s'." % s_type)
            raise ValueError
    dat = rand_zinb(
        mu = mu + 0.0,
        alpha = disp + 0.0,
        infl = infl + 0.0,
        size = n
    )
    return(dat)



def __simu_RD_feature(params, n, s, s_type, index):
    dat = simu_RD_feature(params, n, s, s_type)
    return((dat, index))



def __simu_RD_cell_type_batch(
    params_fn,
    n,
    s = None,
    s_type = None
):  
    params = load_params(params_fn)
    df = params
    p = df.shape[0]
    
    result = []
    for i in range(p):
        par = df.iloc[i, ]           # feature-specific params.
        res = __simu_RD_feature(
            params = par,
            n = n,
            s = s,
            s_type = s_type,
            index = par["index"] + 0
        )
        result.append(res)
        
    del params
    del df
    gc.collect()

    return(result)



def simu_RD_cell_type(
    params_fn, 
    tmp_dir,
    n,
    s = None, dtype = np.int32, ncores = 1, verbose = False
):
    """Simulate RD values for all features in one cell type.
    
    Parameters
    ----------
    params_fn : str
        File storing the cell-type-specific parameters fitted in
        :func:`fit_RD_cell_type`.
    tmp_dir : str
        Folder to store temporary data.
    n : int
        Number of cells to be simulated in this cell type.
    s : numpy.ndarray of float or None, default None
        The size factors for simulated cells of this cell type.
        Its length should be `n`.
        Set to `None` if do not use it.
    dtype
        The dtype of the simulated matrix.
    ncores : int, default 1
        Number of cores.
    verbose : bool, default False
        Whether to show detailed logging information.
    
    Returns
    -------
    numpy.ndarray of int
        Simulated *cell x feature* RD values.
    """
    if verbose:
        info("start ...")

    # check args
    params = load_params(params_fn)
    os.makedirs(tmp_dir, exist_ok = True)

    assert len(params["params_nz"]) == len(params["fet_idx_nz"])
    p_nz = len(params["fet_idx_nz"])
    p_oth = len(params["fet_idx_oth"])
    p = p_nz + p_oth

    s_type = params["size_factor_type"]
    if s_type is None and s is not None:
        error("size factors unused.")
        raise ValueError
    elif s_type is not None and s is None:
        error("size factors missing.")
        raise ValueError

    if s is not None:
        assert len(s) == n
    
    if p_nz <= 0:
        return(np.zeros((n, p)))


    # feature-specific count simulation.
    if verbose:
        info("simulate %d features in %d cells (ncores = %d) ..." %  \
            (p, n, ncores))
        
    # split features into batches.
    # here, max_n_batch to account for max number of open files.
    bd_m, bd_n, bd_indices = split_n2batch(p_nz, ncores, max_n_batch = 1000)

    if verbose:
        info("split features into %d batches for multiprocessing." % bd_m)
               
    params_fn_list = []
    for i, (b, e) in enumerate(bd_indices):
        fn = os.path.join(tmp_dir, "fet.b%d.pickle" % i)
        save_params(params["params_nz"].iloc[b:e, ], fn)
        params_fn_list.append(fn)
    del params
    gc.collect()
    params = None


    # multiprocessing.
    if verbose:
        info("begin multiprocessing with %d cores." % min(ncores, bd_m))

    pool = multiprocessing.Pool(processes = min(ncores, bd_m))
    mp_res = []
    for fn in params_fn_list:
        mp_res.append(pool.apply_async(
            __simu_RD_cell_type_batch,
            kwds = dict(
                params_fn = fn,
                n = n,
                s = s,
                s_type = s_type
            ),
            callback = None,
            error_callback = mp_error_handler
        ))
    pool.close()
    pool.join()
    mp_res = [res.get() for res in mp_res]

    if verbose:
        info("multi-processing finished.")
        info("merge results ...")
               
    result = []
    for res in mp_res:
        result.extend(res)

    # TODO:
    # - construct matrix in each batch and then merge here.
    # - use scipy.sparse.lil_array/matrix or dok_array/matrix instead of
    #   ndarray for efficiently constructing and modifying the sparse 
    #   structure.
    mtx = np.zeros((n, p))
    for dat, index in result:
        mtx[:, index] = dat
    mtx = mtx.astype(dtype)
    
    mtx = array2sparse(mtx, which = "csr")

    return(mtx)



def simu_RD(
    params_fn,
    features,
    tmp_dir,
    cell_type_new = None,
    cell_type_old = None,
    n_cell_each = None,
    s = None,
    cn_fold = None,
    total_count_new = None,
    libsize_ratio = 1.0,
    dtype = np.int32,
    ncores = 1, 
    verbose = False
):
    """Simulate RD values for all features in all cell types.
    
    Parameters
    ----------
    params_fn : str
        File storing the fitted parameters returned by :func:`fit_RD`.
    features : list of str
        A list of feature names. 
        Its order should match feature index in `params_fn`.
        Typically use the value returned by :func:`~cs.marginal.fit_RD`.
    tmp_dir : str
        Folder to store temporary data.
    cell_type_new : list of str or None, default None
        Cell type names for newly simulated cell clusters.
        Set to `None` to use all the old cell types (in training data).
    cell_type_old : list of str or None, default None
        The old cell types whose parameters (in `params`) will be used by 
        `cell_type_new`.
        Its length and order should match `cell_type_new`.
        Set to `None` to use all the old cell types (in training data).
        Note that when `cell_type_new` is not None, `cell_type_old` must be
        specified with valid values.
    n_cell_each : list of int or None, default None
        Number of cells in each new cell type (`cell_type_new`).
        Its length and order should match `cell_type_new`.
        Set to `None` to use #cells of old cell types (in training data).
    s : list of float or None, default None
        Cell-type-specific size factors.
        Its length and order should match `cell_type_new`.
        Its elements are vectors whose lengths matching elements of 
        `n_cell_each`.
        Set to `None` if do not use it.
    cn_fold : dict or None, default None
        The copy number (CN) fold, e.g., 1.0 for copy neutral; >1.0 for copy
        gain; and <1.0 for copy loss.
        Its keys are new cell types (str) and values are vectors of CN folds.
        For each such vector, length and order should be the same with
        `features`.
        Note that you can specify cell types with copy number alterations only,
        since all features are assumed have fold 1.0 unless specified.
        Set to `None` to use fold 1.0 on all features in all cell types.
    total_count_new : int, list of int or None, default None
        The total read counts to be simulated.
        If a int, it is the total libray size of all simulated cells; 
        If a list, it is a list of cell-type-specific total read counts whose 
        length and order should match `cell_type_new`.
        Set to `None` to set the scaling factor of total library size to 1.
    libsize_ratio : float, default 1.0
        Ratio of library size of simulated cells compared to seed cells.
        It will be used only when `total_count_new` is None.
    dtype
        The dtype of the simulated matrix.
    ncores : int, default 1
        Number of cores/sub-processes.
    verbose : bool, default False
        Whether to show detailed logging information.
    
    Returns
    -------
    anndata.AnnData
        Simulated RD values of *cell x feature*. 
        It has one column "cell_type" in `.obs` and one column "feature" 
        in `.var`.
    dict
        The updated `params` incorporating CN-folds, the same length as
        `cell_type_new`, while keeping the input `params` unchanged.
    """
    if verbose:
        info("start ...")

    # check args
    params = load_params(params_fn)
    params = copy.deepcopy(params)
    os.makedirs(tmp_dir, exist_ok = True)

    all_cell_types = list(params.keys())
    p = len(features)
    for c in all_cell_types:
        par = params[c]
        assert len(par["fet_idx_nz"]) + len(par["fet_idx_oth"]) == p

    if cell_type_new is None:
        cell_type_new = all_cell_types
        cell_type_old = all_cell_types
    else:
        assert cell_type_old is not None
        cell_type_new = list(cell_type_new)
    assert len(cell_type_new) == len(set(cell_type_new))
    assert len(cell_type_new) == len(cell_type_old)
    for c in cell_type_old:      # duplicates in `cell_type_old` are allowed.
        assert c in all_cell_types
    n_cell_types = len(cell_type_new)


    if n_cell_each is None:
        n_cell_each = [params[c]["n_cell"] for c in cell_type_old]
    else:
        assert len(n_cell_each) == len(cell_type_new)
    n_cell_new = np.sum(n_cell_each)
               

    info("simulate %d features in %d new cells from %d cell types (ncores = %d) ..." %  \
            (p, n_cell_new, len(cell_type_new), ncores))
    info("number of cells in each simulated cell type:\n\t%s." % \
             str(n_cell_each))


    if s is None:
        s = [None for _ in cell_type_new]
    else:
        assert len(s) == len(cell_type_new)
        for c_s, n_cell in zip(s, n_cell_each):
            if c_s is not None:
                assert len(c_s) == n_cell


    cn_fold_arr = None        # np.array (2d); cell_type x feature
    if cn_fold is None:
        cn_fold_arr = np.array([np.repeat(1.0, p) \
                            for _ in range(len(cell_type_new))])
    else:
        assert isinstance(cn_fold, dict)
        if verbose:
            info("CN folds are specified in %d cell types." % len(cn_fold))

        cn_fold_arr = np.array([np.repeat(1.0, p) \
                            for _ in range(n_cell_types)])
        for c, fold in cn_fold.items():
            if c not in cell_type_new:
                error("invalid cell type '%s' in cn_fold." % c)
                raise ValueError
            assert len(fold) == p
            idx = cell_type_new.index(c)
            cn_fold_arr[idx] = np.array(fold)           


    if total_count_new is None:
        pass
    else:
        assert is_scalar_numeric(total_count_new) or \
            len(total_count_new) == len(cell_type_new)

    total_count_old = np.array([params[c]["n_read"] for c in cell_type_old])
    n_cell_old = np.array([params[c]["n_cell"] for c in cell_type_old])
               

    # simulation
    # TODO: consider copy number fold when scaling to total_count_new.
    r = None                      # scaling factor
    if total_count_new is None:
        r = np.repeat(libsize_ratio, n_cell_types)
    elif is_scalar_numeric(total_count_new):
        r = np.repeat(
            total_count_new/np.sum(total_count_old / n_cell_old * n_cell_each),
            n_cell_types)
    else:
        # scDesign2: r = (total_count_new / n_cell_new) / \
        #                  (total_count_old / n_cell_old)
        r = (total_count_new / n_cell_each) / (total_count_old / n_cell_old)
               

    # split data into batches by cell types.
    params_fn_list = []
    bd_args_fn_list = []
    params_new = dict()
    for i, (c_new, c_old) in enumerate(zip(cell_type_new, cell_type_old)):
        par = copy.deepcopy(params[c_old])
        scaling = r[i] * cn_fold_arr[i][par["params_nz"]["index"]]
        par["params_nz"]["mu"] *= scaling

        par_fn = os.path.join(tmp_dir, "celltype.%s.input.params.pickle" % c_new)
        save_params(par, par_fn)
        params_fn_list.append(par_fn)
        
        c_dir = os.path.join(tmp_dir, c_new)
        os.makedirs(c_dir, exist_ok = True)

        args = dict(
            params_fn = par_fn,
            tmp_dir = c_dir,
            n = n_cell_each[i],
            s = s[i],
            dtype = dtype,
            ncores = ncores,
            verbose = verbose
        )
               
        fn = os.path.join(tmp_dir, "celltype.%s.input.args.pickle" % c_new)
        save_pickle(args, fn)
        bd_args_fn_list.append(fn)
               
        params_new[c_new] = copy.deepcopy(par)
               
    params_fn_new = os.path.join(tmp_dir, "updated.cn_fold.params.pickle")
    save_params(params_new, params_fn_new)
    
    del params_new
    del params
    del s
    del cn_fold
    del cn_fold_arr
    gc.collect()
    params_new = params = s = cn_old = cn_fold_arr = None
    

    # simulation in each cell type.
    mtx = None
    for i, (c_new, c_old, args_fn) in enumerate(zip(
        cell_type_new, cell_type_old, bd_args_fn_list)):
        info("simulate for new cell type '%s' based on '%s' ..." %  \
            (c_new, c_old))
        
        args = load_pickle(args_fn)
        c_mtx = simu_RD_cell_type(**args)
        if mtx is None:
            mtx = c_mtx
        else:
            mtx = sp.sparse.vstack([mtx, c_mtx])
    mtx = mtx.astype(dtype)

    adata = ad.AnnData(
        X = mtx,
        obs = pd.DataFrame(data = {
            "cell_type": np.repeat(cell_type_new, n_cell_each)}),
        var = pd.DataFrame(data = {"feature": features})
    )
    
    params_new = load_params(params_fn_new)

    return((adata, params_new))
