# hapidx.py - haplotype and its index.



def hap2idx(hap):
    """Get the haplotype index given the haplotype string.
    
    Parameters
    ----------
    hap : str
        The haplotype string.
    
    Returns
    -------
    int or None
        The haplotype index. None if the `hap` is invalid.
    """
    if hap in MAP_HAP_TO_IDX:
        return(MAP_HAP_TO_IDX[hap])
    else:
        return(None)
    
    

def idx2hap(idx):
    """Get the haplotype string given the haplotype index.
    
    Parameters
    ----------
    idx : int
        The haplotype index.
    
    Returns
    -------
    str or None
        The haplotype string. None if the `idx` is invalid.
    """
    if idx in MAP_IDX_TO_HAP:
        return(MAP_IDX_TO_HAP[idx])
    else:
        return(None)



# HAP : str
#   Haplotype.
# IDX : int
#   Haplotype index. 
#   Note that haplotype index is fine grained compared to coarse HAP, hence
#   HAP to IDX could be one-to-multi mapping, while IDX to HAP is one-to-one
#   mapping.

MAP_HAP_TO_IDX = {
    "A": [0],
    "B": [1],
    "D": [2],
    "O": [-1],
    "U": [-2]
}

MAP_IDX_TO_HAP = {
    0: "A",
    1: "B",
    2: "D",
    -1: "O",
    -2: "U"
}
