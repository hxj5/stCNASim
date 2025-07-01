# sam.py



import pysam


def detect_pe_mode(sam_fn):
    n = 100        # max number of test reads.
    f = 0.3        # min fraction of PE reads for 'PE' mode.
    i = 0          # number of iterated reads;
    j = 0          # number of test (valid) reads;
    k = 0          # number of test PE reads.
    
    sam = pysam.AlignmentFile(sam_fn, "r")
    for read in sam.fetch():
        i += 1
        if j > n:
            break
        if read.is_unmapped:
            continue
        if read.mapping_quality < 20:
            continue
        j += 1
        if read.is_paired:
            k += 1
    sam.close()
    
    if i == 0:
        return('PE')
    if j == 0:
        return('PE')
    if k / j < f:
        return('SE')
    else:
        return('PE')
