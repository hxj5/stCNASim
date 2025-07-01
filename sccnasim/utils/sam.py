# sam.py - sam alignment processing.


import gc
import multiprocessing
import os
import pysam
import shutil
import subprocess
from logging import info, error



def check_read(read, reg, conf):
    if conf.xf_tag:
        if not read.has_tag(conf.xf_tag):
            return(-101)
        xf = read.get_tag(conf.xf_tag)
        if xf not in (17, 25):
            return(-102)
    if conf.gene_tag:
        if not read.has_tag(conf.gene_tag):
            return(-106)
        gene = read.get_tag(conf.gene_tag)
        if gene != reg.name:
            return(-107)
    if check_strand(read, reg.strand, conf.strandness) < 0:
        return(-111)
    if check_included(read, reg.start, reg.end, conf.min_include) < 0:
        return(-112)
    ret = check_basic(read, conf)
    return(ret)



def check_basic(read, conf):
    """Basic QC of a read.

    This function checks whether a read is valid.
    If invalid, it will be filtered.
    
    Parameters
    ----------
    read : pysam.AlignedSegment
        One alignment read.
    conf : object
        Configuration object whose attributes will be used as filtering
        criterias.
        
    Returns
    -------
    int
        Return code. 0 if read is valid, negative otherwise.
    """
    if read.is_unmapped:
        return(-1)
    if read.mapq < conf.min_mapq:
        return(-2)
    if conf.excl_flag and read.flag & conf.excl_flag:
        return(-3)
    if conf.incl_flag and not read.flag & conf.incl_flag:
        return(-4)
    if conf.no_orphan and read.flag & BAM_FPAIRED and not \
        read.flag & BAM_FPROPER_PAIR:
        return(-5)
    if conf.cell_tag and not read.has_tag(conf.cell_tag):
        return(-11)
    if conf.umi_tag:
        if read.has_tag(conf.umi_tag):
            umi = read.get_tag(conf.umi_tag)
            if umi and 'N' in umi.upper():
                return(-13)
        else:
            return(-12)
    if len(read.positions) < conf.min_len:
        return(-21)
    return(0)



def check_strand(read, feature_strand, strandness = "forward"):
    """Check whether the strand of a read is valid.
    
    Different from the rules of htseq-count `--stranded`, since 10x Genomics
    scRNA-seq platform strandness is "forward", and
    - SE read should be sense; 
    - PE R1 should be antisense and R2 sense.
    
    Parameters
    ----------
    read : pysam.AlignedSegment
        One alignment read.
    feature_strand : str
        DNA strand orientation of the feature, "+" (positive) or 
        "-" (negative).
    strandness : {"forward", "reverse", "unstranded"}
        Strandness of the sequencing protocol.
        "forward" - read strand same as the source RNA molecule;
        "reverse" - read strand opposite to the source RNA molecule;
        "unstranded" - no strand information.
        
    Returns
    -------
    int
        Return code. 0 if read is valid, negative otherwise.
    """
    def __get_expected(read, se, pe_r1, pe_r2):
        if read.is_paired:
            if read.is_read1:
                return(pe_r1)
            else:
                return(pe_r2)
        else:
            return(se)


    if strandness not in ("forward", "reverse"):
        return(0)
    
    read_strand = "+" if read.is_forward else "-"
    sense = feature_strand
    anti_sense = "+" if sense == "-" else "-"
    if strandness == "forward":
        expected_strand = __get_expected(read, sense, anti_sense, sense)
        return(0 if read_strand == expected_strand else -3)
    else:
        expected_strand = __get_expected(read, anti_sense, sense, anti_sense)
        return(0 if read_strand == expected_strand else -5)

    
    
def check_included(read, start, end, min_include):
    """Check whether a read is included within specific feature.
    
    Parameters
    ----------
    read : pysam.AlignedSegment
        One alignment read.
    start : int
        The start genomic position of the feature, 1-based and inclusive.
    end : int
        The end genomic position of the feature, 1-based and exclusive.    
    min_include : int or float
        Minimum length of included part within specific feature.
        
    Returns
    -------
    int
        Return code. 0 if read is included, negative otherwise.
    """
    if 0 < min_include < 1:
        if get_include_frac(read, start, end) < min_include:
            return(-3)
    else:
        if get_include_len(read, start, end) < min_include:
            return(-5)
    return(0)



def get_include_frac(read, s, e):
    """Get the fraction of included part within specific feature.

    Parameters
    ----------
    read : pysam.AlignedSegment
        One alignment read.
    s : int
        The start genomic position of the feature, 1-based and inclusive.
    e : int
        The end genomic position of the feature, 1-based and exclusive.

    Returns
    -------
    float or None
        The fraction of included part within specific feature.
        None if no any part of the read is aligned.
    """
    n = len(read.positions)
    if n <= 0:
        return(None)
    m = get_include_len(read, s, e)
    return(m / float(n))



def get_include_len(read, s, e):
    """Get the length of included part within specific feature.

    Parameters
    ----------
    read : pysam.AlignedSegment
        One alignment read.
    s : int
        The start genomic position of the feature, 1-based and inclusive.
    e : int
        The end genomic position of the feature, 1-based and exclusive.

    Returns
    -------
    int
        The length of included part within specific feature.
    """
    include_pos_list = [x for x in read.positions if s - 1 <= x <= e - 2]
    return(len(include_pos_list))



def get_query_bases(read, full_length = False):
    """Qurey bases that are within the alignment.

    Parameters
    ----------
    read : pysam.AlignedSegment
        One alignment read.
    full_length : bool, default False
        If full_length is True, `None` values will be included for any
        soft-clipped or unaligned positions within the read. 
        The returned list will thus be of the same length as the `read`.
    
    Returns
    -------
    list of str
        A list of bases in qurey sequence that are within the alignment.
    """
    cigar_tuples = read.cigartuples
    if not cigar_tuples:
        return []

    result = []
    pos = 0
    s = read.query_sequence

    for op, l in cigar_tuples:
        if op == BAM_CSOFT_CLIP or op == BAM_CINS:
            if full_length:
                for i in range(0, l):
                    result.append(None)
            pos += l
        elif op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i in range(pos, pos + l):
                result.append(s[i])
            pos += l
        # else: do nothing.
    return result



def get_query_qualities(read, full_length = False):
    """Qurey qualities that are within the alignment.

    Parameters
    ----------
    read : pysam.AlignedSegment
        One alignment read.
    full_length : bool, default False
        If full_length is True, `None` values will be included for any 
        soft-clipped or unaligned positions within the read. 
        The returned list will thus be of the same length as the `read`.
    
    Returns
    -------
    list of int
        A list of qualities of bases that are within the alignment.
        Note that the returned qual values are not ASCII-encoded values 
        typically seen in FASTQ or SAM formatted files, no need to 
        substract 33.
    """
    cigar_tuples = read.cigartuples
    if not cigar_tuples:
        return []

    result = []
    pos = 0
    s = read.query_qualities

    for op, l in cigar_tuples:
        if op == BAM_CSOFT_CLIP or op == BAM_CINS:
            if full_length:
                for i in range(0, l):
                    result.append(None)
            pos += l
        elif op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i in range(pos, pos + l):
                result.append(s[i])
            pos += l
        # else: do nothing.
    return result



BAM_FPAIRED = 1
BAM_FPROPER_PAIR = 2

# Cigar
# reference: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9
