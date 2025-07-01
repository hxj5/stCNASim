# config.py - configuration


import sys
from ..config import Defaults as MainDefaults



COMMAND = "afc"

class Config:
    """Configuration of the `afc` (allele-specific counting) module.

    Attributes
    ----------
    See :func:`~.main.afc_wrapper()`.
    """
    def __init__(self):
        self.defaults = Defaults()

        # command-line arguments/parameters.
        self.sam_fn = None
        self.barcode_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.out_dir = None
        self.debug = self.defaults.DEBUG

        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.ncores = self.defaults.NCORES
        
        # snp filtering.
        self.min_count = self.defaults.MIN_COUNT
        self.min_maf = self.defaults.MIN_MAF
                
        # read assignment.
        self.strandness = self.defaults.STRANDNESS
        self.min_include = self.defaults.MIN_INCLUDE
        self.multi_mapper_how = self.defaults.MULTI_MAPPER_HOW
        self.xf_tag = self.defaults.XF_TAG
        self.gene_tag = self.defaults.GENE_TAG
        
        # read filtering.
        self.min_mapq = self.defaults.MIN_MAPQ
        self.min_len = self.defaults.MIN_LEN
        self.incl_flag = self.defaults.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.NO_ORPHAN
        
        # out_feature_dirs : list of str
        #   A list of output folders for feature-specific results.
        self.out_feature_dirs = None


        # internal parameters.

        # alleles : tuple of str
        #   All alleles.
        self.alleles = ("A", "B", "D", "O", "U")

        # cumi_alleles : tuple of str
        #   Alleles whose CUMIs will be outputed for read sampling.
        self.cumi_alleles = ("A", "B", "U")
        
        # hap_idx_tag : int
        #   Tag for haplotype index.
        self.hap_idx_tag = "HI"

        # out_prefix : str
        #   The prefix of the output files.
        self.out_prefix = COMMAND
        
        # pe_mode : str
        #   Mode of pair-end reads.
        #   One of 
        #   * 'auto' - auto detect;
        #   * 'PE' - pair-end reads;
        #   * 'SE' - single-end reads;
        self.pe_mode = 'auto'
        
        # no_orphan_post_qc : bool
        #   Whether to filter post-QC orphan reads when determining the 
        #   haplotype states of the UMIs.
        self.no_orphan_post_qc = self.defaults.NO_ORPHAN

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%sbarcode_file = %s\n" % (prefix, self.barcode_fn)
        s += "%sfeature_file = %s\n" % (prefix, self.feature_fn)
        s += "%ssnp_file = %s\n" % (prefix, self.snp_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sdebug_level = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%snumber_of_cores = %d\n" % (prefix, self.ncores)
        s += "%s\n" % prefix
        
        # snp filtering.
        s += "%smin_count = %d\n" % (prefix, self.min_count)
        s += "%smin_maf = %f\n" % (prefix, self.min_maf)
        s += "%s\n" % prefix

        # read assignment.
        s += "%sstrandness = %s\n" % (prefix, self.strandness)
        s += "%smin_include = %f\n" % (prefix, self.min_include)
        s += "%smulti_mapper_how = %s\n" % (prefix, self.multi_mapper_how)
        s += "%sxf_tag = %s\n" % (prefix, str(self.xf_tag))
        s += "%sgene_tag = %s\n" % (prefix, str(self.gene_tag))
        s += "%s\n" % prefix

        # read filtering.
        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix
        
        s += "%slen(out_feature_dirs) = %d\n" % (prefix, \
            len(self.out_feature_dirs) if self.out_feature_dirs else 0)
        s += "%s\n" % prefix

        
        # internal parameters.

        s += "%salleles = %s\n" % (prefix, str(self.alleles))
        s += "%scumi_alleles = %s\n" % (prefix, str(self.cumi_alleles))
        s += "%shap_idx_tag = %s\n" % (prefix, self.hap_idx_tag)
        s += "%sout_prefix = %s\n" % (prefix, self.out_prefix)
        s += "%spe_mode = %s\n" % (prefix, self.pe_mode)
        s += "%sno_orphan_post_qc = %s\n" % (prefix, self.no_orphan_post_qc)
        s += "%s\n" % prefix

        fp.write(s)
        
        
    def is_stranded(self):
        return self.strandness in ("forward", "reverse")

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None
    
    def use_xf(self):
        return self.xf_tag is not None
    
    def use_gene(self):
        return self.gene_tag is not None



class Defaults(MainDefaults):
    def __init__(self):
        super().__init__()



if __name__ == "__main__":
    conf = Config()
    conf.show()
