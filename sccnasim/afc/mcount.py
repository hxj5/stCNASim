# mcount.py



from logging import warning as warn
from ..utils.sam import get_query_bases



class SRead:
    """One single Read."""
    
    # - reads info is stored in `snp_data`.
    # - after all reads of the feature are added, the `infer_xxx()` can be
    #   called.
    
    def __init__(self, rid):
        # rid : str
        #   ID of the read. 
        #   Typically "{chrom}:{start}-{end}".
        #   Note that PE read could have different IDs.
        self.rid = rid
        
        # snp_data : dict of {str : tuple}
        #   Data of SNPs covered by this read.
        #   Key is the SNP ID, value is a tuple of (int, str), which is the 
        #   hap_idx and allele of SNP in this read.
        self.snp_data = None
        
        
        # Below are inferred attributes.
        
        # hap_flag : int
        #   The haplotype flag of this read. 
        #   Each bit:
        #   * 1 (A): some SNP is in SNP-hap status 'A'.
        #   * 1<<1 (B): some SNP is in SNP-hap status 'B'.
        #   * 1<<2 (O): some SNP is in SNP-hap status 'O'.
        self.hap_flag = None

        # hap_idx : {-3, -2, -1, 0, 1, 2}
        #   The haplotype index of this read.
        #   * 0 (A): some SNPs are in SNP-hap status 'A' but no SNPs are in
        #     SNP-hap status 'B'.
        #   * 1 (B): some SNPs are in SNP-hap status 'B' but no SNPs are in
        #     SNP-hap status 'A'.
        #   * 2 (D; both): some SNPs are in SNP-hap status 'A' and some SNPs
        #     are in SNP-hap status 'B'.
        #   * -1 (O; oth): no SNPs are in SNP-hap status 'A' or 'B', but some
        #     SNPs are in SNP-hap status 'O'.
        #   * -2 (U; unknown): all SNPs are in SNP-hap status 'U'.
        #   * -3 (U; unknown): the read is not fetched by any SNPs.
        self.hap_idx = None

        
    def add_read(self, read):
        return(0)

    
    def add_read_from_snp(self, read, snp_cnt):
        rid = get_read_id(read)
        assert rid == self.rid
        if self.snp_data is None:
            self.snp_data = dict()
        hap_idx, allele = get_hap_idx(read, snp_cnt.snp)
        snp_id = snp_cnt.get_id()
        if snp_id in self.snp_data:
            #warn("read '%s' (%s) fetched more than once." % \
            #    (read.query_name, rid))
            pass
        else:
            self.snp_data[snp_id] = (hap_idx, allele)
        return(0)
    
    
    def del_snp(self, snp_id):
        if self.snp_data is not None and snp_id in self.snp_data:
            del self.snp_data[snp_id]
    
    
    def __infer_hap_idx_of_sread(self, x):
        """x: haplotype indexes of SNPs."""
        n = {h:0 for h in (-2, -1, 0, 1)}
        for h in x:
            assert h in n
            n[h] += 1
            
        flag = 0
        if n[0] > 0:
            flag |= (1<<0)
        if n[1] > 0:
            flag |= (1<<1)
        if n[-1] > 0:
            flag |= (1<<2)
        
        hap_idx = None
        if n[0] > 0:
            if n[1] > 0:
                hap_idx = 2
            else:
                hap_idx = 0
        elif n[1] > 0:
            hap_idx = 1
        elif n[-1] > 0:
            hap_idx = -1
        else:
            hap_idx = -2
        return(hap_idx, flag, n)
    
    
    def infer_haplotype(self):
        """This function infers the haplotype of the read.
        
        It infers the haplotype state of one read by checking all the phased
        SNPs covered by this read, since each covered SNP carries some 
        information about the haplotype state of the read.
        """
        self.hap_flag = 0
        if self.snp_data is None or len(self.snp_data) <= 0:
            self.hap_idx = -3
            return
        snp_hidx = [v[0] for v in self.snp_data.values()]
        self.hap_idx, self.hap_flag, _ = self.__infer_hap_idx_of_sread(snp_hidx)



class UMICount:
    # - reads info is stored in `greads`.
    # - after all reads of the feature are added, the `infer_xxx()` can be
    #   called.
    
    def __init__(self, umi, conf):
        self.umi = umi
        self.conf = conf
        
        # greads : dict
        #   Two layers of dict storing `SRead` objects belonging to this UMI.
        #   The first layer: key is the read QNAME shared by this group, 
        #   values are dict;
        #   the second layer: key is read id, i.e., "chrom:start-end", and 
        #   the value is the `SRead` object.
        self.greads = dict()
        

        # Below are inferred attributes.
        
        # hap_flag : int
        #   The haplotype flag of this UMI. 
        #   Each bit:
        #   * 1 (A): some gread is in gread-hap status 'A'.
        #   * 1<<1 (B): some gread is in gread-hap status 'B'.
        #   * 1<<2 (O): some gread is in gread-hap status 'O'.
        self.hap_flag = None
        
        # hap_idx : int
        #   The haplotype index of this UMI.
        #   - A (Haplotype-A; internal index: 0)
        #       No greads are in gread-hap status 'B' or 'D' &&
        #       some greads are in gread-hap status 'A'.
        #   - B (Haplotype-B; internal index: 1)
        #       No greads are in gread-hap status 'A' or 'D' &&
        #       some greads are in gread-hap status 'B'.
        #   - D (Duplicate; internal index: 2)
        #       Some greads are in gread-hap status 'D' ||
        #       some greads are in gread-hap status 'A' and some greads are
        #       in gread-hap status 'B'.
        #   - O (Others; internal index: -1)
        #       No greads are in gread-hap status 'A', 'B' or 'D' &&
        #       some greads are in gread-hap status 'O'.
        #   - U (Unknown; internal index: -2)
        #       All greads are in gread-hap status 'U'.
        self.hap_idx = None
        
        # snp_data : dict of {str : str}
        #   Data of SNPs covered by this UMI.
        #   Key is the SNP ID, value is the **inferred** allele of the SNP
        #   given the SNP alleles in all greads in this UMI.
        self.snp_data = None

        
    def add_read(self, read):
        gid = get_read_gid(read)
        if gid not in self.greads:
            self.greads[gid] = dict()
        rid = get_read_id(read)
        if rid not in self.greads[gid]:
            self.greads[gid][rid] = SRead(rid)
        sr = self.greads[gid][rid]
        ret = sr.add_read(read)
        return(ret)
        
        
    def add_read_from_snp(self, read, snp_cnt):
        gid = get_read_gid(read)
        if gid not in self.greads:
            self.greads[gid] = dict()
        rid = get_read_id(read)
        if rid not in self.greads[gid]:
            self.greads[gid][rid] = SRead(rid)
        sr = self.greads[gid][rid]
        ret = sr.add_read_from_snp(read, snp_cnt)
        return(ret)
    
    
    def del_snp(self, snp_id):
        for gid, gr_dat in self.greads.items():
            for rid, sr in gr_dat.items():
                sr.del_snp(snp_id)
    
    
    def __infer_allele_of_gread(self, x):
        """x: alleles (bases) of one SNP in several sreads. 
        It should not contain None.
        """
        if len(x) <= 0:
            return(None)
        y = set(x)
        if len(y) > 1:
            return(None)
        return(x[0].upper())
    
    
    def __infer_allele_of_umi(self, x):
        n = dict()
        for a in x:
            if a is None:
                continue
            a = a.upper()
            assert a in ('A', 'C', 'G', 'T', 'N')
            if a not in n:
                n[a] = 0
            n[a] += 1
        if len(n) <= 0:
            return(None)
        t = sum(n.values())
        m = max(n.values())
        if m / t < 0.9:
            return(None)
        for a, k in n.items():
            if k == m:
                return(a.upper())

    
    def __infer_hap_idx_of_gread(self, x):
        """x: haplotype indexes of sreads."""
        # hap_idx : int
        #   The haplotype index of this gread.
        #   - A (Haplotype-A; internal index: 0)
        #       No sreads are in sread-hap status 'B' or 'D' &&
        #       some sreads are in sread-hap status 'A'.
        #   - B (Haplotype-B; internal index: 1)
        #       No sreads are in sread-hap status 'A' or 'D' &&
        #       some sreads are in sread-hap status 'B'.
        #   - D (Duplicate; internal index: 2)
        #       Some sreads are in sread-hap status 'D' ||
        #       some sreads are in sread-hap status 'A' and some sreads are
        #       in sread-hap status 'B'.
        #   - O (Others; internal index: -1)
        #       No sreads are in sread-hap status 'A', 'B' or 'D' &&
        #       some sreads are in sread-hap status 'O'.
        #   - U (Unknown; internal index: -2)
        #       All sreads are in sread-hap status 'U'.
        
        n = {h:0 for h in (-3, -2, -1, 0, 1, 2)}
        for h in x:
            assert h in n
            n[h] += 1
            
        flag = 0
        if n[0] > 0:
            flag |= (1<<0)
        if n[1] > 0:
            flag |= (1<<1)
        if n[-1] > 0:
            flag |= (1<<2)
        
        hap_idx = None
        if n[2] > 0 or (n[0] > 0 and n[1] > 0):
            hap_idx = 2
        elif n[0] > 0:
            hap_idx = 0
        elif n[1] > 0:
            hap_idx = 1
        elif n[-1] > 0:
            hap_idx = -1
        else:
            hap_idx = -2
        return(hap_idx, flag, n)
    
    
    def __infer_hap_idx_of_umi(self, x):
        return self.__infer_hap_idx_of_gread(x)
    
    
    def infer_haplotype(self):
        """Infer haplotype of the UMI based on the haplotypes of its reads.
        """
        self.qc()
        hap_n_greads = {h:0 for h in (0, 1, 2, -1, -2)}     # number of supporting greads for each haplotype.
        gr_hidx = []
        for gid, gr_dat in self.greads.items():
            sr_hidx = []
            for rid, sr in gr_dat.items():
                sr.infer_haplotype()
                sr_hidx.append(sr.hap_idx)
            hap_idx, flag, _ = self.__infer_hap_idx_of_gread(sr_hidx)
            gr_hidx.append(hap_idx)
        self.hap_idx, self.hap_flag, _ = self.__infer_hap_idx_of_umi(gr_hidx)
        
     
    def infer_allele(self):
        """Infer the allele of each SNP covered by this UMI."""
        self.qc()
        gr_ale_data = dict()             # {snp_id : a list of gread-level alleles of that SNP}
        for gid, gr_dat in self.greads.items():
            sr_ale_data = dict()         # {snp_id : a list of sread-level alleles of that SNP}
            for rid, sr in gr_dat.items():
                if sr.snp_data is None:
                    continue
                for snp_id, (hap_idx, sr_ale) in sr.snp_data.items():
                    if sr_ale is not None:
                        if snp_id not in sr_ale_data:
                            sr_ale_data[snp_id] = []
                        sr_ale_data[snp_id].append(sr_ale)
            for snp_id, sr_alleles in sr_ale_data.items():
                gr_ale = self.__infer_allele_of_gread(sr_alleles)
                if gr_ale is not None:
                    if snp_id not in gr_ale_data:
                        gr_ale_data[snp_id] = []
                    gr_ale_data[snp_id].append(gr_ale)
                    
        snp_data = dict()                # {snp_id : UMI-level allele of that SNP}
        for snp_id, gr_alleles in gr_ale_data.items():
            umi_ale = self.__infer_allele_of_umi(gr_alleles)
            if umi_ale is not None:
                snp_data[snp_id] = umi_ale
        self.snp_data = snp_data
            
            
    def qc(self):
        conf = self.conf
        qcfail_gids = []
        for gid, gr_dat in self.greads.items():
            if conf.pe_mode == "PE" and len(gr_dat) != 2:
                if conf.no_orphan_post_qc:
                    qcfail_gids.append(gid)
        for gid in qcfail_gids:            
            del self.greads[gid]


 
class MCount:
    def __init__(self, samples, conf):
        self.cells = samples
        self.conf = conf
        
        # umi_cnts : dict
        #   Two layers of dict storing `UCount` objects belonging to the
        #   feature.
        #   The first layer is the sample/cell level: where the key is the
        #   cell barcodes and values are dict;
        #   the second layer is the UMI level, where the key is the UMI
        #   barcode, and the value is the `UCount` object.
        self.umi_cnts = {cell:dict() for cell in self.cells}
        
        
    def add_read(self, read):
        """Count one fetched read given specific feature.
        
        Parameters
        ----------
        read : :class:`~pysam.AlignedSegment`
            A fetched BAM read of specific feature.

        Returns
        -------
        int
            Return code. 0 if success, -1 error, 1 QC fail.
        UMICount or None
            The UMI object that the `read` belongs to.
        """
        conf = self.conf
        cell, umi = get_cell_and_umi(read, conf)
        if (not cell) or (cell not in self.umi_cnts) or (not umi):
            return(1, None)
        if umi not in self.umi_cnts[cell]:
            self.umi_cnts[cell][umi] = UMICount(umi, conf)
        ucnt = self.umi_cnts[cell][umi]
        ret = ucnt.add_read(read)
        return(ret, ucnt)
        
        
    def add_read_from_snp(self, read, snp_cnt):
        """Add one read from SNPCount.
        
        Parameters
        ----------
        read : :class:`~pysam.AlignedSegment`
        cell : str
            Cell barcode.
        umi : str
            UMI barcode.
        snp_cnt : :class:`SNPCount`
        
        Returns
        -------
        int
            Return code. 0 if success, -1 error, 1 QC fail.
        UMICount or None
            The UMI object that the `read` belongs to.
        """
        conf = self.conf
        cell, umi = get_cell_and_umi(read, conf)
        if (not cell) or (cell not in self.umi_cnts) or (not umi):
            return(1, None)
        if umi not in self.umi_cnts[cell]:
            self.umi_cnts[cell][umi] = UMICount(umi, conf)
        ucnt = self.umi_cnts[cell][umi]
        ret = ucnt.add_read_from_snp(read, snp_cnt)
        return(ret, ucnt)
    
    
    def del_snp(self, snp_id):
        for cell, cdat in self.umi_cnts.items():
            for umi, ucnt in cdat.items():
                ucnt.del_snp(snp_id)
                
                
    def get_hap_idx(self, read):
        # this function should be called after `stat()` is called.
        cell, umi = get_cell_and_umi(read, self.conf)
        if cell not in self.umi_cnts:
            return(None)
        if umi not in self.umi_cnts[cell]:
            return(None)
        ucnt = self.umi_cnts[cell][umi]
        gid = get_read_gid(read)
        if gid not in ucnt.greads:
            return(None)
        return(ucnt.hap_idx)
    
    
    def stat(self):
        # stat_hap : dict
        #   Two layer of dict.
        #   First layer: key is hap_idx, value is dict.
        #   Second layer: key is cell barcode, value is number of UMIs.
        stat_hap = dict()
        
        # stat_umi : dict
        #   Two layer of dict.
        #   First layer: key is hap_idx, value is dict.
        #   Second layer: key is cell barcode, value is set of UMIs.
        stat_umi = dict()
        
        for cell, cdat in self.umi_cnts.items():
            for umi, ucnt in cdat.items():
                ucnt.infer_haplotype()
                assert ucnt.hap_idx is not None
                if ucnt.hap_idx not in stat_hap:
                    stat_hap[ucnt.hap_idx] = dict()
                    stat_umi[ucnt.hap_idx] = dict()
                if cell not in stat_hap[ucnt.hap_idx]:
                    stat_hap[ucnt.hap_idx][cell] = 0
                    stat_umi[ucnt.hap_idx][cell] = set()
                stat_hap[ucnt.hap_idx][cell] += 1
                stat_umi[ucnt.hap_idx][cell].add(umi)
        return(stat_hap, stat_umi)



class SNPCount:
    def __init__(self, snp, mcnt, samples, conf):
        # snp : :class:`~..utils.gfeature.SNP`
        self.snp = snp
        
        # mcnt : :class:`MCount`
        self.mcnt = mcnt
        
        self.cells = samples
        self.conf = conf

        # umi_cnts : dict
        #   Two layers of dict storing `UCount` objects covering this SNP.
        #   The first layer is the sample/cell level: where the key is the
        #   cell barcodes and values are dict;
        #   the second layer is the UMI level, where the key is the UMI
        #   barcode, and the value is the `UCount` object.
        self.umi_cnts = {cell:dict() for cell in self.cells}

        
    def add_read(self, read):
        """Count one fetched read covering the SNP.
        
        Parameters
        ----------
        read : :class:`~pysam.AlignedSegment`
            A fetched BAM read covering the SNP.

        Returns
        -------
        int
            Return code. 0 if success, -1 error, 1 QC fail. 
        """
        conf = self.conf
        cell, umi = get_cell_and_umi(read, conf)
        if (not cell) or (cell not in self.umi_cnts) or (not umi):
            return(1)
        ret, ucnt = self.mcnt.add_read_from_snp(read, self)
        if ret < 0:
            return(-1)
        elif ret > 0:
            return(1)
        if umi not in self.umi_cnts[cell]:
            self.umi_cnts[cell][umi] = ucnt
        return(0)
    
    
    def get_id(self):
        return self.snp.get_id()
    
    
    def stat(self):
        alleles = ('A', 'C', 'G', 'T', 'N')
        stat_allele = {a:0 for a in alleles}
        snp_id = self.get_id()
        for cell, cdat in self.umi_cnts.items():
            for umi, ucnt in cdat.items():
                ucnt.infer_allele()
                if snp_id not in ucnt.snp_data:
                    continue
                ale = ucnt.snp_data[snp_id]
                if ale is None:
                    continue
                assert ale in alleles
                stat_allele[ale] += 1
        return(stat_allele)



def get_cell_and_umi(read, conf):
    cell = umi = None
    if conf.use_barcodes():
        cell = read.get_tag(conf.cell_tag)
    else:
        raise ValueError

    if conf.use_umi():
        umi = read.get_tag(conf.umi_tag)
    else:
        raise ValueError
    return(cell, umi)



def get_hap_idx(read, snp):
    """Get haplotype index of the read and allele of SNP.
    
    It infers the haplotype state of the read by comparing the fetched SNP 
    allele to the phased ones, e.g., 
    if the fetched allele of this SNP is 'A' in this read, and the phased 
    (REF and ALT) alleles of the SNP are 'A' and 'C', respectively, 
    then the read would be inferred as from the REF haplotype.
    
    Returns
    -------
    int
        Haplotype of the `read`. 
        * 0 (A; ref): the fetched SNP allele is on the reference haplotype.
        * 1 (B; alt): the fetched SNP allele is on the alternative haplotype.
        * -1 (O; oth): some allele is fetched but is on neither the
          reference nor alternative haplotype.
        * -2 (U; unknown): no allele is extracted (allele is None), e.g.,
          in split read.
    str or None
        Allele (base) of the SNP.
    """
    allele = hap_idx = None
    try:
        idx = read.positions.index(snp.pos - 1)
    except:
        hap_idx = -2
    else:
        bases = get_query_bases(read, full_length = False)
        allele = bases[idx].upper()
        hap_idx = snp.get_hap_idx(allele)
    return(hap_idx, allele)
    

    
def get_read_id(read):
    start, end = read.positions[0], read.positions[-1]
    rid = "%s:%d-%d" % (read.reference_id, start, end)
    return(rid)



def get_read_gid(read):
    return read.query_name
