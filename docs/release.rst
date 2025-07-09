
..
   History
   =======
   
   
Release v0.6.0 (09/07/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version mainly restructure the ``mcount`` object in ``afc`` module.
Specifically, it stores intermediate read information in a hierarchical
structure, from sread (single read), gread (group of reads, compatible with
both SE and PE reads), to ucnt (all greads of one UMI).
The haplotype information is firstly inferred in sread level, based on the
phased SNPs covered by the sread, and then passed to higher levels.

The ``afc`` workflow for one feature is:

* fetch reads for every phased SNPs covered by this feature, and store the
  SNP info into the sread objects.
* iterate all reads of this feature and assign the reads into sread, gread,
  and ucnt objects.
* reads QC, e.g., removing orphan reads if ``no_orphan_post_qc`` is True.
* infer UMI-level allele (base) for every SNP of this feature, and do SNP QC.
  If QC failed, then remove the SNP from sread, gread, and ucnt objects.
* infer haplotype info of sread, gread, and ucnt objects based on the post-QC
  phased SNPs.

Compared to previous versions, this new version:

1. when determing the UMI-level allele (base) for one SNP, it checks whether
   the the allele (base) is identical in every covered gread (UMI collapsing).

2. it explicitly checks whether the haplotype information is homogeneous 
   within a group in various levels, i.e.,

* if multiple SNPs are located in one sread, whether these SNPs give same
  haplotype information for this read.
* when gread contains multiple sreads (e.g., for PE reads), whether these
  sreads have identical haolotype info.
* when UMI contains multiple greads, whether these greads have identical 
  haolotype info.
  
Previously, homogeneity of haplotype info is only checked in UMI level, by
checking the phased SNPs covered by this UMI.


Others:

* afc: add option ``no_orphan_post_qc`` to control whether to remove the
  orphan read if its mate read is missing or filtered.
* cs: add option ``barcode_whitelist_fn`` to enable inputting candidate 
  cell barcodes for simulated data.
* change default ``min_include`` from 0.9 to 0.5.
* rename ``cn_fold`` to ``cn_ratio``.
   


Release v0.5.3 (16/06/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* cs: change value of ``min_nonzero_num`` from single int to a tuple of int,
  default (1, 1, 3), 
  which is the minimum number of cells that have non-zeros in one feature,
  for alleles 'A', 'B', and 'U', respectively.
* cs: update cell-wise QC parameters, i.e., ``qc_min_features`` from 0.01 
  (fraction) to 100 (count); ``qc_cw_low_quantile`` from 0.005 to 0.0; and
  ``qc_cw_up_quantile`` from 0.995 to 1.0.



Release v0.5.2 (03/06/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* add option ``libsize_ratio`` (default 1.0) which is the ratio of 
  library size of simulated cells compared to seed cells.
  This option enables users to simulate counts with various overall coverage,
  or total library size.
   


Release v0.5.1 (22/05/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version produces the same results as v0.5.0, while it has a few minor
updates, including:

* config: fix typo that default ``min_count = 1`` and ``min_maf = 0``.
* update some doc files.



Release v0.5.0 (13/05/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version addresses the issue of multi-gene reads.
Here, the multi-gene reads are reads that

* map uniquely to single genomic locus where two or more genes overlap.
* map to multiple genomic loci, with each locus annotated to a different gene.
  
We wrap the pipeline for processing multi-gene reads into ``afc.mfu`` 
(multi-feature UMI) submodule, and add one option ``multi_mapper_how`` to 
select which strategy to use.

Specifically, in the ``afc.mfu`` submodule, we merge all ``afc`` output CUMIs 
(combination of cell and UMI barcodes) from all genes into one file, 
and define multi-gene UMIs as those shared by different genes or different 
alleles of the same gene.
By default, the multi-gene UMIs will be discarded 
(multi_mapper_how = "discard").
Furthermore, the allele-specific count matrices will be re-calculated based on
the updated list of CUMIs.

Additionally, we add two new options ``xf_tag`` and ``gene_tag`` to make full
use of the optimized read assignment from CellRanger or SpaceRanger, while 
the two tools consider at least gene structures (intron and exon etc) and 
the multi-gene reads.
These two options are useful for making the feature-counting results close to
the CellRanger or STAR counting, by providing almost-the-same input (reads).


Implementation:

* cs: change default ``min_nonzero_num`` from 5 to 1.
* pp: check duplicates in each input file.
* rename option ``merge_features_how`` to ``overlap_features_how``;
  rename its value "none" to "raw" and set as default.

Others:

* afc: reduce number of batches.
* docs: update manual.



Release v0.4.1 (16/04/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version improves its computational efficiency in terms of running speed
and memory usage, mainly in the cs and rs modules.
Specifically,

* try to use sparse matrix instead of numpy ndarray to store counts.
* cs: improve batch assignment in fit_RD() and simu_RD() to make full use of 
  multiprocessing pool.
* rs: optimize ``max_mem`` to 4G in pysam.sort to speedup.
  Now the peak memory of the whole framework is roughly ``4G * ncores`` for
  scRNA-seq data of typical size.
* utils.xbarcode: use standard random.sample() instead of numpy.random.choice().
* improve multiprocessing memory footprint via gc.collect().

Implementation:

* main: add interface to ``loss_allele_freq``.
  Note that setting the ``loss_allele_freq`` from 0.01 to 0.001 leads to
  failure of detecting LOHs by Numbat, possibly because of allelic signal loss
  resulting from the homozygous SNP filtering in Numbat.
* cs: add interface to ``cna_mode``.
* set default ``min_count = 20`` and ``min_maf = 0.1``.
* cs: restructure the function for allele-specific count simulation, to make
  it more modular.
* minor restructure all four modules, wrapping the steps of loading data and
  checking arguments into xxx_init() or xxx_pp().

Bug fix:

* replace loc with iloc in DataFrame.
* utils: use samtools cat when pysam.cat raises error as some version of pysam
  does not recognize "-@/--threads" option.
* afc: mark one issue of anndata that it cannot save adata whose ".X" is None.

Others:

* delete deprecated source codes.



Release v0.4.0 (30/03/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version restructure the afc and rs modules by using the feature (gene)
as unit.
Specifically,

* afc: output the allele-specific alignments for each feature.
* rs: read sampling is performed in each feature instead of in each chromosome.
* all feature-specific data and results are stored inside feature-specific
  folders.

Implementation:

* rs: only use post-filtering SNPs (by min_count and min_maf) for read masking.
* main: add interface to ``min_count`` and ``min_maf``.

Others:

* improve the batch assignment to take full advantage of multiprocessing pool,
  e.g., when ncores=10, split features into more batches (e.g., 150) than the
  usual 10 batches.

Warning:

* significant API changes: quite a few functions, classes, and modules have
  been renamed.



Release v0.3.0 (18/03/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version takes strandness into account.
Both CellRanger and STARsolo account for strandness, by default they only
count the sense reads.

* update the input feature annotation file, adding the fifth column "strand",
  in which "+" and "-" stand for the positive and negative strand of the 
  feature, respectively.
* the strand information is considered in resolving gene overlaps.
  If strandness is "forward" or "reverse" (i.e., strand-specific data), then
  the gene overlaps have to be on the same strand;
  otherwise, the gene overlaps are classified purely based on their genomic
  range, regardless of which strands they are in.
* add option ``strandness`` for specifying the strandness of the sequencing
  protocol.
  Three possible values are
  (1) "forward" (default): read strand is same as source RNA; 
  (2) "reverse": read strand is opposite to the source RNA;
  (3) "unstranded": the protocol is not strand-specific, i.e., read strand
  could be same as or opposite to its source RNA.
* strand information is considered in read filtering and assignment.
  (1) when strandness=forward, for single end (SE) read, the read has to be
  sense read; for pair end (PE) reads, R1 has to be antisense and R2 sense.
  (2) when strandness=reverse, the rules in (1) are reversed.
  (3) when strandness=unstranded, no check on sense or antisense.
  

Implementation:

* discard the reads with 'N' in its UMI.
* clean unused strategies for resolving gene overlaps.
  Only keep "quantile" (alias to "quantile2"), "union", and "none".

Output:

* update output BAM tags, mainly make the CR and UR values match the newly
  generated CB and UB values.
  It is useful for STARsolo feature counting because STARsolo requires 
  specifying tags of raw cell and UMI barcodes (default CR and UR) for UMI
  grouping (collapsing).
* compress output h5ad files with option gzip.
  It can greatly reduce the file size.

Others:

* cs: suppress statsmodels RuntimeWarning messages.
* detect whether the BAM index file exists when using pysam fetch().
  If not exist, the simulator will report error.
  For some (or all?) pysam version, while BAM index is required for using
  fetch() method, it does not report error when the index file is missing.
* simplify docstrings of Config.
* docs: add ref genome into section "input" in manual.
* README: add potential issues of installation related to pysam installation.
  pysam can be installed via conda when pip install failed.



Release v0.2.0 (05/03/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Improve the quality of simulated cells, to avoid generating some noisy clones.

* cs: add QC step to filter low-quality seed cells, e.g., 
  with small library size or small number of expressed features.
  The filtered cells are outputted for potential further analysis.
* cs: use more stringent up and low bound of simulated library size, e.g.,
  the minimum simulated library size allowed is 1000.

Input:

* update file format of CNA profile, removing the ``region`` column.

For library size simulation:

* cs: add lognormal and swr (sampling with replacement) strategies for
  library size fitting and simulation.
* cs: add interface for default kwargs_fit_sf and kwargs_fit_rd.
  Set ``lognormal`` as default strategy for library size (size factor)
  fitting and simulation.

For fitting read depth:

* cs: use Poisson as default when distribution fitting is not converged.
* cs: set default max_iter to 1000 when fitting read depth.

Others:

* cs: add small epsilon value to mean when calculating cv.
* mark module or folder ``tests`` deprecated.
* pp: rename the filename of features after resolving overlapping features.
  Specifically, suffix changed from "merged.tsv" to "resolve_overlap.tsv".
* better support processing SNP file names, no matter the suffix is in
  lower or upper case.


Release v0.1.2 (11/02/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version mainly aims to reduce the variation in the simulated BAF signals
of normal features/regions.

* afc: set default min_count=20, min_maf=0.1.
  It may filter some input phased SNPs whose expression levels in the seed
  data are low.
  Motivation: if one gene contains mainly lowly-expressed SNPs, then its
  haplotype-specific counts (Hap-A and Hap-B) will be small, and its simulated
  Hap-A and Hap-B counts are probably also small, hence AF may be biased
  towards 0 or 1.
* cs: use Poisson distribution when fitting NB failed.
  Previously, empirical parameters of NB were used when fitting NB failed.
  We expect Poisson to produce lower variation level in simulated counts, 
  compared to NB, especially for lowly-expressed features.
* pp: set "quantile2" as default option of ``merge_features_how``.
  Both "quantile2" and "quantile2_union" strategies can remove features that
  overlap large number of other features, while "quantile2" seem to produce
  stronger CNA signals.
* main: logging APP and VERSION.


Release v0.1.1 (03/02/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* pp: add ``merge_features_how`` - How to merge overlapping features.
* Support both INT and FLOAT as value of ``--minINCLUDE``.
  If float between (0, 1), it is the minimum fraction of included length.
* Set default value of ``--minINCLUDE`` or ``min_include`` as 0.9.
* docs: add TODO.


Release v0.1.0 (06/12/2024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Add ``--minINCLUDE`` option for read filtering.

* ``--minINCLUDE`` is the minimum length of included part within specific
  feature. 
* For example, if the genomic range of a feature is chr1:1000-3000, and one
  fetched read (100bp) aligned to two locus, chr1:601-660 (60bp) and 
  chr1:3801-3840 (40bp), then no any part of the read is actually included 
  within the feature, hence it will be filtered by ``--minINCLUDE=30``, 
  whereas older versions of scCNASim may keep the read.
  Note, when features are processed independently, one read filtered by
  --minINCLUDE in one feature may still be fetched and counted by other 
  features.
* Previously, there is noise present in inferCNV heatmap that both signals 
  of duplication and deletion present in a strip of genes, even in the
  reference cells.
  By using ``--minINCLUDE`` (default 30), the noise is largely removed.
  
Others

* rs: do not output sampled reads of multi-feature UMIs for non-overlapping
  features.
  If one multi-read UMI is sampled by specific feature (in rs module), and
  some of its reads are not included within the feature (``--minINCLUDE``),
  then those reads will not be outputted to BAM for this feature.
  Without this step, there will be inflation of UMI counts in rs BAM, compared
  to the simulated counts in cs module, considering the non-included reads may
  be counted by other features.
* rs: output sampled UMIs aligned to distinct alleles in different features.
  Assume there is a multi-feature UMI (due to error in UMI collapse?) 
  aligned to distinct alleles in different features, e.g., Hap-B in one 
  feature and Hap-U in another feature.
  If the UMI is sampled by both features, then the UMI is outputted for both
  features, while mimicking the real scRNA-seq BAM (error in UMI collapse?).
  Previously, this UMI is only outputted once for one (first iterated) 
  feature, which may result in the decrease of UMI counts in rs BAM, compared
  to the simulated counts in cs module.
* pp: filter features by chromosomes.
  Filter features whose chromosomes are not in the input chrom list.
* convert column chrom astype str in anndata.
  Previously, the chrom column will be of int dtype if all chromosome names are
  numeric strings, e.g., "1", "2", etc.
* init setting random seed.
  Currently the whole simulation results are not reproducible with a seed,
  possibly due to the parallel computing.
* cs: also output the counts into sparse matrices, in addition to the
  ``h5ad`` file.
* pp and afc: rename ``utils`` to ``io``.


Bug fix:

* utils: fix bug in ``xbarcode.str2int()``.


Release v0.0.2 (12/10/2024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* rename CNV to CNA.
* allow input empty CNA profile file.
* require Python>=3.11.
* fix typos.


Release v0.0.1 (17/09/2024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Implement a pipeline wrapping four modules:

#. ``pp``: preprocessing.
#. ``afc``: allele-specific feature counting.
#. ``cs``: count simulation.
#. ``rs``: read simulation.
