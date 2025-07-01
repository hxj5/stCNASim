
Manual
======

.. contents:: Contents
   :depth: 2
   :local:



Quick Usage
-----------
First, please look at `Input`_ to prepare the input data.

Then call the ``main_wrapper()`` function to run the simulation pipeline.

For typical single-cell or spatial transcriptomics data from 10x Genomics
platform, an example is:

.. code-block:: python

   from sccnasim import main_wrapper

   main_wrapper(
       sam_fn = "{sample}.bam",
       cell_anno_fn = "cell_anno.tsv", 
       feature_fn = "hg38.features.tsv",
       phased_snp_fn = "phased.snp.vcf.gz",
       clone_anno_fn = "clone_anno.tsv",
       cna_profile_fn = "cna_profile.tsv", 
       refseq_fn = "hg38.fa",
       out_dir = "./simu_result",
       umi_len = 10,
       strandness = "forward",      # set to "reverse" for typical 10x PE 5' scRNA-seq data.
       ncores = 10
   )

The full parameters can be found at section `Full Parameters`_.

You may also run each step (module) explicitly by calling corresponding 
wrapper functions (see ``Tutorial`` page).

See `Implementation`_ for details of the four modules.



Full Parameters
---------------

.. code-block:: python

    main_wrapper(
        sam_fn,
        cell_anno_fn, feature_fn, phased_snp_fn,
        clone_anno_fn, cna_profile_fn, 
        refseq_fn,
        out_dir,
        sam_list_fn = None, sample_ids = None, sample_id_fn = None,
        overlap_features_how = "raw",
        size_factor = "libsize",
        marginal = "auto",
        libsize_ratio = 1.0,
        loss_allele_freq = 0.01,
        kwargs_fit_sf = None,
        kwargs_fit_rd = None,
        chroms = "human_autosome",
        cell_tag = "CB", umi_tag = "UB", umi_len = 10,
        ncores = 1, seed = 123, verbose = False,
        min_count = 1, min_maf = 0,
        strandness = "forward", min_include = 0.5, multi_mapper_how = "discard",
        xf_tag = "xf", gene_tag = "GN",
        min_mapq = 20, min_len = 30,
        incl_flag = 0, excl_flag = -1,
        no_orphan = True,
        debug_level = 0
    )

    
The details are listed below:

sam_fn : str or None
    Comma separated indexed BAM file.
    Note that one and only one of `sam_fn` and `sam_list_fn` should be
    specified.

cell_anno_fn : str
    The cell annotation file. 
    It is a header-free TSV file and its first two columns are:
    
    - "cell" (str): cell barcodes.
    - "cell_type" (str): cell type.

feature_fn : str
    A TSV file listing target features. 
    It is header-free and its first 5 columns shoud be: 
    
    - "chrom" (str): chromosome name of the feature.
    - "start" (int): start genomic position of the feature, 1-based
      and inclusive.
    - "end" (int): end genomic position of the feature, 1-based and
      inclusive.
    - "feature" (str): feature name.
    - "strand" (str): feature strand, either "+" (positive) or 
      "-" (negative).

phased_snp_fn : str
    A TSV or VCF file listing phased SNPs.
    If TSV, it is a header-free file containing SNP annotations, whose
    first six columns should be:
    
    - "chrom" (str): chromosome name of the SNP.
    - "pos" (int): genomic position of the SNP, 1-based.
    - "ref" (str): the reference allele of the SNP.
    - "alt" (str): the alternative allele of the SNP.
    - "ref_hap" (int): the haplotype index of `ref`, one of {0, 1}.
    - "alt_hap" (int): the haplotype index of `alt`, one of {1, 0}.
    
    If VCF, it should contain "GT" in its "FORMAT" field.

clone_anno_fn : str
    A TSV file listing clonal anno information.
    It is header-free and its first 3 columns are:
    
    - "clone" (str): clone ID.
    - "source_cell_type" (str): the source cell type of `clone`.
    - "n_cell" (int): number of cells in the `clone`. If negative, 
      then it will be set as the number of cells in `source_cell_type`.

cna_profile_fn : str
    A TSV file listing clonal CNA profiles. 
    It is header-free and its first 6 columns are:
    
    - "chrom" (str): chromosome name of the CNA region.
    - "start" (int): start genomic position of the CNA region, 1-based
      and inclusive.
    - "end" (int): end genomic position of the CNA region, 1-based and
      inclusive.
    - "clone" (str): clone ID.
    - "cn_ale0" (int): copy number of the first allele.
    - "cn_ale1" (int): copy number of the second allele.

refseq_fn : str
    A FASTA file storing reference genome sequence.

out_dir : str
    The output folder.

sam_list_fn : str or None, default None
    A file listing indexed BAM files, each per line.

sample_ids : str or None, default None
    Comma separated sample IDs.
    It should be specified for well-based or bulk data.
    When `barcode_fn` is not specified, the default value will be
    "SampleX", where "X" is the 0-based index of the BAM file(s).
    Note that `sample_ids` and `sample_id_fn` should not be specified
    at the same time.

sample_id_fn : str or None, default None
    A file listing sample IDs, each per line.

overlap_features_how : str, default "raw"
    How to process overlapping features.
    
    - "raw": 
       Leave all input gene annotations unchanged.
    - "quantile": remove highly overlapping genes.
       Remove genes with number of overlapping genes larger than a given
       value (default is the 0.99 quantile among all genes that have 
       overlaps).
    - "union": keep the union range of gene overlaps.
       Replace consecutive overlapping genes with their union genomic 
       range, i.e., aggregate overlapping genes into non-overlapping
       super-genes.

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

kwargs_fit_sf : dict or None, default None
    The additional kwargs passed to function 
    :func:`~.marginal.fit_libsize_wrapper` for fitting size factors.
    The available arguments are:
    
    - dist : {"lognormal", "swr", "normal", "t"}
        Type of distribution.
    If None, set to `{}`.

kwargs_fit_rd : dict or None, default None
    The additional kwargs passed to function 
    :func:`~.marginal.fit_RD_wrapper` for fitting read depth.
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
        
    If None, set to `{}`.

chroms : str, default "human_autosome"
    Comma separated chromosome names.
    Reads in other chromosomes will not be used for sampling and hence
    will not be present in the output BAM file(s).
    If "human_autosome", set to `"1,2,...22"`.

cell_tag : str or None, default "CB"
    Tag for cell barcodes, set to None when using sample IDs.

umi_tag : str or None, default "UB"
    Tag for UMI, set to None when reads only.

umi_len : int, default 10
    Length of output UMI barcode.

ncores : int, default 1
    Number of cores.

seed : int or None, default 123
    Seed for random numbers.
    None means not using a fixed seed.

verbose : bool, default False
    Whether to show detailed logging information.

min_count : int, default 1
    Minimum aggragated count for SNP.

min_maf : float, default 0
    Minimum minor allele fraction for SNP.

strandness : {"forward", "reverse", "unstranded"}
    Strandness of the sequencing protocol.
    
    - "forward": SE sense; PE R1 antisense and R2 sense;
        e.g., 10x 3' data.
    - "reverse": SE antisense; PE R1 sense and R2 antisense;
        e.g., 10x 5' data.
    - "unstranded": no strand information.

min_include : int or float, default 0.5
    Minimum length of included part within specific feature.
    If float between (0, 1), it is the minimum fraction of included length.

multi_mapper_how : {"discard", "duplicate"}
    How to process the multi-feature UMIs (reads).
    
    - "discard": discard the UMI.
    - "duplicate": count the UMI for every mapped gene.

xf_tag : str or None, default "xf"
    The extra alignment flags set by tools like CellRanger or SpaceRanger.
    If set, only reads with tag's value 17 or 25 will count.
    If `None`, turn this tag off.

gene_tag : str or None, default "GN"
    The tag for gene name set by tools like CellRanger or SpaceRanger.
    If `None`, turn this tag off.

min_mapq : int, default 20
    Minimum MAPQ for read filtering.

min_len : int, default 30
    Minimum mapped length for read filtering.

incl_flag : int, default 0
    Required flags: skip reads with all mask bits unset.

excl_flag : int, default -1
    Filter flags: skip reads with any mask bits set.
    Value -1 means setting it to 772 when using UMI, or 1796 otherwise.

no_orphan : bool, default True
    If `False`, do not skip anomalous read pairs.



Input
-----
The inputs to the simulator include:

* Alignment file of seed data (BAM file).
* Cell annotations of seed data (TSV file).
* Feature annotations of seed data (TSV file).
* Phased SNPs of seed data (TSV or VCF file).
* Reference genome sequence of seed data (FASTA file).
* Clone annotations of simulated data (TSV file).
* Clonal CNA profiles of simulated data  (TSV file).


Alignment file of seed data (BAM file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The aligned reads stored in either one single BAM file (from droplet-based 
sequencing platform) or a list of BAM files (from well-based sequencing 
platform).


Cell annotations of seed data (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The cell annotation stored in a header-free TSV file.
Its first two columns are ``cell`` and ``cell_type``, where

cell : str
    Cell barcodes (droplet-based data) or sample ID (well-based data).

cell_type : str
    Cell type.

An example is as follows:

.. code-block::

   AAAGATGGTCCGAAGA-1    immune
   AACCATGTCTCGTATT-1    immune
   AACGTTGTCTCTTGAT-1    epithelial
   AACTCAGAGCCTATGT-1    immune
   AAGACCTAGATGTAAC-1    epithelial
   AAGCCGCTCCTCAATT-1    epithelial


Feature annotations of seed data (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The feature annotation stored in a header-free TSV file.
Its first five columns are ``chrom``, ``start``, ``end``, ``feature``,
and ``strand``, where

chrom : str
    Chromosome name of the feature.

start : int
    Start genomic position of the feature, 1-based and inclusive.

end : int
    End genomic position of the feature, 1-based and inclusive.

feature : str
    Feature name.
    
strand : str
    DNA strand orientation of the feature, "+" (positive) or "-" (negative).

An example is as follows:

.. code-block::

   chr1       29554   31109   MIR1302-2HG     +
   chr1       34554   36081   FAM138A -
   chr1       65419   71585   OR4F5   +
   chr2       38814   46870   FAM110C -
   chr2       197569  202605  AC079779.1      +
   chr3       23757   24501   LINC01986       +


Phased SNPs of seed data (TSV or VCF file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The phased SNPs stored in either a TSV file or a VCF file.

Phased SNPs in TSV format
+++++++++++++++++++++++++
If it is in a TSV file, it should be header-free and its first 6 columns
should be ``chrom``, ``pos``, ``ref``, ``alt``, ``ref_hap``, and 
``alt_hap``, where

chrom : str
    The chromosome name of the SNP.

pos : int
    The genomic position of the SNP, 1-based.

ref : str
    The reference (REF) allele of the SNP, one of ``{'A', 'C', 'G', 'T'}``.

alt : str
    The alternative (ALT) allele of the SNP, one of ``{'A', 'C', 'G', 'T'}``.

ref_hap : int
    The haplotype index of ``ref``, one of ``{0, 1}``.

alt_hap : int
    The haplotype index of ``alt``, one of ``{1, 0}``.
 
An example is as follows:

.. code-block::

   chr1    986336   C       A   0   1
   chr1    1007256  G       A   1   0
   chr1    1163041  C       T   1   0
   chr2    264895   G       C   0   1
   chr2    277003   A       G   0   1
   chr2    3388055  C       T   1   0


Phased SNPs in VCF format
+++++++++++++++++++++++++
If it is in VCF format, the file should contain the ``GT`` in its
``FORMAT`` field (i.e., the 9th column).
The corresponding phased genotype could be delimited by either ``'/'`` or
``'|'``, e.g., "0/1", or "0|1".

.. note::
   * As reference phasing, e.g., with Eagle2, is not perfect, one UMI may 
     cover two SNPs with conflicting haplotype states.
   * Reference phasing tends to have higher rate in longer distance.
     Therefore, further local phasing (e.g., in gene level) and global phasing
     (e.g., in bin level) could be used to reduce error rate, e.g., with the
     3-step phasing used by CHISEL_ in scDNA-seq data and XClone_ in scRNA-seq
     data.
     

Reference genome sequence of seed data (FASTA file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The sequence of reference genome, e.g., the human genome version hg38, 
should be stored in a FASTA file.
Its version should match the one used for generating the alignment (BAM)
file of seed data.


Clone annotations of simulated data (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone annotation stored in a header-free TSV file.
Its first 3 columns should be ``clone``, ``source_cell_type``, and ``n_cell``,
where

clone : str
    The clone ID.

source_cell_type : str
    The source cell type of ``clone``.

n_cell : int
    Number of cells in the ``clone``.
    If negative, then it will be set as the number of cells in 
    ``source_cell_type``.
 
An example is as follows:

.. code-block::

   clone1_normal    immune  -1
   clone2_normal    epithelial  -1
   clone3_cancer    epithelial  -1
   clone4_cancer    epithelial  -1
   clone5_cancer    epithelial  -1

.. note::
   The simulator is designed for diploid genome.
   Generally, it is recommended to use normal cells as ``source_cell_type``
   for simulation of somatic CNAs.


Clonal CNA profiles of simulated data (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The clonal CNA profile stored in a header-free TSV file.
Its first 6 columns should be ``chrom``, ``start``, ``end``,
``clone``, ``cn_ale0``, and ``cn_ale1``, where

chrom : str
    The chromosome name of the CNA region.

start : int
    The start genomic position of the CNA region, 1-based and inclusive.

end : int or "Inf"
    The end genomic position of the CNA region, 1-based and inclusive.
    To specify the end of the whole chromosome, you can use either the actual
    genomic position or simply ``Inf``.

clone : str
    The clone ID.

cn_ale0 : int
    The copy number of the first allele (haplotype).

cn_ale1 : int
    The copy number of the second allele (haplotype).
 
One clone-specific CNA per line.
An example is as follows:

.. code-block::

   chr8 1   Inf clone3_cancer   1   2
   chr6 1   Inf clone4_cancer   0   1
   chr8 1   Inf clone4_cancer   1   2
   chr6 1   Inf clone5_cancer   1   0
   chr8 1   Inf clone5_cancer   1   2
   chr11    1   Inf clone5_cancer   2   0


**Support all three major CNA types**

By specifying different values for ``cn_ale0`` and ``cn_ale1``, you may
specify various CNA types, including copy gain (e.g., setting ``1, 2``), 
copy loss (e.g., setting ``0, 1``), LOH (e.g., setting ``2, 0``).

**Support allele-specific CNA**

This format fully supports allele-specific CNAs.
For instance, to simulate the scenario that two subclones have copy loss in
the same region while on distinct alleles, setting ``cn_ale0, cn_ale1``
to ``0, 1`` and ``1, 0`` in two subclones, respectively, as the example of
copy loss in chr6.

**Support whole genome duplication (WGD)**

It also supports whole genome duplication (WGD), e.g., by setting 
``cn_ale0, cn_ale1`` of all chromosomes to ``2, 2``.
Generally, detecting WGD from scRNA-seq data is challenging, as it is hard
to distinguish WGD from high library size.
One scenario eaiser to detect WGD is that a balanced copy loss occurred 
after WGD, e.g., setting ``cn_ale0, cn_ale1`` of chr3 to ``1, 1``, while
``2, 2`` for all other chromosomes.
In this case, chr3 may have signals of balanced BAF while copy-loss RDR,
which should not happen on normal diploid genome.

**Notes**

* All CNA clones ``clone`` in this file must be in the clone annotation file.
* Only the CNA clones are needed to be listed in this file. Do not list normal
  clones in this file.



Output
------
The final output is available at folder ``{out_dir}/4_rs``.
It contains

* Simulated alignment file (BAM file).
* Simulated cell annotation (TSV file).


Simulated alignment file (BAM file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The simulated reads stored in either one single BAM file (from droplet-based
sequencing platform) or a list of BAM files (from well-based sequencing 
platform).
The BAM file(s) are available at folder ``{out_dir}/4_rs/bam``.


Simulated cell annotation (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The simulated cell annotation stored in a header-free TSV file, located at
``{out_dir}/4_rs/rs.cell_anno.tsv``.
It has two columns ``cell`` and ``clone``, where

cell : str
    The cell barcode (droplet-based data) or sample ID (well-based).

clone : str
    The clone ID.

Note that there is a one-column TSV file storing ``cell`` (cell barcodes or
sample ID) only, located at ``{out_dir}/4_rs/rs.samples.tsv``.



Implementation
--------------
The simulator outputs simulated haplotype-aware alignments for clonal single 
cells based on user-specified CNA profiles, by training on input BAM files.

It mainly includes four modules:

#. ``pp``: preprocessing.
#. ``afc``: allele-specific feature counting.
#. ``cs``: count simulation.
#. ``rs``: read simulation.


The ``pp`` (preprocessing) module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module is implemented in the function ``pp.main.pp_wrapper()``.
The results of this module are stored in the folder ``{out_dir}/1_pp``.

It preprocesses the inputs, including:

* Check and merge overlapping features in the input feature annotation file.
* Check and merge overlapping CNA profiles in the input clonal CNA profile 
  file.


The ``afc`` (allele-specific feature counting) module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module extracts and counts allele-specific UMIs/reads in single cells.

It is implemented in the function ``afc.main.afc_wrapper()``.
The results of this module are stored in the folder ``{out_dir}/2_afc``.

To speedup, features are splitted into batches for multi-processing.
In one feature, the haplotype state of each UMI/read is inferred by
integrating haplotype information from all SNPs covered by the UMI/read.

The output allele-specific *feature x cell* count matrices are at folder 
``{out_dir}/2_afc/counts``.

Additionally, all the count matrices are also saved into one anndata ".h5ad"
file, ``{out_dir}/2_afc/afc.counts.cell_anno.h5ad``, which will be used by 
downstream ``cs`` module.


The ``cs`` (count simulation) module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module simulates new allele-specific *cell x feature* count matrices
based on existing matrices.

It is implemented in the function ``cs.main.cs_wrapper()``.
The results of this module are stored in the folder ``{out_dir}/3_cs``.

This module processes the count matrices of haplotypes "A", "B", "U",
separately, mainly following three steps:

#. Fit feature-specific counts with a specific distribution.
#. Update the fitted feature-specific parameters based on the CNA profile.
#. Generate new feature-specific counts based on the updated parameters.


The ``rs`` (read simulation) module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module simulates new reads for new clonal single cells by sampling reads
from the input BAM file(s) according to the simulated counts.

It is implemented in the function ``rs.main.rs_wrapper()``.
The results of this module are stored in the folder ``{out_dir}/4_rs``.

Specifically, it includes following steps:

#. Sample *cell x feature* CUMIs based on simulated counts.
#. Extract output reads according to the sampled CUMIs.

The output reads of all chromosomes will be merged into new BAM file(s) and
stored in folder ``{out_dir}/4_rs/bam``.



.. _CHISEL: https://www.nature.com/articles/s41587-020-0661-6
.. _XClone: https://www.biorxiv.org/content/10.1101/2023.04.03.535352v2

