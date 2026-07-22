Implementation
==============

afc module
----------

Read Assignment
~~~~~~~~~~~~~~~
* (Optional) assign a UMI to some gene only if certain fraction (e.g., 90%) 
  of the UMI’s post-QC reads mapped to the gene. 
  It could be useful to help address the issue of multi-gene UMIs, 
  considering one UMI that 70% of its reads mapped to one gene, and the rest 
  30% reads mapped on a different gene.
  * add option ``--minRAF`` (minimum read assignment fraction)?
  * Other strategies can be used to process the above example 70%-30% UMI, 
    e.g., assign the UMI to the gene with more supporting reads and discard 
    the non-supporting reads, instead of using a hard cutoff of fraction.


Multi-gene UMIs
~~~~~~~~~~~~~~~
* Use pandas to load whole file in submodule ``afc.mfu``, to speedup. 
  The divide-conquer strategy, which splits the cell or feature files into 
  small batches for parallel processing and then merge, can be time-consuming
  in large datasets, such as the HCC3 dataset.


Haplotype inference
~~~~~~~~~~~~~~~~~~~
While the current four-layer hierarchy could be conceptually simplified by 
directly aggregating SNP-level calls to the UMI level, 
we retain the full layered structure to preserve modularity and extensibility. 
Specifically, the single-read (sread) level constitutes a technically 
meaningful analytical unit that provides an informative intermediate tier for 
haplotype assignment. 
Furthermore, this framework delivers extensibility by supporting stage-specific
threshold tuning, such as setting minimum consensus voting rates, 
at each step of the inference hierarchy.

TODO:

* SNP-level: use UMI-specific consensus allele (base) of the SNP, instead of
  read-specific base at SNP position, for haplotype assignment;
* sread, gread, and UMI level: use state-specific min consensus rate for
  haplotype assignment.


For reference, current (v0.6.0) haplotype inference workflow is:

To resolve the haplotype state of each UMI, allelic assignment proceeds 
through four nested hierarchical aggregation layers: 
individual phased SNPs – single reads – read groups - full UMIs.

1. SNP level: For any read spanning a given SNP coordinate, 
   we first assign a SNP-level haplotype state by comparing the observed base 
   at the SNP position against the reference and alternate alleles linked to 
   each parental haplotype.
   Four mutually exclusive SNP-level states are defined:

   - A: Observed base matches the reference haplotype (haplotype index 0).
   - B: Observed base matches the alternative haplotype (haplotype index 1).
   - O (Other): Observed base matches neither the reference nor alternative 
     allele and thus lacks support for either parental haplotype.
   - U (Unknown): No base is fetched, e.g., the SNP lies within a split read.

2. Single read (sread) level:
   A single read can cover multiple SNPs, each contributing an independent 
   haplotype vote.
   The read's overall haplotype state is determined by the consensus of all 
   SNP votes.
   The five possible states are:

   - A: No SNPs are in state B; and at least one SNP is in state A.
   - B: No SNPs are in state A; and at least one SNP is in state B.
   - D (Dual): The sread carries both A-state and B-state SNPs.
   - O (Other): No SNPs are in state A or B; and at least one SNP is in 
     state O.
   - U (Unknown): Either (1) the read covers no phased SNPs; or 
     (2) No SNPs are in state A, B or O; and at least one SNP is in state U.

3. Read group (gread) level: All reads sharing the same query name are grouped
   into a gread.
   For paired-end sequencing, a gread consists of matching R1 and R2 mate 
   reads; 
   for single-end sequencing, each sread constitutes its own gread.
   The gread haplotype state is the consensus of its constituent reads.

   - A: No sreads are in state B or D; and at least one sread is in state A.
   - B: No sreads are in state A or D; and at least one sread is in state B.
   - D: Either (1) at least one sread is in state D, or (2) the UMI group 
     contains both A-state and B-state sreads.
   - O: No sreads are in state A, B, or D; and at least one sread is in 
     state O.
   - U: All sreads are in state U.

4. UMI level: All greads sharing the same cell barcode and UMI sequence form 
   one UMI-level read group.
   The final UMI haplotype state is assigned via consensus across all its 
   constituent greads.
   The five UMI-level states are:

   - A: No greads are in state B or D; and at least one gread is in state A.
   - B: No greads are in state A or D; and at least one gread is in state B.
   - D: Either (1) at least one gread is in state D, or (2) the UMI group 
     contains both A-state and B-state greads.
   - O: No greads are in state A, B, or D; and at least one gread is in 
     state O.
   - U: All greads are in state U.



cs module
---------
* Write a wrapper function for cs module to make it an independent tool for
  count(RDR)-based CNA simulatioin, to take CNA profile and clone annotation
  files as input.



Tests
=====
* Prepare notebooks & scripts to test each module.
