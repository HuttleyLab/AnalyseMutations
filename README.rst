##########################################
Analysing mutations with log-linear models
##########################################

**************************
Sampling the germline data
**************************

#. ``sample_ensembl_germline.py``, produces a summary dump of germline SNPs whose flanks match the reference genome. This produces a massive pickle format file (``data/ensembl_snps_79/GERMLINE_flanks_match_ref.txt.gz``).
#. ``sample_snp.py`` using the ``generate_sampled_snps.ipy`` script. ``sample_snp.py`` reads pickled records of the correct consequence type (e.g. missense) to produce a summary table containing alleles, ancestral base, chromosome location and reported allele frequencies for non-somatic mutations. The results of this were dumped into a ``raw`` directory.
#. ``generate_germline_raw_w_flanks.ipy`` calls ``germline_flanks.py``, combines the above results and slices from the entire chromosome sequence to produce the completed SNP dump records which are written into ``raw_with_flanks``.
#. ``generate_germline_alignments.ipy`` followed by ``generate_germline_counts.ipy``, produces separate counts tables for each mutation direction.
#. ``generate_germline_all_counts.ipy``, produces combined counts tables (all mutation directions) for strand asymmetric and strand symmetric conditions.
#. ``generate_germline_long_flanks.ipy``, reads the germline alignments for the CtoT and AtoG mutations, producing counts from larger flanks. Data written to ``../data/ensembl_snps_79/counts/long_flanks``.

*******************************
Sampling the COSMIC cancer data
*******************************

#. Downloaded `CosmicMutantExport.tsv.gz <sftp://sftp-cancer.sanger.ac.uk/files/grch38/cosmic/v72/CosmicMutantExport.tsv.gz>`_ from COSMIC.
#. ``extract_cosmic_snps.py`` script, which limited output to SNPs and relevant fields (like primary histology).
#. ``dump_chromosomes.py`` script to obtain the human chromosome sequences, then used ``cosmic_flanks.py`` to generate the data format consistent with that produced for the analysis of germline mutations.
#. ``generate_cancer_alignments.ipy`` and ``generate_cancer_counts.ipy``.

*********************************
Analysis of neighbourhood effects
*********************************

#. ``do_nbr.ipy``: Statistical analysis of a single group for contributions of neighbouring bases to point mutations. Results from this written to ``../results/ensembl_snps_79/<chrom class>/<seq class>/directions/<direction>``. The malignant melanoma results are written to ``../results/ensembl_snps_79/malignant_melanoma``.
#. ``do_nbr_long_flanks.ipy``: Statistical analysis of a single group for contributions of neighbouring bases to point mutations. Results written to ``../results/ensembl_snps_79/long_flanks/A/intergenic/directions/<direction>``
#. ``do_nbr_compare.ipy``: Statistical analysis comparing neighbourhood effects between groups (tissue types, chromosome classes, sequence classes, sequence strands). Results are written to ``<group_vs_group>`` directories within ``../results/ensembl_snps_79/``, depending on the level of comparison, e.g. ``../results/ensembl_snps_79/A_vs_X``.

****************************
Analysis of mutation spectra
****************************

Specified in the ``do_spectra_compare.ipy`` script. Output files (``spectra_analysis.json``, ``spectra_analysis.log`` and ``spectra_analysis.txt``) are written to the same directories as the corresponding neighbour analysis.

*********************
Generating grid plots
*********************

These encoded in in ``do_grid_draw.ipy``. Includes both malignant melanoma, germline intergenic plus the mutation spectra plots.


***************************************
Producing manuscript tables and figures
***************************************

This is all encoded in a jupyter notebook file (``manuscript_figs_tables.ipynb``). The resulting tables for the manin manuscript are written to ``../results/ensembl_snps_79/manuscript_tables.tex``. Tables and figures for supplementary material are written to ``../results/ensembl_snps_79/supp_materials_tables.tex``.

Figures are copied to the ``figures`` directory.
