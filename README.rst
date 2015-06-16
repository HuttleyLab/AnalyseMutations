##########################################
Analysing mutations with log-linear models
##########################################

**************************
Sampling the germline data
**************************

Ran ``sample_snp.py`` to produce a summary table containing relevant flanks, ancestral base, chromosome location and reported allele frequencies for non-somatic mutations.

Then ran ``generate_germline_alignments.ipy`` followed by ``generate_germline_counts.ipy`` to produce separate counts tables for each mutation direction.

*******************************
Sampling the COSMIC cancer data
*******************************

Downloaded `CosmicMutantExport.tsv.gz <sftp://sftp-cancer.sanger.ac.uk/files/grch38/cosmic/v72/CosmicMutantExport.tsv.gz>`_ from COSMIC.

Ran the ``extract_cosmic_snps.py`` script, which limited output to SNPs and relevant fields (like primary histology).

Ran the ``dump_chromosomes.py`` script to obtain the human chromosome sequences, then used ``cosmic_flanks.py`` to generate the data format consistent with that produced for the analysis of germline mutations.

Then ran ``generate_cancer_alignments.ipy`` and ``generate_cancer_counts.ipy``.
