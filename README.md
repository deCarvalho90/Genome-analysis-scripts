# Genome-analysis-scripts
Python scripts to help comparative genomics analyses

Scripts in this folder:


a) get_gene_equivalent_corresponding.py - Compares the genomes (CDS) of the same species and searches for single nucleotide variations between corresponding genes.


b) recip_finder_tetraploids - Uses a LASTZ output file between 2 species where one of them is an allotetraploid and looks for the 2 best hits of the allotetraploid (compared to the other species).


c) permutation_test_gene_expression - Performs a permutation test using gene expression data to check if the expression of a certain gene is higher or lower than expected.

d) remove_overlap_hmmscan.py - Assembles domain architecture for proteins, removing any overlapping domains and domains with i-evalue < 0.001. Input: hmmscan domtblout output format.
