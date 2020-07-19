# Script to find gene equivalent pairs based on the first line of fasta sequence
# Use: Comparison between two CDS from the same species to check for single nucleotide variation
# Input: ffn files
# Output: list of gene names and their equivalents from genomes 311 and 314

import sys
from pyfasta import Fasta

g311 = sys.argv[1] # 311 ffn file
g314 = sys.argv[2] # 314 ffn file

f311 = Fasta(g311)
f314 = Fasta(g314)

snvs = [] # list of number of snvs between pairs of corresponding genes
percent = [] # List of % of SNV in pairs of genes

count = 0
for g1 in f311:
	seq1 = f311[g1]
	st1 = seq1[:60]
#"""
	for g2 in f314:
		seq2 = f314[g2]
		st2 = seq2[:60]
		if st1 == st2: # If the first 60 nucleotides are the same, check if the genes are the same
			if len(seq1) == len(seq2):
				count += 1
				print(g1+';'+g2+';'+str(len(seq1))+';'+str(len(seq2)))
#print(count)