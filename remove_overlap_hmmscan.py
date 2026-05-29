# Remove overlapping domains
# Usage:
# > python3 remove_overlap_hmmscan.py domtblout_output_file_name

import sys
import re

hmm = sys.argv[1] # hmmscan domtblout output
fh = open(hmm)

domains = {} # Key = gene name, value [domains i-eval < 0.001]
for x in fh: # Reads the hmmscan output
	if '#' in x: continue 
	y = re.sub('\s+','\t',x) # Parse hmmscan output. Replace empty spaces for tabs
	z = y.strip().split('\t')
	
	sum_dom = z[0] # Summarized domain name
	gene = z[3]
	start = int(z[17])
	stop = int(z[18])
	ieval = float(z[12])
	func = ' '.join(z[22::])
	#domain = [func, start, stop, ieval]
	domain = [sum_dom, start, stop, ieval]
	if not ieval <= 0.001: continue # Filter out domains with evalue > 0.001
	if not gene in domains:
		domains[gene] = []
	domains[gene].append(domain)
	
for g in domains:
	all_doms = domains[g]
	sorted_doms = sorted(all_doms, key=lambda x: x[1])
	nonover = [] # List of nonoverlapping domains
	
	for dom in sorted_doms:
		s_dom = dom[0] # Summarized domain function
		start = dom[1]
		end = dom[2]
		ev = dom[3] # i-evalue
		
		if len(nonover) < 1: # Adds the first domain when nonover is empty
			nonover.append(dom)
			continue
		if start >= nonover[-1][2]: # Keeps the domains that starts after the end of the previous domain
			nonover.append(dom)
			continue
		if start < nonover[-1][2]: # Overlapping domains. Keep the one with lower i-evalue
			if ev < nonover[-1][3]:
				nonover[-1] = dom
			else:
				continue
	if len(nonover) < 1: continue # Skip lines without significant domains
	ds = []
	for n in nonover:
		s_dom = n[0]
		ds.append(s_dom)
	ds = '||'.join(ds)
	print(g + ';' + ds)
