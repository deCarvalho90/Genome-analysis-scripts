# Script to find best reciprocal lastz hits using sorghum as query
# This script retrieves the best 2 hits in the case of an allotetraploid species
# Ex: speciesA gene1 is a hit with both speciesB gene1 and gene2
# Input: LASTZ output file
import sys

fh = open(sys.argv[1])

sb_dict = {}
pm_dict = {}

for x in fh:
	if x[0] == '#': continue
	y = x.strip().split('\t')
	sbname = y[1]
	pmname = y[6]
	if not sbname in sb_dict:
		sb_dict[sbname] = {}
	if not pmname in sb_dict[sbname]: 
		sb_dict[sbname][pmname] = int(y[0])
	else:
		if int(y[0]) > sb_dict[sbname][pmname]:
			sb_dict[sbname][pmname] = int(y[0])

	if not pmname in pm_dict: 
		pm_dict[pmname] = {}
	if not sbname in pm_dict[pmname]: 
		pm_dict[pmname][sbname] = int(y[0])
	else:
		if int(y[0]) > pm_dict[pmname][sbname]:
			pm_dict[pmname][sbname] = int(y[0])

for sbgene in list(sb_dict):
	pmgenelist = list(sb_dict[sbgene])
	pmgenelist.sort(key=lambda a:sb_dict[sbgene][a])
	pmbesthit1 = pmgenelist[-1]
	pmbesthit2 = 'No Gene'
	if len(pmgenelist) > 1:
		pmbesthit2 = pmgenelist[-2]
	sbgenelist = list(pm_dict[pmbesthit1])
	sbgenelist.sort(key=lambda a:pm_dict[pmbesthit1][a])
	sbbesthit = sbgenelist[-1]
	if sbbesthit == sbgene:
		print sbgene,pmbesthit1,pmbesthit2 
