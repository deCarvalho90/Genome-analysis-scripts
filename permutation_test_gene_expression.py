# Script to run a permutation test of gene expression values
# Jun/2016 - C4 plants: Sorghum, maize and Setaria
import random
import numpy as np
import statistics
import matplotlib.pyplot as plt
import sys

syndoms = sys.argv[1] # Synteny spreadsheet - syntenic genes with same domain
sorg_exp = sys.argv[2] # Sorghum Gene expression
maiz_exp = sys.argv[3] # Maize gene expression
set_exp = sys.argv[4] # Setaria gene expression
synonym = sys.argv[5] # Synonym - setaria (convert "Si" format to "Seita" format)

fh1 = open(syndoms)
fh2 = open(sorg_exp)
fh3 = open(maiz_exp)
fh4 = open(set_exp)
fh5 = open(synonym)

sorg = {} # Key = sorghum gene; value = gene expression value
maiz = {} # Key = maize gene; value = gene expression value
sita = {} # Key = new annotation setaria gene; value = gene expression value
set_syn = {} # Key = old annotation setaria; value = new annotation setaria
all_exps = [] # List with all expression values

for line in fh5: # Builds set_syn dictionary
	y = line.strip().split("\t")
	new_anno = y[0]
	new_anno = new_anno[:-2]
	old_anno = y[1]
	if old_anno in set_syn: continue
	set_syn[old_anno] = new_anno

for line in fh2: # Builds sorg dictionary
	y = line.strip().split("\t")
	if "tracking" in y[0]: continue
	gene = y[4]
	exp = y[9]
	sorg[gene] = exp
	#all_exps.append(float(exp))

for line in fh3: # Builds maiz dictionary
	y = line.strip().split("\t")
	if "tracking" in y[0]: continue
	gene = y[4]
	if "." in gene:
		gene = gene.replace(".", "_")
	exp = y[9]
	maiz[gene] = exp
        #all_exps.append(float(exp))

for line in fh4: # Builds sita dictionary
	y = line.strip().split("\t")
	if "tracking" in y[0]: continue
        gene = y[4]
	if gene in set_syn:
		gene = set_syn[gene]
		exp = y[9]
	else:
		continue
        sita[gene] = exp
        #all_exps.append(float(exp))

all_exps = []
for line in fh1: # Creates new lines with the average expression value of C4 plants and syntenic genes.
        y = line.strip().split(";")
        if y[0][0] == "#": continue
        line1 = ""
        sorghum = y[1]

        maiz1 = y[2]

        maiz2 = y[3]

        seita = y[4]
        if sorghum in sorg:
                exp_sorg = float(sorg[sorghum])
        else:
                exp_sorg = 0.0

        if maiz1 in maiz:
                exp_maiz1 = float(maiz[maiz1])
        else:
                exp_maiz1 = 0.0

        if maiz2 in maiz:
                exp_maiz2 = float(maiz[maiz2])
        else:
                exp_maiz2 = 0.0

        if seita in sita:
                exp_sita = float(sita[seita])
        else:
                exp_sita = 0.0
        #print "exp_sorg=", exp_sorg, type(exp_sorg), "exp_maiz1=", exp_maiz1, type(exp_maiz1), "exp_maiz2=", exp_maiz2, type(exp_maiz2), "exp_sita=", exp_sita, type(exp_sita)
        avg_exp = (exp_sorg + exp_maiz1 + exp_maiz2 + exp_sita)/3.0 # Considering 3 species
	all_exps.append(float(avg_exp))

pop_mean = statistics.mean(all_exps)
pop_median = statistics.median(all_exps)
pop_std = statistics.pstdev(all_exps)

#print len(all_exps)

samp_mean = []
samp_median = []
samp_std = []
for x in range(0, 1000):
	perm = np.random.choice(all_exps, 385, replace=False)
	s_mean = np.mean(perm)
	s_median = np.median(perm)
	s_std = np.std(perm)
	samp_mean.append(s_mean)
	samp_median.append(s_median)
	samp_std.append(s_std)
	#print s_stdev
	#pop_median = statistics.median(perm)
	#stdev = statistics.pstdev(perm)

samp_mean.sort()
samp_median.sort()
samp_std.sort()

print "Population mean", pop_mean, "\n", "Sample mean 25th position", samp_mean[24]
print "Sample mean 975th position", samp_mean[974]

print "\n", "Population median", pop_median
print "Sample median 25th position", samp_median[24]
print "Sample median 975th position", samp_median[974]

print "\n", "Population standard deviation", pop_std
print "Sample median 25th position", samp_std[24]
print "Sample median 975th position", samp_std[974]

samp_mean_line1 = samp_mean[24] - 1.3
samp_mean_line2 = samp_mean[974] + 0.5
### Mean histogram ###
fig = plt.figure()
ax = fig.add_subplot('111')
ax.hist(samp_mean, bins=20, color='w')
plt.axvline(20.2717605829, color='b', linestyle='dashed', linewidth=2) #  Actual mean 2 domains
plt.text(20.7, 100, 'Actual mean 2 domains', rotation=90, color='b', size=15)
plt.axvline(28.0328805334, color='g', linestyle='dashed', linewidth=2) #  Actual mean same domains
plt.text(28.5, 100, 'Actual mean same domains', rotation=90, color='g', size=15)
plt.axvline(samp_mean[24], color='r', linestyle='dashed', linewidth=2) # p = 0.025 same domains
plt.text(samp_mean_line1, 100, 'p=0.025 same domain', rotation=90, color='r', size=15)
plt.axvline(samp_mean[974], color='r', linestyle='dashed', linewidth=2) # p = 0.975 same domains
plt.text(samp_mean_line2, 100, 'p=0.975 same domain', rotation=90, color='r', size=15)
plt.title("Sample mean - same domains")
plt.show()

samp_median_line1 = samp_median[24] - 0.18
samp_median_line2 = samp_median[974] + 0.07
### Median histogram ###
fig = plt.figure()
ax = fig.add_subplot('111')
ax.hist(samp_median, bins=20, color='w')
plt.axvline(7.39788333333, color='b', linestyle='dashed', linewidth=2) #  Actual median 2 domains
plt.text(7.45, 100, 'Actual median 2 domains', rotation=90, color='b', size=15)
plt.axvline(8.01254333333, color='g', linestyle='dashed', linewidth=2) #  Actual median same domains
plt.text(8.07, 100, 'Actual median same domains', rotation=90, color='g', size=15)
plt.axvline(samp_median[24], color='r', linestyle='dashed', linewidth=2) # p = 0.025 same domains
plt.text(samp_median_line1, 100, 'p=0.025 same domain', rotation=90, color='r', size=15)
plt.axvline(samp_median[974], color='r', linestyle='dashed', linewidth=2) # p = 0.975 same domains
plt.text(samp_median_line2, 100, 'p=0.975 same domain', rotation=90, color='r', size=15)
plt.title("Sample median - same domains")
plt.show()

samp_std_line1 = samp_std[24] - 9.4
samp_std_line2 = samp_std[974] + 3.5
### Standard deviation histogram ###
fig = plt.figure()
ax = fig.add_subplot('111')
ax.hist(samp_std, bins=20, color='w')
plt.axvline(65.7766179879, color='b', linestyle='dashed', linewidth=2) #  Actual standard deviation 2 domains
plt.text(69.2, 120, 'Actual St Dev 2 domains', rotation=90, color='b', size=15)
plt.axvline(124.937459224, color='g', linestyle='dashed', linewidth=2) #  Actual standard deviation same domains
plt.text(128.2, 120, 'Actual St Dev same domains', rotation=90, color='g', size=15)
plt.axvline(samp_std[24], color='r', linestyle='dashed', linewidth=2) # p = 0.025 same domains
plt.text(samp_std_line1, 120, 'p=0.025 same domain', rotation=90, color='r', size=15)
plt.axvline(samp_std[974], color='r', linestyle='dashed', linewidth=2) # p = 0.975 same domains
plt.text(samp_std_line2, 120, 'p=0.975 same domain', rotation=90, color='r', size=15)
plt.title("Sample Standard deviation - same domains")
plt.show()

#print len(all_exps)
#print perm, len(perm), type(perm)
