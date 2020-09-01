import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import scipy.stats as stats
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42

import pybedtools


###########################################

def get_coordinates(fileName):
    loci_dict ={}
    with open(fileName, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            loci_dict[values[3]] = [values[0], values[1], values[2], values[4]]
    return loci_dict


def read_up_file(up_file):
# Function to read the peak file and return a dict with each chromosome as the key of a list wit the lengths of the peaks in that chromosome.
    up_list = []
    with open(up_file, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            up_list.append(values[3])
    f.close()
    return up_list


def read_bed(bed_file):
# Function to read the peak file and return a dict with each chromosome as the key of a list wit the lengths of the peaks in that chromosome.
    chrom_dict = {}
    with open(bed_file, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            p_len = int(values[2]) - int(values[1])
            if values[0] in chrom_dict:
                chrom_dict[values[0]].append(p_len)
            else:
                chrom_dict[values[0]]= [p_len]
    f.close()
    return chrom_dict


def read_genome_file(genome_file):
# Function to read the peak file and return a dict with each chromosome as the key of a list wit the lengths of the peaks in that chromosome.
    chrom_dict = {}
    with open(genome_file, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            chrom_dict[values[0]]= int(values[1])
    f.close()
    return chrom_dict


###############################################################################
# Genome_file
genome_file = "hg38.genome"

prefix = "shM-shC.L1a"
prefix2 = "shM-shC"

# File with the upregulated loci in bed format
upLoci_file = "log2FC1/up_"+prefix+".bed"
upLoci_BT = pybedtools.BedTool(upLoci_file).sort()


# Files with the annotation all genes
gene_file = "genes_tss.Chr.sorted.bed"
gene_BT = pybedtools.BedTool(gene_file).sort()

# ISGs
isg_file = "../gene_counts/log2FC1/up_"+prefix2+".isgs.tss.bed"
isg_BT = pybedtools.BedTool(isg_file).sort()

############################################################
# Distance for ISGs:
temp2 = upLoci_BT.closest(isg_file, d=True)
temp2.saveas("up_"+prefix+".closest_to_up-ISGs.txt")

distance_isg = np.array([])
for t in temp2:
    dist = int(t.fields[10])
    if dist == -1:
        dist = np.nan
    distance_isg = np.append(distance_isg, dist)

median_isgs = np.nanmedian(distance_isg)

############################################################ Tenorios91Metepec85
### For randomizations with coordinates
loci_coords = get_coordinates(gene_file)
# Make array with loci names
gene_ids = np.array([r for r in loci_coords])

# How many genes in total?
g_len = len(gene_ids)

# How many were ISG?
choose_n = len(upLoci_BT)

# Initializing array for iterations:
n_iter = 1000
flat_median = np.array([])
for i in range(n_iter):
    # random genes
    temp_string = ''
    for n in range(choose_n):
        idx = np.random.randint(0, g_len)
        random_id = gene_ids[idx]
        r = loci_coords[random_id]
        temp_string = temp_string + r[0] +'\t'+ r[1] +'\t'+ r[2] + '\t.\t' + r[3] + '\n'
    # Creating the bedtool object and getting the number of overlapping intervals
    temp_BT = pybedtools.BedTool(temp_string, from_string=True).sort().saveas("temp.bed")
    # Distance to random genes:
    temp3 = upLoci_BT.closest("temp.bed", d=True)
    distance_rand = np.array([])
    for t in temp3:
        dist = int(t.fields[10])
        if dist == -1:
            dist = np.nan
        distance_rand = np.append(distance_rand, dist)
    flat_median = np.append(flat_median, np.nanmedian(distance_rand))


p_val = float(len(flat_median[flat_median<=median_isgs])+1)/(n_iter+1)


############################################################
# Making DF for boxplot:

distance_isg = np.log10(distance_isg+0.0001)
distance_rand = np.log10(distance_rand+0.0001)
flat_median = np.log10(flat_median+0.0001)

gene_type_list = ["up-ISGs" for i in range(len(distance_isg))] + ["example\nrandomization" for i in range(len(distance_rand))] + ["median of distances over\n1,000 randomizations"]
distance_combined = np.append(np.append(distance_isg, distance_rand), np.nan)
# Making dataframe
isg_df = pd.DataFrame(np.column_stack([gene_type_list, distance_combined]), columns=["gene_type", "distance"])
isg_df['distance'] = isg_df['distance'].astype('float')
isg_df['gene_type']= isg_df['gene_type'].astype(basestring)

# making df for swarmplot:
type_list = ["up-ISGs"]+ ["example\nrandomization"] + ["median of distances over\n1,000 randomizations" for i in range(len(flat_median))]
distance_combined = np.append(np.append(np.nan, np.nan), flat_median)
# Making dataframe
rand_df = pd.DataFrame(np.column_stack([type_list, distance_combined]), columns=["gene_type", "distance"])
rand_df['distance'] = rand_df['distance'].astype('float')
rand_df['gene_type']= rand_df['gene_type'].astype(basestring)


flatui = ["#9b59b6", "#3498db", "#95a5a6", "#2ecc71"]
f, ax = plt.subplots()
cru = sns.boxplot(x="gene_type", y="distance", data=isg_df, notch=True,  palette=sns.color_palette(flatui), fliersize=0, )
cru = sns.stripplot(x="gene_type", y="distance", data=rand_df, size=3.5, alpha=0.6,  palette=sns.color_palette(flatui))
cru.text(2, 8, ("%.4f" % p_val),  horizontalalignment='left', size='medium', color='black')
plt.xticks(rotation=90)
plt.xlabel("", fontsize=12)
plt.ylabel("log10(distance)", fontsize=12)
plt.ylim(0,10)
f.subplots_adjust(bottom=0.45)
f.savefig(prefix+"_distance_loci_to_up_ISGs.median_of_rand.pdf")
plt.close()



