import os
import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc
matplotlib.rcParams['pdf.fonttype'] = 42
import statsmodels.stats.multitest as multitest

import xenaPython as xena

##########################################################################################

def get_pheno_data(file):
	sample_dict = {}
	FILE = open(file, "r")
	for line in FILE:
		values = line.strip().split("\t")
		sample_dict[values[0]] = [values[1], values[2],values[3], values[4], values[5], values[6]]
	return sample_dict


def get_cancer_tissue(file):
	sample_list = []
	FILE = open(file, "r")
	for line in FILE:
			sample_list.append(line.strip())
	return sample_list

##########################################################################################


# I will get the data from Xena:
hub = "https://toil.xenahubs.net"
dataset = "TcgaTargetGtex_RSEM_Hugo_norm_count"

# which gene am I using?
gene_list = ["MPHOSPH8", "IL17A"]

# Getting the phenotype data:
pheno_file = "../TcgaTargetGTEX_phenotype.common.txt"
pheno_data = get_pheno_data(pheno_file)

# Immune subtype analysis:

immune_file = "Subtype_Immune_Model_Based.common.txt"
sampleImm_dict = {}
immune_types = [[], [], [], [], [], []]
FILE = open(immune_file, "r")
for line in FILE:
	values = line.strip().split("\t")
	charL = len(values[1])
	iType = values[1][charL-2]
	sampleImm_dict[values[0]] = ["C"+iType]
	immune_types[int(iType)-1].append(values[0])
FILE.close()



flat_immune_type = []
flat_expression_value = np.array([])
flat_cancer_type = []
# flat_cancer_type = []
for c in range(len(immune_types)):
# Getting list of samples for that cancer
	cancer_samp_list = immune_types[c]
	cancer_vals = xena.dataset_probe_values(hub, dataset,cancer_samp_list,gene_list)
	cancer_np = np.array(cancer_vals[1][0], dtype=float)
	# Removing fields with nan in them:
	cancer_np = cancer_np[~np.isnan(cancer_np)]
	# combining data
	immune_type_sym = "C"+str(c+1)
	immune_type = [immune_type_sym for i in range(len(cancer_np))]
	# cancer type:
	for i in cancer_samp_list:
		if i in pheno_data:
			c_subtype = pheno_data[i][0]
			flat_cancer_type.append(c_subtype )
	# Adding other values
	flat_immune_type = flat_immune_type + immune_type
	flat_expression_value = np.append(flat_expression_value, cancer_np)
	# flat_cancer_type = flat_cancer_type + cancer_type


# Normal from TCGA
# Getting the phenotype data:
pheno_file = "../Tcga_normal_phenotype.txt"
pheno_data2 = get_pheno_data(pheno_file)

normal_samps = [i for i in pheno_data2]
normal_vals = xena.dataset_probe_values(hub, dataset,normal_samps,gene_list)
normal_np = np.array(normal_vals[1][0], dtype=float)
# Removing fields with nan in them:
normal_np = normal_np[~np.isnan(normal_np)]
immune_type = ["normal_tissue" for i in range(len(normal_np))]
# Adding values
flat_immune_type = flat_immune_type + immune_type
flat_expression_value = np.append(flat_expression_value, normal_np)
flat_cancer_type = flat_cancer_type + immune_type


# Making dataframe 1
immune_df = pd.DataFrame(np.column_stack([flat_immune_type, flat_expression_value, flat_cancer_type]), columns=["immune_type", "expression_value", "cancer_type"])
immune_df['expression_value'] = immune_df['expression_value'].astype('float')
immune_df['immune_type']= immune_df['immune_type'].astype(basestring)
immune_df['cancer_type']= immune_df['cancer_type'].astype(basestring)


# Making dataframe 2
immune_df2 = pd.DataFrame(np.column_stack([flat_immune_type, flat_expression_value]), columns=["immune_type", "expression_value"])
immune_df2['expression_value'] = immune_df2['expression_value'].astype('float')
immune_df2['immune_type']= immune_df2['immune_type'].astype(basestring)

immune_df2.to_csv("expression_values.by_subtype.txt", sep="\t", index=False)


median_all = np.median(immune_df2['expression_value'])

plot_dims = (11, 6)
f, ax = plt.subplots()# figsize=plot_dims
cru = sns.boxplot(x="immune_type", y="expression_value", data=immune_df2, notch=True, color='lightgray',)
ax.set_ylim((7,15))
# ax.hlines(median_all, -2, 6, color='gray', linewidth=0.5, linestyle='dashed')
plt.xticks(rotation=90)
plt.subplots_adjust(bottom=0.3)
f.savefig("MPP8_by_immune-subtype.boxplot.pdf")
plt.close()


# Using dataframe to get a general overview of the immune subtypes:
# Plotting
c1 = np.array(immune_df2['expression_value'][immune_df['immune_type']=="C1"])
c2 = np.array(immune_df2['expression_value'][immune_df['immune_type']=="C2"])
c3 = np.array(immune_df2['expression_value'][immune_df['immune_type']=="C3"])
c4 = np.array(immune_df2['expression_value'][immune_df['immune_type']=="C4"])
c5 = np.array(immune_df2['expression_value'][immune_df['immune_type']=="C5"])
c6 = np.array(immune_df2['expression_value'][immune_df['immune_type']=="C6"])
normal = np.array(immune_df2['expression_value'][immune_df['immune_type']=="normal_tissue"])

# Doing the comparison to normal:
subtypes_ls = [c1,c2,c3,c4,c5,c6]
subtypes_labs = ["C1", "C2", "C3", "C4", "C5", "C6"]

pval_lst = []
for s in subtypes_ls:
	pval_lst.append(stats.ttest_ind(s,normal)[1])

stats_to_normal = pd.DataFrame(np.column_stack([subtypes_labs, pval_lst, multitest.fdrcorrection(np.array(pval_lst))[1]]), columns=["immune_type", "pval", "qval"])

stats_to_normal.to_csv("stats_to_normal.by_subtype.txt", sep="\t", index=False)

