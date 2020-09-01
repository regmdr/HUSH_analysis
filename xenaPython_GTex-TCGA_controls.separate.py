
import os
import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
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
		values = line.strip().split("\t")
		values2 = values[1].split(" ; ")
		if len(values2) > 1:
			sample_list.append([values[0], values2[0], values2[1]])	
		else:
			sample_list.append([values[0], values[1]])
	return sample_list


##########################################################################################


# I will get the data from Xena:
hub = "https://toil.xenahubs.net"
dataset = "TcgaTargetGtex_RSEM_Hugo_norm_count"

# which gene am I using?
gene_name = "MPHOSPH8"

# Getting the phenotype data:
pheno_file = "TcgaTargetGTEX_phenotype.txt"
pheno_data = get_pheno_data(pheno_file)

normal_tissue_dict = {}
for samp in pheno_data:
	info = pheno_data[samp]
	if info[3] == "Normal Tissue":
		if info[1] in normal_tissue_dict:
			normal_tissue_dict[info[1]] += 1
		else:
			normal_tissue_dict[info[1]] = 1

# Cancer types to use:
cancer_samps = get_cancer_tissue("TCGA_GTEx.two_controls.sample_relations.txt")

# Initializing lists
flat_labels = []
flat_values = np.array([])
flat_type = []
cancer_median_list = []
matched_median_list = []
normal_median_list = []
ctrl1_mannU_values = []
ctrl2_mannU_values = []
# Going through all cancers
for c in range(len(cancer_samps)):
	cancer_samp_list = []
	healthy_samp_list = []
	matched_samp_list = []
	for samp in pheno_data:
		info = pheno_data[samp]
		if info[1] == cancer_samps[c][0]:
			if info[3] != "Normal Tissue" and info[3] != "Solid Tissue Normal":
				cancer_samp_list.append(samp)
			elif info[3] == "Solid Tissue Normal":
				matched_samp_list.append(samp)
		elif info[1] == cancer_samps[c][1]:
			if info[3] == "Normal Tissue": # I think this is redundant, but it's better to check
				healthy_samp_list.append(samp)
		elif len(cancer_samps[c]) > 2:
			if info[1] == cancer_samps[c][2]:
				if info[3] == "Normal Tissue": # I think this is redundant, but it's better to check
					healthy_samp_list.append(samp)
	genes =["MPHOSPH8", "MORC2"]
	# Getting values for the different sample sublists
	# Cancer
	cancer_vals = xena.dataset_probe_values(hub, dataset,cancer_samp_list,genes)
	cancer_np = np.array(cancer_vals[1][0], dtype=float)
	cancer_np = cancer_np[np.logical_not(np.isnan(cancer_np))]
	cancer_no = len(cancer_np) # cancer_vals[1][0]
	# Matched Normal
	matched_vals = xena.dataset_probe_values(hub, dataset,matched_samp_list,genes)
	matched_np = np.array(matched_vals[1][0], dtype=float)
	matched_np = matched_np[np.logical_not(np.isnan(matched_np))]
	matched_no = len(matched_np) # healthy_vals[1][0]
	# GTEx normal
	healthy_vals = xena.dataset_probe_values(hub, dataset,healthy_samp_list,genes)
	healthy_np = np.array(healthy_vals[1][0], dtype=float)
	healthy_np = healthy_np[np.logical_not(np.isnan(healthy_np))]
	healthy_no = len(healthy_np) # healthy_vals[1][0]
	###
	print cancer_samps[c][0], cancer_no, healthy_no, matched_no
	# Doing some statistics:
	ctrl1_mannU_stat = stats.mannwhitneyu(cancer_np, matched_np, alternative="two-sided")
	ctrl1_mannU_values.append([ctrl1_mannU_stat.pvalue, ctrl1_mannU_stat.statistic])
	ctrl2_mann_stat = stats.mannwhitneyu(cancer_np, healthy_np, alternative="two-sided")
	ctrl2_mannU_values.append([ctrl2_mann_stat.pvalue, ctrl2_mann_stat.statistic])
	# Getting lists for dataframe
	temp_labels = [cancer_samps[c][0] for i in range(cancer_no+matched_no+healthy_no)]
	temp_type = ["cancer" for i in range(cancer_no)]+ ["matched" for i in range(matched_no)] +["normal" for i in range(healthy_no)]
	temp_vals = np.append(cancer_np, matched_np)
	temp_vals = np.append(temp_vals,healthy_np)
	flat_labels = flat_labels + temp_labels
	flat_type = flat_type + temp_type
	flat_values = np.append(flat_values,temp_vals)
	# getting median values:
	test = healthy_vals[1][0]
	cancer_median_list.append(np.nanmedian(cancer_np))
	matched_median_list.append(np.nanmedian(matched_np))
	normal_median_list.append(np.nanmedian(healthy_np))

# Cells - Transformed Fibroblasts

# Making dataframe
gene_df = pd.DataFrame(np.column_stack([flat_labels, flat_values, flat_type]), columns=["labels", "fpkm_values", "type"])
gene_df['fpkm_values'] = gene_df['fpkm_values'].astype('float')
gene_df['labels']= gene_df['labels'].astype(basestring)
gene_df['type']= gene_df['type'].astype(basestring)

# Plotting
m_width = 0.12
s_val = 0.2 # shift 

f, ax = plt.subplots(figsize=(10,6))
cru = sns.stripplot(x="labels", y="fpkm_values", hue="type", data=gene_df, size=3.5, alpha=0.5, dodge=True,  palette={"cancer": "purple", "normal": "gray", "matched" : "magenta"}, zorder=1)
ticks = cru.get_xticks()
for n in range(len(ticks)):
	tick = ticks[n]-(s_val+0.5*s_val)
	cru.plot([tick-m_width, tick+m_width], [cancer_median_list[n], cancer_median_list[n]], lw=1, color='k')
	tick = ticks[n]
	cru.plot([tick-m_width, tick+m_width], [matched_median_list[n], matched_median_list[n]], lw=1, color='k')
	tick = ticks[n]+(s_val+0.5*s_val)
	cru.plot([tick-m_width, tick+m_width], [normal_median_list[n], normal_median_list[n]], lw=1, color='k')
ax.get_legend().set_visible(False)
ax.set_ylim((7,15))
plt.xticks(rotation=90)
plt.xlabel("", fontsize=12)
f.subplots_adjust(bottom=0.5)
f.savefig(gene_name+".stripplot_TCGA-GTEx.two_controls.separated.May20.pdf")
plt.close()




##########################################################################################
# Writting test statistic results:
# ctrl1 = matched
# ctrl2 = gtex

pvals_mannU = np.array([i[0] for i in ctrl1_mannU_values+ctrl2_mannU_values])
padjs_mannU = multitest.fdrcorrection(pvals_mannU)

out_file_name = gene_name+".MannU_TCGA-GTEx.two_controls.separated.May20.txt"
oFILE = open(out_file_name, "w")
n=0
i=0
for entry in ctrl1_mannU_values:
	oFILE.write( "matched\t%s\t%s\t%s\t%s\n" % (cancer_samps[i][0], entry[0], entry[1], str(padjs_mannU[1][n])))
	i += 1
	n += 1

i=0
for entry in ctrl2_mannU_values:
	oFILE.write( "GTex\t%s\t%s\t%s\t%s\n" % (cancer_samps[i][0], entry[0], entry[1], str(padjs_mannU[1][n])))
	i +=1
	n += 1
oFILE.close()