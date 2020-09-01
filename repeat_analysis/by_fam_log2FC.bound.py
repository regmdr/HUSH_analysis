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





##########################################################################################


def read_input(fileName):
    fam_dict = {}
    with open(fileName, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            fam_dict[values[0]] = values[1]
    return fam_dict



def get_repeat_fam_dict(fileName, boundList):
    rep_family = {}
    with open(fileName, "r") as f:
        next(f)
        for line in f:
            values = line.strip().split("\t")
            values2 = values[0].split("_")
            if len(values2) == 2:
                rfam = values2[0]
            elif len(values2) == 3:
                rfam = values2[0]+"_"+values2[1]
            else:
                rfam = values2[0]+"_"+values2[1]+"_"+values2[2]
            # Will store:
            #	l2FC, p-adj, pval, baseMean
            #  values[2], values[5], values[4], values[1]
            if values[5] == "NA":
            	values[5] = np.nan
            else:
            	values[1] = np.float(values[1])
            	values[2] = np.float(values[2])
            	values[4] = np.float(values[4])
            	values[5] = np.float(values[5])
            # Checking if instance is bound
            if values[0] in boundList:
                b_flag = 1
            else:
                b_flag = 0
            rvals = np.array([values[2], values[5], values[4], values[1], b_flag])
            if rfam in rep_family:
            	rep_family[rfam] = np.vstack((rep_family[rfam], rvals))
            else:
            	rep_family[rfam] = rvals
    return rep_family


def get_bound(fileName):
    rep_dict = {}
    with open(fileName, "r") as f:
        next(f)
        for line in f:
            values = line.strip().split("\"") 
            rep_dict[values[1]] = 0 
    return rep_dict

##########################################################################################

"""
First argument is a file with list of repeat families to be plotted in the first column and the corresponding label in the second column.
Second argument is the output from deseq2, with pvalues and log2FCs.
Third argument are repeats bound by factor, in gtf. Naming of repeat integrants should be consistent between DEA output and this file.
Fourth argument is the suffix for the output file.
"""


famN_file = sys.argv[1]
file_name = sys.argv[2]
bound_file = sys.argv[3]
suffix = sys.argv[4]

significant_padj = 0.05


# reading families to be plotted
famN_dict = read_input(famN_file)

# Getting bound repeats
bound_dict = get_bound(bound_file)


# Reading data
dea_rep = get_repeat_fam_dict(file_name, bound_dict)

# Initializinf lists of values to be used for plotting
flat_labels = []
flat_log2FC = []
flat_type = [] # now type will be bound vs not

# For each repeat to be plotted
for r in famN_dict:
    if r in dea_rep:
        temp_array = dea_rep[r]
        if len(dea_rep[r].shape) > 1:
            # for each locus
            for l in temp_array:
                flat_log2FC.append(l[0])
                flat_labels.append(famN_dict[r])
                if np.float(l[4]) == 1:
                    flat_type.append("bound")
                else:
                    flat_type.append("not_bound")
        else:
            print("family ", r, " with only one locus")
    else:
        print(r, "missing from locus analysis")



# Making dataframe
rep_df = pd.DataFrame(np.column_stack([flat_labels, flat_log2FC, flat_type]), columns=["labels", "log2FC", "type"])
rep_df['log2FC'] = rep_df['log2FC'].astype('float')
rep_df['labels']= rep_df['labels'].astype(basestring)
rep_df['type']= rep_df['type'].astype(basestring)


sorted_df = rep_df.sort_values(by='type', ascending=False)

# Plotting
l_order = ["L1PA17", "L1PA16", "L1PA15-16", "L1PA15", "L1PA14", "L1PA13", "L1PA12", "L1PA11", "L1PA10", "L1PA8A","L1PA8","L1PA7","L1PA6","L1PA5","L1PA4", "L1PA3","L1PA2", "L1HS"]
# l_order = ["HERV1 (LTRd)", "HERV9 (LTR12)", "HERVH (LTR7)", "HERVK (LTR22A)", "MER48 (LTR)", "MER92A (LTR)", "ERVL-MaLR (MLT1)", "MamGyp (LTR2b)", "ERVL (MER54A)", "ERVL (LTR82A)", "ERV1 (MER101)", "ERV1 (LTR15)", "LTR69", "LOR1 (LTR)", "HERV3"]


# For numbers
counts_df = rep_df.groupby('labels')['type'].value_counts().unstack().fillna(0)
# reorder by l_order
counts_df = counts_df.reindex(l_order)
tot_locus = counts_df.sum(axis=1)


percent_lab = [str(int(counts_df.bound[i]))+"/"+str(int(tot_locus[i])) for i in range(counts_df.shape[0])]



# Plotting parameters

y_pos = -6
f, ax = plt.subplots(figsize=(5,3))
cru = sns.boxplot(x="labels", y="log2FC", data=sorted_df, zorder=1, color="lightgrey", fliersize=0, order=l_order)
cru2 = sns.stripplot(x="labels", y="log2FC", hue="type", data=sorted_df, size=3, alpha=0.4, dodge=True,  palette={"bound": "orange", "not_bound": "darkgrey"}, zorder=1, order=l_order, jitter=0.2)
for line in range(0,len(percent_lab)):
     cru.text(line, y_pos, percent_lab[line], horizontalalignment='center', size='small', color='black', rotation=45)
plt.xticks(rotation=90)
plt.xlabel("", fontsize=11)
plt.ylim(-10,10)
f.subplots_adjust(bottom=0.2)
cru2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
f.savefig(suffix+"_bound-loci.stripplot_dodged.pdf")
plt.close()


