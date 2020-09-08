import os
import sys
import re
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42


def read_values(fileName):
    families = []
    with open(fileName, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            values[0] = values[0].replace("-int", "")
            values[0] = re.sub("_$", "", values[0])
            if values[0] in families:
                next
            else:
                families.append(values[0])
    return families


def repeat_clade(inFile):
    # Function to get the species that contain a given repeat family.
    # Use the order to assign the broadest category
    order=np.array(["Homo_sapiens", "Homo", "Homininae", "Hominidae", "Hominoidea", "Catarrhini", "Simiiformes", "Haplorrhini", "Primates", "Rodentia", "Glires", "Scandentia", "Euarchontoglires", "Carnivora", "Cetartiodactyla", "Laurasiatheria", "Boreoeutheria", "Xenarthra", "Afrotheria", "Eutheria", "Metatheria", "Theria", "Monotremata", "Mammalia", "Amniota", "Tetrapoda"])
    # The file was parsed from the dfam hmm file from RepeatMasker
    clade_dict = {}
    with open(inFile, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            values[0] = values[0].replace("_3end", "")
            values[0] = values[0].replace("_5end", "")
            values[0] = values[0].replace("-int", "")
            values[0] = values[0].replace("_orf2", "")
            # I changed the value stored to only the tax name, for visualization purposes.
            values2 = values[1].split(":")
            tax_name = values2[1]
            if tax_name == "Theria_Mammalia":
                tax_name = "Theria"
            # Since it is only 2 instances, I will add a conditional specifically for ALR and BRS
            if values[0] == "ALR":
                values[0] = "ALR/Alpha"
            elif values[0] == "BSR":
                values[0] = "BSR/Beta"
          # choose the most "ancient" species specificity
            if values[0] in clade_dict:
                if tax_name in order and clade_dict[values[0]] in order:
                    temp_order = np.where(order==tax_name)
                    current_order = np.where(order==clade_dict[values[0]])
                    if temp_order > current_order:
                        clade_dict[values[0]] = tax_name
            else:
                clade_dict[values[0]] = tax_name
    return clade_dict



def match_family_clade(clade_dict, fam_list):
    subset_clade = {}
    for f in fam_list:
        if f in clade_dict:
            subset_clade[f] = clade_dict[f]
            next
        else:
            print f
    return subset_clade



######################################################################


# Family / species relationship from Dfam:

classification_file = "dfam.speciesID.txt"
repeats_by_clade = repeat_clade(classification_file)


# Retrieving human-spacific repeat families:

human_repeat_names_file = "repeat_names.txt"
human_list = read_values(human_repeat_names_file)


# Getting the clades corresponding to human families:

human_clades = match_family_clade(repeats_by_clade, human_list)
fam_list = [human_clades[r] for r in human_clades]

# Making dataframe
counts_df = pd.DataFrame(fam_list, columns=["families"])

sum_df = counts_df["families"].value_counts(normalize=False).rename("sum").reset_index()

order=["Homo_sapiens", "Homo", "Homininae", "Hominidae", "Hominoidea", "Catarrhini", "Simiiformes", "Haplorrhini", "Primates", "Rodentia", "Glires", "Scandentia", "Euarchontoglires", "Carnivora", "Cetartiodactyla", "Laurasiatheria", "Boreoeutheria", "Xenarthra", "Afrotheria", "Eutheria", "Metatheria", "Theria_Mammalia", "Monotremata", "Mammalia", "Amniota", "Tetrapoda"]


df = sum_df.set_index('index').loc[order].reset_index()
df.to_csv("tree_weights.txt", header=False, index=False)

fract_df = counts_df["families"].value_counts(normalize=True).rename("sum").reset_index()
fract_df['sum'] = fract_df['sum']*100
# df = fract_df.set_index('index').loc[order].reset_index()
fract_df.to_csv("tree_weights.proportion.txt", header=False, index=False)


##################################################

### For upregulated families:
# fname="shM-shC"
fname="M-P"

up_repeat_names_file = "up_"+fname+".repeat_names.txt"

up_list = read_values(up_repeat_names_file)


# Getting the clades corresponding to human families:

up_clades = match_family_clade(repeats_by_clade, up_list)
up_fam_list = [up_clades[r] for r in up_clades]


up_counts_df = pd.DataFrame(up_fam_list, columns=["families"])
up_fract_df = up_counts_df["families"].value_counts(normalize=True).rename("sum").reset_index()
up_fract_df['sum'] = up_fract_df['sum']*100
up_fract_df.to_csv("tree_weights.up_proportion."+fname+".txt", header=False, index=False)
