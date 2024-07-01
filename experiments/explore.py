from enum import Enum

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from biosynfoni.subkeys import get_names

# Load the data
chebi = pd.read_csv(
    "~/thesis/input/chebi/chebi_overlap_biosynfoni.csv", dtype=int, header=None
)
# classes = np.loadtxt('~/thesis/input/chebi/chebi_classes.tsv', delimiter=',', dtype='str')
classes = pd.read_csv("~/thesis/input/chebi/chebi_classes.tsv", sep="\t", header=None)
classes = classes[0].str.split(",", expand=True)
classes.loc[classes[1] == "isoprenoid", [0, 1]] = (
    classes.loc[classes[1] == "isoprenoid"][[1, 0]]
).values
classes = classes.astype("category")

chebi = pd.DataFrame(chebi)
chebi.columns = get_names()
substructure_names = chebi.columns

chebi[[f"class_{i}" for i in classes.columns]] = classes
chebi

chebi_fp = chebi[substructure_names]
chebi_fp
chebi


fig, axes = plt.subplots(13, 3, figsize=(12, 54))
for i, col in enumerate(chebi_fp.columns):
    # plot the distribution of col per class in chebi[class_0], normalised
    ax = axes[i // 3, i % 3]
    legend = True if i == 1 else False
    sns.kdeplot(data=chebi, x=col, hue="class_0", fill=True, ax=ax, legend=legend)
    ax.set_title(col)
    ax.set_xlabel("")
    # ax.set_ylabel('')
