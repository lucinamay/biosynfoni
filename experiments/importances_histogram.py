from sys import argv

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from biosynfoni.subkeys import fpVersions, defaultVersion, get_values, get_pathway
from utils.figuremaking import cat_to_colour
from utils.colours import colourDict
from utils import set_style


bsf_name = defaultVersion
substructure_names = get_values("name", version=bsf_name)
pathways = get_pathway(version=bsf_name)
colors = cat_to_colour(pathways, colourDict["pathways"])
set_style()

importances = np.loadtxt(argv[1], delimiter="\t", dtype=float)
means = np.mean(importances, axis=0)
barplot = plt.bar(substructure_names, means)  # , color="#8C8")
# add standard deviations as error bars
stds = np.std(importances, axis=0)
print(stds.shape)
e1 = plt.errorbar(substructure_names, means, yerr=stds, fmt="o", color="#606060")
e2 = plt.errorbar(substructure_names, means, yerr=stds, fmt="none", color="#606060")

plt.xticks(range(len(substructure_names)), substructure_names, rotation=90)
# set bar colours
for i, bar in enumerate(barplot):
    bar.set_color(colors[i])

plt.xticks(rotation=90)
plt.ylabel("feature importance")
plt.xlabel("substructure")
plt.suptitle("Feature importances for random forest model", weight="bold")
plt.title("(averaged over k-fold cross validation with k=5)", weight="light")
plt.tight_layout()
plt.savefig("importances.png")
