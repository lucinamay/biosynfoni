from sys import argv

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from biosynfoni.subkeys import fpVersions, defaultVersion, get_values, get_pathway
from experiments.figure_producing.utils.figures import cat_to_colour
from utils.colours import colourDict
from utils import set_style


def main():
    """
    Plot the feature importances from a random forest model as a barplot
    """
    importances = np.loadtxt(argv[1], delimiter="\t", dtype=float)
    bsf_name = defaultVersion
    substructure_names = get_values("name", version=bsf_name)
    pathways = get_pathway(version=bsf_name)
    colors = cat_to_colour(pathways, colourDict["pathways"])
    if len(substructure_names) != importances.shape[1]:
        print("WARNING: substructure names not equal to importances")
        substructure_names = [f"{i}" for i in range(importances.shape[1])]
        colors = ["#888888" for _ in range(importances.shape[1])]
    set_style()

    means = np.mean(importances, axis=0)

    # set plot size
    # default: 6.4, 4.8
    ratio = importances.shape[1] / 39
    plt.figure(figsize=(ratio * 6.4, 4.8))

    # plot barplot
    barplot = plt.bar(substructure_names, means)
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
    plt.savefig(argv[1].replace(".tsv", ".png"), bbox_inches="tight")
    return None


if __name__ == "__main__":
    main()
