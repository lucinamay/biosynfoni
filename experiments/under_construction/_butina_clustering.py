import argparse, os, sys
from sys import argv
from enum import Enum

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
import matplotlib.pyplot as plt

from biosynfoni.fingerprints import array_to_bitvect


def distance_matrix(fps: list[DataStructs.ExplicitBitVect]) -> list[float]:
    # create a similarity matrix
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])
    return dists


# # cluster the molecules using Butina clustering
# def butina_clustering(molecules, cutoff=0.35):
#     # generate fingerprints
#     fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in molecules]

#     # create a similarity matrix
#     dists = []
#     nfps = len(fps)
#     for i in range(1, nfps):
#         sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
#         dists.extend([1 - x for x in sims])

#     # cluster the molecules
#     picker = MaxMinPicker()
#     clusters = list(picker.ClusterData(dists, nfps, cutoff, isDistData=True))

#     return clusters


# get ten representative molecules that are central to cluster


# get random molecules from each cluster
def get_random(molecules, clusters, n=10):
    from random import sample

    # get the random molecules
    random_molecules = []
    for i, cluster in enumerate(clusters):
        if len(cluster) >= n:
            indices = sample(cluster, n)
        else:
            indices = cluster
        random_molecules.extend([molecules[j] for j in indices])

    return random_molecules


def plot_cluster_size(clusters):
    fig = plt.figure(1, figsize=(18, 5))
    plt1 = plt.subplot(111)
    plt.axis([0, 151, 0, len(clusters[0]) + 1])
    plt.xlabel("Cluster index")
    plt.ylabel("Number of molecules")
    plt1.bar(range(1, 151), [len(c) for c in clusters[:150]], lw=0)
    plt.show()


def main():
    fps_path = argv[1]
    sdf_path = argv[2]
    names = argv[3]
    fps = np.loadtxt(fps_path, delimiter=",", dtype=int)
    print(fps.shape)
    fps = array_to_bitvect(fps)
    # dist_matrix = distance_matrix(fps)

    # compounds = [
    #     (Chem.MolFromSmiles(x), x)
    #     for x in np.loadtxt(sdf_path, usecols=0, dtype="str")
    # ]
    compounds = [mol for mol in Chem.SDMolSupplier(sdf_path)]

    fps = [AllChem.GetMACCSKeysFingerprint(m, 2, 2048) for m in compounds]
    dist_matrix = distance_matrix(fps)
    num_fps = len(fps)

    c_names = np.loadtxt(names, usecols=0, dtype="str")
    compounds = zip(compounds, c_names)

    clusters = Butina.ClusterData(
        dist_matrix, num_fps, 0.2, isDistData=True
    )  # distance cutoff = 0.5
    print("number of clusters =", len(clusters))
    num_clust_g5 = len([c for c in clusters if len(c) > 5])
    print("number of clusters with more than 5 compounds =", num_clust_g5)

    print("Ten molecules from first 10 clusters:")
    # Draw molecules
    grid_image = Draw.MolsToGridImage(
        [compounds[clusters[i][0]][0] for i in range(10)],
        legends=[compounds[clusters[i][0]][1] for i in range(10)],
        molsPerRow=5,
    )
    svg_str = grid_image.replace("svg:", "")
    print(svg_str)
    # save
    with open("first_10_clusters.svg", "w") as f:
        f.write(svg_str)


if __name__ == "__main__":
    main()
