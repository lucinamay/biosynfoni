import sys, os

import numpy as np

from biosynfoni.fingerprints import counted_tanimoto_sim

def distance_matrix(data: np.array, metric: str = "euclidean") -> np.array:
    """returns distance matrix of array"""
    if metric == "euclidean":
        return np.array([[np.linalg.norm(i - j) for j in arr] for i in arr])
    elif metric == "tanimoto":
        return np.array([[1.0 - counted_tanimoto_sim(i, j) for j in arr] for i in arr])
    elif metric == "manhattan":
        print("warning: not checked!! check formula for manhattan distance")
        return np.array([[np.linalg.norm(i - j, ord=1) for j in arr] for i in arr]))

def main():
    filename = sys.argv[1]
    data = np.loadtxt(filename, delimiter=",", dtype=int)
    distances = distance_matrix(data, metric="tanimoto")
    np.savetxt("distances.csv", distances, delimiter=",", fmt="%.8f")


if __name__ == "main":
    main()
