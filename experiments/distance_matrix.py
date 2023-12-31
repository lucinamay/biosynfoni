import sys, os

import numpy as np
from tqdm import tqdm

from biosynfoni.fingerprints import counted_tanimoto_sim
from biosynfoni.inoutput import outfile_namer

def distance_matrix(data: np.array, metric: str = "euclidean") -> np.array:
    """returns distance matrix of array"""
    if metric == "euclidean":
        return np.array([[np.linalg.norm(i - j) for j in data] for i in tqdm(data)])
    elif metric == "tanimoto":
        return np.array([[1.0 - counted_tanimoto_sim(i, j) for j in data] for i in tqdm(data)])
    elif metric == "manhattan":
        print("warning: not checked!! check formula for manhattan distance")
        return np.array([[np.linalg.norm(i - j, ord=1) for j in data] for i in tqdm(data)])

def main():
    filename = sys.argv[1]
    metric = "tanimoto"
    print(f"loading data from {filename}...")
    data = np.loadtxt(filename, delimiter=",", dtype=int)
    print(f"calculating distance matrix with {metric} metric...")
    distances = distance_matrix(data, metric=metric)
    print(f"saving distance matrix to {outfile_namer(filename, metric)}...")
    np.savetxt(outfile_namer(filename, metric), distances, delimiter=",", fmt="%.8f")
    print(f"succesfully saved distance matrix to {outfile_namer(filename, metric)}")


if __name__ == "__main__":
    main()
