import sys, os
import argparse

import numpy as np
from scipy.spatial.distance import cdist
from tqdm import tqdm

# from ..biosynfoni.fingerprints import counted_tanimoto_sim
# from ..biosynfoni.inoutput import outfile_namer


# def euclidean_distance_matrix(x: np.array, y: np.array) -> float:
#     """custom metric for distance matrix"""
#     return np.sum((X[:, None, :] - Y[None, :, :]) ** 2, axis=-1) ** 0.5
def cli():
    """Command line interface for distance matrix calculation"""
    parser = argparse.ArgumentParser(
        description="Calculate distance matrix from fingerprint file"
    )
    parser.add_argument("fingerprint", help="Fingerprint file")
    parser.add_argument(
        "-m",
        "--metric",
        required=False,
        type=str,
        help="distance metric",
        default="hamming",
    )
    return parser.parse_args()


def counted_tanimoto_sim(fp1: np.array, fp2: np.array) -> float:
    """
    Tanimoto similarity for two fingerprints

        Args:
            fp1: np.array
            fp2: np.array

        Returns:
            float: Tanimoto similarity
    """
    nom = sum(np.minimum(fp1, fp2))  # overlap
    denom = float(sum(np.maximum(fp1, fp2)))  # all bits that are 'on'
    if denom == 0:
        return -1  # fingerprint does not describe either molecule
    else:
        return nom / denom


def distance_matrix(data: np.array, metric: str = "euclidean") -> np.array:
    """returns distance matrix of array"""
    if metric in ["euclidean", "cosine", "manhattan", "hamming"]:
        array = cdist(data, data, metric=metric)
        # print(array.shape)
        # print(array[0])
        # return array
        return cdist(data, data, metric=metric)
    elif metric == "tanimoto":
        return np.array(
            [[1.0 - counted_tanimoto_sim(i, j) for j in data] for i in tqdm(data)]
        )
    # elif metric == "manhattan":
    #     print("warning: not checked!! check formula for manhattan distance")
    #     return np.array([[np.linalg.norm(i - j, ord=1) for j in data] for i in tqdm(data)])


def main():
    args = cli()
    filename = args.fingerprint
    metric = args.metric
    print(f"loading data from {filename}...")
    data = np.loadtxt(filename, delimiter=",", dtype=int)
    print(f"calculating distance matrix with {metric} metric...")
    distances = distance_matrix(data, metric=metric)
    # outfile = outfile_namer(filename, metric)
    outfile = f"distances_{metric}.csv"
    print(f"saving distance matrix to {outfile}...")
    np.savetxt(outfile, distances, delimiter=",", fmt="%.8f")
    print(f"succesfully saved distance matrix to {outfile}")


if __name__ == "__main__":
    main()
