# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:53:10 2023

@author: lucina-may
"""

import sys
from sys import argv

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd

from figuremaking import df_scatterplot

# parameters


def read_coco_biosynfoni(fingerprint_file: str) -> np.array:
    #fingerprint_file = "../arch/0816_BioSynFoni_fpsfull_3.csv"
    data = []
    with open(fingerprint_file, 'r') as f:
        for line in f:
            data.append([int(x) for x in line.strip().split(',')])
    data = np.array(data)
    return data


def pcaer(
        data: np.array,
        n_components: int = 50,
        random_state=1):
    "random_state Pass an int for reproducible results across multiple function calls."
    pca = PCA(n_components=n_components, random_state=random_state)
    pca.fit(data)
    pcaed = pca.transform(data)
    var = pca.explained_variance_
    return pcaed, var


def tsner(
        data: np.array,
        n_components: int = 2,
        perplexity: int = 40,
        n_iter: int = 300,
        *args,
        **kwargs
) -> pd.DataFrame:
    # tsne = TSNE(n_components=n_components, *args, **kwargs)
    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        n_iter=n_iter,
        *args,
        **kwargs
    )
    z = tsne.fit_transform(data)
    df = pd.DataFrame()
    df["comp-1"] = z[:, 0]
    df["comp-2"] = z[:, 1]
    return df


def get_first_ncps_est(annotfile: str):
    npcs = []
    # with open('../arch/0914_COCONUT_DB_rdk_npcs.tsv','r') as cl:
    with open(annotfile, 'r') as cl:
        for line in cl:
            npcs.append(line.strip().split('\t')[0].split(',')[0])
    return npcs


def annotate_df(df: pd.DataFrame, col: str, annot: list) -> pd.DataFrame:
    df[col] = annot
    df[col] = df[col].replace('', "None")
    return df


def get_subset(df: pd.DataFrame, n: int = 10000) -> pd.DataFrame:  # remove
    subset_df = df.sample(n=n)
    # print(subset_df.index().tolist)
    return subset_df


def pca_plot(arr, annotfile='../arch/0914_COCONUT_DB_rdk_npcs.tsv'):
    # pca
    print('starting pca...')
    n_components = len(arr[0])
    n_components = 2
    print('running pca function...')
    pcaed, var = pcaer(arr, n_components=n_components, random_state=1)
    print('making dataframe...')
    pca_df = pd.DataFrame(pcaed)
    pca_df.columns = [f"{_a} ({_b})" for _a, _b in zip(pca_df.columns, var)]
    print('annotating dataframe...')
    df = annotate_df(
        pca_df,
        'npclassifier',
        get_first_ncps_est(annotfile)[:len(arr)])
    print("plotting pca...")
    df_scatterplot(
        df,
        df.columns[0],
        df.columns[1],
        figtitle=f"pca of biosynfoni, coloured on npclassifier prediction {var[0]}{var[1]}",
        marginal_x='box',
        marginal_y='box',
        color='npclassifier',
        filename='pca_bsf',
        auto_open=False,
        #color_discrete_map = fm.COLOUR_DICT['pathways']
    )
    return None


def pcaed_tsne(arr, annotfile='../arch/0914_COCONUT_DB_rdk_npcs.tsv'):
    n_components = len(arr[0])
    n_components = 10
    pcaed, _ = pcaer(arr, n_components=n_components, random_state=1)
    print("running tsne...")
    tsne_comp = tsner(pcaed, n_components=2, perplexity=50, n_iter=500)
    df = annotate_df(
        tsne_comp,
        'npclassifier',
        get_first_ncps_est(annotfile)[:len(arr)])
    df_scatterplot(
        df,
        df.columns[0],
        df.columns[1],
        figtitle="tsne of pca'ed biosynfoni, coloured on npclassifier prediction",
        marginal_x='box',
        marginal_y='box',
        color='npclassifier',
        filename='tsne_pcaed_bsf',
        auto_open=False
        #color_discrete_map = fm.COLOUR_DICT['pathways']
    )
    return None


def main():
    print("hello")
    fingerprintfile = argv[1]
    fingerprintfile = '1008_coconut_bsf/1008_0814_COCONUT_DB_rdk_bsf.bsf'
    arr = read_coco_biosynfoni(fingerprintfile)
    
    annotfile = '../arch/0914_COCONUT_DB_rdk_npcs.tsv'
    annotfile = argv[2]

    pca_plot(arr, annotfile=annotfile)
    pcaed_tsne(arr, annotfile=annotfile)

    print("done")
    sys.exit()
    return None


if __name__ == "__main__":
    main()
