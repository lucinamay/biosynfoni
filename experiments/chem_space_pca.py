import argparse
import os, sys

import numpy as np
import tqdm as tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdchem
from sklearn.decomposition import PCA
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


def cli():
    parser = argparse.ArgumentParser()

    parser.add_argument("sdf_path", type=str, help="path to data sdf file")

    args = parser.parse_args()

    return args


def mol_to_maccs(mol: Chem.Mol) -> np.ndarray:
    return rdMolDescriptors.GetMACCSKeysFingerprint(mol)


def mol_to_morgan(mol: Chem.Mol) -> np.ndarray:
    return np.array(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2))


def mol_to_morgan_count(mol: Chem.Mol) -> np.ndarray:
    return rdMolDescriptors.GetHashedMorganFingerprint(mol, 2)


def mol_to_morgan_feat(mol: Chem.Mol) -> np.ndarray:
    return rdMolDescriptors.GetMorganFingerprint(mol, 2)


def fp_to_pca(fp: np.ndarray, n_components: int = 2) -> np.ndarray:
    pca = PCA(n_components=n_components)
    return pca.fit_transform(fp)


def set_mpl_style(style_path: str) -> mpl.style:
    return mpl.style.use(style_path)


def main():
    args = cli()

    sdf_path = os.path.abspath(args.sdf_path)
    folder = os.path.dirname(sdf_path)
    name = os.path.basename(sdf_path).split(".")[0]

    iwd = os.getcwd()
    os.chdir(folder)

    suppl = Chem.SDMolSupplier(sdf_path)
    fps = []
    for i, mol in tqdm.tqdm(enumerate(suppl), total=len(suppl)):
        if mol is None:
            continue
        fps.append(mol_to_morgan(mol))

    fps = np.array(fps)
    fps = fps.reshape(fps.shape[0], -1)

    pca = PCA(n_components=2)
    fps = pca.fit_transform(fps)

    plt.scatter(fps[:, 0], fps[:, 1], s=1, c="black")
    plt.savefig(f"{name}_pca.png")
    plt.close()

    os.chdir(iwd)
