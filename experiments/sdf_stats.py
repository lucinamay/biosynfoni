import argparse
import os, sys
from enum import Enum

from tqdm import tqdm
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdchem


def cli():
    parser = argparse.ArgumentParser()

    parser.add_argument("sdf", type=str, help="input sdf file")
    args = parser.parse_args()

    return args


def get_all_properties(suppl: Chem.SDMolSupplier) -> dict[str, np.ndarray]:
    properties = {}

    nones = np.zeros(len(suppl))
    nones.fill(np.nan)

    properties["molecular_weight"] = np.copy(nones)
    properties["n_atoms"] = np.copy(nones)
    properties["carbons"] = np.copy(nones)
    properties["nitrogens"] = np.copy(nones)
    properties["oxygens"] = np.copy(nones)
    properties["halogens"] = np.copy(nones)
    properties["heteroatoms"] = np.copy(nones)

    properties["double_bonds"] = np.copy(nones)
    properties["triple_bonds"] = np.copy(nones)
    properties["aromatic_bonds"] = np.copy(nones)
    properties["rings"] = np.copy(nones)
    properties["ring_atoms"] = np.copy(nones)
    properties["ring_bonds"] = np.copy(nones)
    properties["rotatable_bonds"] = np.copy(nones)

    for i, mol in tqdm(enumerate(suppl), total=len(suppl)):
        if mol is None:
            continue
        properties["molecular_weight"][i] = Descriptors.ExactMolWt(mol)
        properties["n_atoms"][i] = mol.GetNumAtoms()
        # get number of nitrogens:

    return properties


def main():
    args = cli()

    sdf_path = os.path.abspath(args.sdf)
    folder = os.path.dirname(sdf_path)
    name = os.path.basename(sdf_path).split(".")[0]

    iwd = os.getcwd()
    os.chdir(folder)

    suppl = Chem.SDMolSupplier(sdf_path)
    properties = get_all_properties(suppl)
    for key, val in properties.items():
        np.savetxt(f"{name}_{key}.tsv", val, fmt="%s", delimiter="\t")

    os.chdir(iwd)
    exit(0)
    return None


if __name__ == "__main__":
    main()
