import os
from rdkit import Chem
from tqdm import tqdm
import pandas as pd


def get_mol_props_df(sdf_path) -> pd.DataFrame:
    queries = ["coconut_id", "inchi", "SMILES", "rdk_inchi", "rdk_SMILES"]
    props = []

    supplier = Chem.SDMolSupplier(sdf_path)
    for i, mol in tqdm(enumerate(supplier), total=len(supplier)):
        this_mol_props = {q: "" for q in queries}
        if mol is None:
            print("molecule {} could not be read".format(i))
            continue
        else:
            for prop in queries:
                if mol.HasProp(prop):
                    this_mol_props[prop] = mol.GetProp(prop)
        props.append(this_mol_props)

    return pd.DataFrame(props)


def main():
    coco_mols_path = "/Users/lucina-may/thesis/metacyc/coconut-links.tsv"
    meta_mols_path = "/Users/lucina-may/thesis/metacyc/compound-links.dat"

    if not os.path.exists(coco_mols_path):
        coco_mols = get_mol_props_df(coco_mols_path)
        coco_mols.to_csv(coco_mols_path, sep="\t", index=False)

    coco_mols = pd.read_csv(coco_mols_path, sep="\t")

    meta_mols = np.loadtxt(meta_mols_path, dtype=str, delimiter="\t", usecols=(0, 1, 2))
    meta_mols = pd.DataFrame(meta_mols, columns=["unique_id", "inchi", "SMILES"])
