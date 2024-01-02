import argparse
import os

import numpy as np
from rdkit import Chem


def cli():
    parser = argparse.ArgumentParser(description="Pick mols from a file")

    parser.add_argument("input_file", type=str, help="Input file of mols (sdf)")
    parser.add_argument(
        "-n", "--number", type=int, default=10, help="Number of mols to pick"
    )
    parser.add_argument("-s", "--seed", type=int, default=None, help="Random seed")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=os.path.join(os.getcwd(), "testmols/"),
        help="Output file",
    )

    args = parser.parse_args()
    return args


def pick_mols(supplier: Chem.SDMolSupplier, number: int, seed=None):
    if seed is not None:
        np.random.seed(seed)
    picked_indexes = np.random.choice(len(supplier), number, replace=False)
    picked_mols = [supplier[i] for i in picked_indexes]
    return picked_mols


def main():
    args = cli()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    supplier = Chem.SDMolSupplier(args.input_file)
    picked_mols = pick_mols(supplier, args.number, args.seed)
    picked_smiles = [Chem.MolToSmiles(mol) for mol in picked_mols]
    picked_inchi = [Chem.MolToInchi(mol) for mol in picked_mols]

    with open(os.path.join(args.output, "picked_smiles.smi"), "w") as f:
        f.write("\n".join(picked_smiles))
    with open(os.path.join(args.output, "picked_inchi.inchi"), "w") as f:
        f.write("\n".join(picked_inchi))
    with Chem.SDWriter(os.path.join(args.output, "picked_mols.sdf")) as writer:
        for mol in picked_mols:
            writer.write(mol)

    exit(0)
    return None


if __name__ == "__main__":
    main()
