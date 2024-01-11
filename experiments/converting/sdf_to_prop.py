import argparse, os
from enum import Enum

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMolDescriptors as rdmd
from tqdm import tqdm
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sdf_path",
        type=str,
        help="path to sdf file",
    )
    parser.add_argument(
        "prop_name",
        type=str,
        help="name of property to extract",
    )
    parser.add_argument(
        "-o",
        type=str,
        default="sdf_prop.tsv",
        help="path to write output to",
    )
    parser.add_argument(
        "-i",
        "--index",
        action="store_true",
        help="if flag is set, write index of molecule to file",
    )
    parser.add_argument(
        "-c",
        "--calculate",
        type=str,
        help="give function of (mol) to calculate property",
    )
    args = parser.parse_args()
    return args


def handle_none(prop_name, mol):
    try:
        if prop_name.lower() == "smiles":
            return Chem.MolToSmiles(mol)
        elif prop_name.lower() == "inchi":
            return Chem.MolToInchi(mol)
        elif prop_name.lower() == "inchikey":
            return Chem.MolToInchiKey(mol)
        elif prop_name.lower() == "formula":
            return rdmd.CalcMolFormula(mol)
        elif prop_name.lower() == "molwt":
            return Descriptors.MolWt(mol)
        elif prop_name.lower() == "numatoms":
            return mol.GetNumAtoms()
        elif prop_name.lower() == "numheavyatoms":
            return mol.GetNumHeavyAtoms()
        elif prop_name.lower() == "numrings":
            return rdmd.CalcNumRings(mol)
        elif prop_name.lower() == "numaromaticrings":
            return rdmd.CalcNumAromaticRings(mol)
        elif prop_name.lower() == "numaliphaticrings":
            return rdmd.CalcNumAliphaticRings(mol)
        elif prop_name.lower() == "numheterocycles":
            return rdmd.CalcNumHeterocycles(mol)
        elif prop_name.lower() == "numaromaticheterocycles":
            return rdmd.CalcNumAromaticHeterocycles(mol)
        elif prop_name.lower() == "numaliphaticheterocycles":
            return rdmd.CalcNumAliphaticHeterocycles(mol)
        elif prop_name.lower() == "numsaturatedheterocycles":
            return rdmd.CalcNumSaturatedHeterocycles(mol)
        elif prop_name.lower() == "numaromaticcarbocycles":
            return rdmd.CalcNumAromaticCarbocycles(mol)
        elif prop_name.lower() == "numaliphaticcarbocycles":
            return rdmd.CalcNumAliphaticCarbocycles(mol)
        elif prop_name.lower() == "numsaturatedcarbocycles":
            return rdmd.CalcNumSaturatedCarbocycles(mol)
        elif prop_name.lower() == "numhba":
            return rdmd.CalcNumHBA(mol)
        else:
            return None
    except:
        return None


def sdf_prop_getter(supplier_path: str, prop_name: str = "coconut_id") -> list:
    with Chem.SDMolSupplier(supplier_path) as supl:
        for i, mol in tqdm(enumerate(supl), total=len(supl)):
            if mol is not None:
                if mol.HasProp(prop_name):
                    prop = mol.GetProp(prop_name)
                    yield i, prop
                else:
                    prop = handle_none(prop_name, mol)
                    yield i, prop
            else:
                yield i, None


def sdf_prop_calc(supplier_path: str, function: str = "Chem.MolToSmiles") -> list:
    with Chem.SDMolSupplier(supplier_path) as supl:
        for i, mol in tqdm(enumerate(supl), total=len(supl)):
            if mol is not None:
                prop = function(mol)
                yield i, prop
            else:
                yield i, None


def main():
    args = cli()

    # get absolute path
    sdf_path = os.path.abspath(args.sdf_path)

    # read sdf
    props = []

    # if args.calculate:
    #     print("calculating property")
    #     sdf_prop_getter = sdf_prop_calc
    #     function = args.calculate
    #     sdf_prop_getter = lambda x: sdf_prop_calc(x, function=function)

    for i, prop in sdf_prop_getter(sdf_path, prop_name=args.prop_name):
        if prop is None:
            with open(args.o.replace(".tsv", ".err"), "a") as f:
                f.write(f"{i}\n")
            continue
        if args.index:
            props.append([str(i), prop])
        else:
            props.append(prop)
    # its faster to append to list and then convert to array than to append to array
    props = np.array(props)
    np.savetxt(args.o, props, delimiter="\t", fmt="%s")

    return None


if __name__ == "__main__":
    main()
