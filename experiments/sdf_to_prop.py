import argparse, os
from enum import Enum

import numpy as np
from rdkit import Chem
from tqdm import tqdm
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")


def cli():
    parser = argparse.ArgumentParser(
        description="Extracts a property from an sdf file and writes it to a tsv file."
    )
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
    args = parser.parse_args()

    args.sdf_path = os.path.abspath(args.sdf_path)
    return args


def sdf_prop_getter(supplier_path: str, prop_name: str = "coconut_id") -> list:
    """
    Get a property from an sdf file.
    """
    with Chem.SDMolSupplier(supplier_path) as supl:
        for i, mol in tqdm(enumerate(supl), total=len(supl)):
            if mol is not None:
                if mol.HasProp(prop_name):
                    prop = mol.GetProp(prop_name)
                    yield i, prop
                else:
                    yield i, None
            else:
                yield i, None


def main():
    args = cli()

    # read sdf
    props = []
    for i, prop in sdf_prop_getter(args.sdf_path, prop_name=args.prop_name):
        if prop is None:
            with open(args.o.replace(".tsv", ".err"), "a") as f:
                f.write(f"{i}\n")
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
