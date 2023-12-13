"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: cOnVERTure (cdk_rdk_parser)  |
language: python                    |
author: Lucina-May Nollen           | 
institute: WUR Bioinformatics       |
____________________________________

||||||||||||  ()()()  |||||||||||||||

description:    parses (COCONUT) CDK-style sdf's for its CNP-ID, SMILEs,
                InChi and molecular formula, returning 

"""
# --------------------------------- IMPORTS-------------------------------------
import argparse
from datetime import date
import os
from enum import Enum

from tqdm import tqdm
from rdkit import Chem
from rdkit import RDLogger  # for muting warnings

from biosynfoni.inoutput import readr, entry_parser


# ================================ FUNCTIONS ===================================


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to SDF file containing molecules.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to SDF file containing rdkit-loadable molecular represenations in the properties.",
    )
    parser.add_argument(
        "-x",
        "--exclusive",
        type=str,
        default=None,
        choices=[None, "smiles", "inchi"],
        help="Represenation of intererst. Default chooses inchi and smiles, with preference of inchi",
    )
    return parser.parse_args()


# ================================ converter ====================================


def propertifier(sdf_entries) -> dict[int, dict[str, str]]:  # 22 sec
    """lst of properties in entry as follows:
    every entry is [[coconut_id],[inchi],[smiles],[molecular_formula],[NPLS]]
    #for later:     [molecular weight], [name],]
    """
    entries_properties = {}
    for ind, entry in tqdm(enumerate(sdf_entries)):
        entries_properties[ind] = {}
        i = 0
        # coconut_id, inchi, smiles, molecular_formula = "", "", "", ""
        # NPLS = ""
        for i in range(len(entry)):
            if entry[i].startswith("> <"):
                key = entry[i].split("<")[-1].strip(">")
                entries_properties[ind][key] = entry[i + 1]
    return entries_properties


def name_change(
    entries_properties: dict[int, dict[str, str]], pre: str, post: str
) -> dict[int, dict[str, str]]:
    """changes names of keys in dictionary if present, usable for typos in sdf annotations"""
    for ind, property_dict in entries_properties.items():
        if key in property_dict.keys():
            property_dict[post] = property_dict.pop(pre)
    return entries_properties


# def molify_lists(entries_properties: dict) -> tuple[list]:
#     """converts molecule representations to Chem.Mol objects.
#     using this function keeps all mols in memory so repr_to_annotated_sdf is
#     recommended over this+sdf_writer
#     input:  (dict) entries_properties -- dict of annotations per entry
#                                         should contain 'inchi' and 'SMILES'
#                                         in keys
#     returns:    (list) all_mols -- mol list of all successful mols in either
#                                 inchi or smiles
#                 (list) inchi_mols -- mol list of all successful inchi mols
#                 (list) smiles_mols -- mol list of all successful smiles mols
#     """
#     any_mols = []
#     inchi_mols = []
#     smiles_mols = []

#     for ind, properties in entries_properties.items():
#         inchi_mol = Chem.MolFromInchi(properties["inchi"])
#         smiles_mol = Chem.MolFromSmiles(properties["SMILES"])
#         if inchi_mol is not None:
#             for key, val in properties.items():
#                 inchi_mol.SetProp(key, str(val))
#             inchi_mol.SetProp("rec_inchi", Chem.MolToInchi(inchi_mol))
#             inchi_mol.SetProp("rec_smiles", Chem.MolToSmiles(inchi_mol))
#             inchi_mol.SetProp("rec_inchikey", Chem.MolToInchiKey(inchi_mol))
#             inchi_mol.SetProp("coco_index", str(ind))
#             inchi_mols.append(inchi_mol)
#         if smiles_mol is not None:
#             for key, val in properties.items():
#                 smiles_mol.SetProp(key, str(val))
#             smiles_mol.SetProp("rec_inchi", Chem.MolToInchi(smiles_mol))
#             smiles_mol.SetProp("rec_smiles", Chem.MolToSmiles(smiles_mol))
#             smiles_mol.SetProp("rec_inchikey", Chem.MolToInchiKey(smiles_mol))
#             smiles_mol.SetProp("coco_index", str(ind))
#             smiles_mols.append(smiles_mol)
#         else:  # these do not end up in the mol collection
#             with open("converture.err", "a") as errors:
#                 errors.write(ind, properties)
#         if inchi_mol:
#             all_mols.append(inchi_mol)
#         elif smiles_mol:
#             all_mols.append(smiles_mol)
#     # output (yes or no stats)
#     return all_mols, inchi_mols, smiles_mols


def repr_to_annotated_sdf(
    entries_properties: dict, new_sdf_name: str, exclusive_repr: str = None
):
    wr = Chem.SDWriter(new_sdf_name)
    # start new file
    with open(new_sdf_name.replace(".sdf", ".err"), "w") as errors:
        errors.write("{\n")

    count = 0
    for ind, properties in tqdm(entries_properties.items()):
        mol = None  # init.

        if not exclusive_repr == "inchi":
            smiles_mol = Chem.MolFromSmiles(properties["SMILES"])
        if not exclusive_repr == "SMILES":
            inchi_mol = Chem.MolFromInchi(properties["inchi"])

        # get any properly converted version of the molecule
        if inchi_mol is not None:
            mol = inchi_mol
        elif smiles_mol is not None:
            mol = smiles_mol

        # add properties
        if mol:
            for key, val in properties.items():
                mol.SetProp(key, str(val))

            # see if rdkit reading changes representations
            recalculations = {
                "inchi": Chem.MolToInchi(mol),
                "SMILES": Chem.MolToSmiles(mol),
                "inchikey": Chem.MolToInchiKey(mol),
            }
            for key, val in recalculations.items():
                if recalculations[key] != properties[key]:
                    mol.SetProp(f"rdk_{key}", val)
            # add corresponding entry index from input sdf
            mol.SetProp("input_sdf_index", str(ind))
            wr.write(mol)
            count += 1
        else:
            # write unconverted mol properties to error file
            with open(new_sdf_name.replace(".sdf", ".err"), "a") as errors:
                errors.write(f"{ind}:{properties},\n")
    wr.close()
    with open(new_sdf_name.replace(".sdf", ".err"), "a") as errors:
        errors.write("}\n")
    return count


# ================================= output =====================================


def sdf_writr(mols, outfile):
    """writes sdf of mols"""
    writer = Chem.SDWriter(outfile)
    for mol in mols:
        writer.write(mol)


# ++++++++++++++++++++++++++++++++++ main ++++++++++++++++++++++++++++++++++++++


def main():
    # mute RDKit warnings
    RDLogger.DisableLog("rdApp.*")
    # input
    args = cli()
    print("\n", 10 * "=", "\n", "cOnVERTURE", "\n", 10 * "=", "\n")
    print(args)
    # handling
    lines = readr(args.input)
    sdf_entries = entry_parser(lines)
    properties = propertifier(sdf_entries)
    # change murko_framework to murcko_framework
    properties_spellchecked = name_change(
        properties, "murko_framework", "murcko_framework"
    )
    # all_mols, ids_props, stats = molifier(properties, representation_info=True)
    count = repr_to_annotated_sdf(
        properties_spellchecked, args.output, exclusive_repr=args.exclusive
    )
    print(f"wrote {count} mols to {args.output}")
    print("~~~bye~~~")

    return None


if __name__ == "__main__":
    main()
