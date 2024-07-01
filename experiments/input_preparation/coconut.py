"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: (c)O(n)VERTure               |
language: python                    |
author: Lucina-May Nollen           | 
institute: WUR Bioinformatics       |
____________________________________

||||||||||||  ()()()  |||||||||||||||

description:    parses SDF files that have more errors in their Mol objects 
                compared to their annotation's InChI and SMILES representations.
                This script uses rdkit to convert the representations to mol objects
                and writes them to a new SDF file. InChI is preferred over SMILES.

"""

# --------------------------------- IMPORTS-------------------------------------
import logging, os, re

from tqdm import tqdm
from rdkit import Chem
from rdkit import RDLogger  # for muting warnings

from biosynfoni.inoutput import readr, entry_parser


# ================================ FUNCTIONS ===================================
# ================================ converter ====================================


def get_mol_props(sdf_entries) -> dict[int, dict[str, str]]:  # 22 sec
    """converts sdf entries to a dictionary of properties per entry
    input:  (list) sdf_entries -- list of sdf entries as parsed by biosynfoni.inoutput.entry_parser
    returns:    (dict) entries_properties -- dict
                    [key] (int) entry index
                    [val] (dict) dictionary of properties
    """
    entries_properties = {}
    for ind, entry in tqdm(enumerate(sdf_entries)):
        entries_properties[ind] = {}

        # coconut_id, inchi, smiles, molecular_formula = "", "", "", ""
        for i, line in enumerate(entry):
            prop_name = re.findall(r"^>  <(.*?)>", line)
            if prop_name:
                prop_name = prop_name[0]
                prop_val = entry[i + 1]
                entries_properties[ind][prop_name] = prop_val

            # if entry[i].startswith(">  <"):
            #     key = entry[i].split("<")[-1].strip(">")
            #     entries_properties[ind][key] = entry[i + 1]

    return entries_properties


def rename_keys(
    entries_properties: dict[int, dict[str, str]], pre: str, post: str
) -> dict[int, dict[str, str]]:
    """changes names of keys in entry-property dictionary if present, usable for typos in sdf annotations"""
    for ind, property_dict in entries_properties.items():
        if pre in property_dict.keys():
            property_dict[post] = property_dict.pop(pre)
    return entries_properties


def entry_props_to_sdf(
    entries_properties: dict,
    new_sdf_name: str,
    exclusive_repr: str = None,
    rem_chir: bool = False,
):
    """converts molecule representations to Chem.Mol objects and writes them to a new sdf file"""
    wr = Chem.SDWriter(new_sdf_name)
    # start new file
    with open(new_sdf_name.replace(".sdf", ".err"), "w") as errors:
        errors.write("{\n")

    count = 0
    for ind, properties in tqdm(entries_properties.items()):
        mol = None  # init.
        inchi_mol, smiles_mol = None, None

        if exclusive_repr != "inchi" and "SMILES" in properties.keys():
            smiles_mol = Chem.MolFromSmiles(properties["SMILES"])
        if exclusive_repr != "SMILES" and "inchi" in properties.keys():
            inchi_mol = Chem.MolFromInchi(properties["inchi"])

        # get any properly converted version of the molecule
        if inchi_mol:
            mol = inchi_mol
        elif smiles_mol:
            mol = smiles_mol

        if mol:
            if rem_chir:
                mol = Chem.RemoveStereochemistry(mol)
                # clean mol
                Chem.SanitizeMol(mol)
            # add properties
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
                    if rem_chir:
                        mol.SetProp(f"nonchir_{key}", val)
                    else:
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


# ++++++++++++++++++++++++++++++++++ main ++++++++++++++++++++++++++++++++++++++


def main():
    # mute RDKit warnings
    RDLogger.DisableLog("rdApp.*")
    # input
    input_sdf = "raw_data/COCONUT_DB.sdf"
    output_sdf = "coconut.sdf"
    assert input_sdf != output_sdf
    assert os.path.exists(input_sdf)

    logging.info("\n", 10 * "=", "\n", "cOnVERTURE", "\n", 10 * "=", "\n")
    logging.info(f"extracting mols from {input_sdf}")
    # handling
    lines = readr(input_sdf)
    sdf_entries = entry_parser(lines)
    props = get_mol_props(sdf_entries)
    # change murko_framework to murcko_framework
    props = rename_keys(props, "murko_framework", "murcko_framework")

    count_number_written = entry_props_to_sdf(
        props,
        output_sdf,
        exclusive_repr=None,
        rem_chir=False,
    )
    logging.info(f"wrote {count_number_written} mols to {output_sdf}")
    logging.info("~~~bye~~~")

    return None


if __name__ == "__main__":
    main()
