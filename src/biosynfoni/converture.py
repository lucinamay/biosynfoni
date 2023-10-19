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
                [0] a sdf of mol objects
                [1] a tsv with related info for further use in fingerprinter 
                    program
                    [0] CNP-ID
                    [1] InChi
                    [2] SMILES
                    [3] molecular formula
                    [4] NPLS
                [2] a csv with stats on errors
                    [0] (str) -- CNP-ID of entry causing InChi2mol errors
                    [1] (str) -- CNP-ID of entry causing SMILES2mol errors
                    [2] (list) -- CNP-ID of entries causing errors in both
                    [3] (int)  -- total number of entries

style:          attempting to follow PEP8 styleguide

"""
# --------------------------------- IMPORTS-------------------------------------
from sys import argv
from datetime import date
import os

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger  # for muting warnings

from biosynfoni.inoutput import readr, csv_writr, sdf_writr
from biosynfoni.inoutput import outfile_namer, output_direr


# ================================ FUNCTIONS ===================================

# ================================ converter ====================================


def propertifier(sdf_entries):  # 22 sec
    """lst of properties in entry as follows:
    every entry is [[coconut_id],[inchi],[smiles],[molecular_formula],[NPLS]]
    #for later:     [molecular weight], [name],]
    """
    entries_properties = []
    for entry in sdf_entries:
        # print(entry)
        i = 0
        coconut_id, inchi, smiles, molecular_formula = "", "", "", ""
        NPLS = ""
        for i in range(len(entry)):
            if entry[i] == "> <coconut_id>":
                coconut_id = entry[i + 1]
            elif entry[i] == "> <inchi>":
                inchi = entry[i + 1]
            elif entry[i] == "> <SMILES>":
                smiles = entry[i + 1]
            elif entry[i] == "> <molecular_formula>":
                molecular_formula = entry[i + 1]
            elif entry[i] == "> <NPL_score>":
                NPLS = entry[i + 1]
        entries_properties.append([coconut_id, inchi, smiles, molecular_formula, NPLS])

    return entries_properties


def molifier(NP_property_list, representation_info=False, stats_out=True):
    """makes mol class out of the NPs from database, using InChi formula,
    if InChi formula is not correct, it tries the SMILES. Keeps track of errors
    in InChi and SMILES, by adding the error-causing entries to the respective
    error-collection
    input:  (list) NP_property_list -- list of entries where
                    [0] CNP-ID
                    [1] InChi
                    [2] SMILES
                    [3] molecular formula
                    [4] NPLS (NaPLeS)
            (bool) stats_out = True -- yes or no return stats

    returns:    (list) all_mols -- mol classes of all successful inchi/smiles
                (list) prop_list --
                    [0] CNP-ID
                    [1] InChi
                    [2] SMILES
                    [3] molecular formula
                    [4] NPLS (NaPLeS)
                (list) stats -- if set to true, returns
                    [0] (str) -- entry CNP causing InChi2mol errors
                    [1] (str) -- entry CNP causing SMILES2mol errors
                    [2] (list) -- entries causing errors in both
                    [3] (int)  -- total number of entries
    """
    inchiNone = []
    smilesNone = []
    bothNone = []
    all_mols = []
    prop_list = []

    for entry in NP_property_list:
        smiles_mol = Chem.MolFromSmiles(entry[2])
        inchi_mol = Chem.MolFromInchi(entry[1])
        if inchi_mol is not None:
            all_mols.append(inchi_mol)
            if not representation_info:
                prop_list.append([entry[0]] + entry[3:])
                # avoiding faulty represent.
            elif representation_info:
                prop_list.append(entry)
        elif smiles_mol is not None:
            all_mols.append(smiles_mol)
            if not representation_info:
                prop_list.append([entry[0]] + entry[3:])
                # avoiding faulty represent.
            elif representation_info:
                prop_list.append(entry)
        else:  # these do not end up in the mol collection
            print("{} not successfully molified".format(entry[0]))
        # stats
        if stats_out:
            if inchi_mol is None:
                inchiNone.append(entry[0])
            if smiles_mol is None:
                smilesNone.append(entry[0])
            if (smiles_mol is None) and (inchi_mol is None):
                bothNone.append(entry[0])
    print("done with molifying  all entries")
    # output (yes or no stats)
    if stats_out:
        stats = [inchiNone, smilesNone, bothNone, len(NP_property_list)]
        return all_mols, prop_list, stats
    else:
        return all_mols, prop_list


# ================================= output =====================================


def sdf_writr(mols, outfile):
    """writes sdf of mols"""
    writer = Chem.SDWriter(outfile)
    for mol in mols:
        writer.write(mol)


# ++++++++++++++++++++++++++++++++++ main ++++++++++++++++++++++++++++++++++++++


def main():
    # input
    cdk_sdf = argv[1]

    print("\n", 10 * "=", "\n", "cOnVERTURE", "\n", 10 * "=", "\n")
    # handling
    lines = readr(cdk_sdf)
    sdf_entries = entry_parser(lines)
    properties = propertifier(sdf_entries)
    all_mols, ids_props, stats = molifier(properties, representation_info=True)

    # output
    actual_name = cdk_sdf.strip(".sdf").split("/")[-1]
    outfile_name = outfile_namer(actual_name)  # str

    out_rdk_sdf_name = "".join(outfile_name + "_rdk.sdf")
    out_info_name = "".join(outfile_name + "_info.tsv")
    out_stats_name = "".join(outfile_name + "_stats.csv")

    output_direr("./output")

    sdf_writr(all_mols, out_rdk_sdf_name)
    csv_writr(ids_props, out_info_name, sep="\t")
    csv_writr(stats, out_stats_name)

    return None


if __name__ == "__main__":
    main()
