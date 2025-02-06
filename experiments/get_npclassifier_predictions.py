"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: NPC_ifier                    ||
language: python                    ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

description:    using RDK Chem mols and NPClassifier API, requests
                classification information with time delay
                saves the classfications every 1000 molecules
                at the end, takes the 400 files and combines them
                
                *checks which files are present in the related folder
                and starts from where it left off (rewriting only the last file
                in case it was not completed with a 1000 molecules)

                *uses mol to smiles to make sure the smiles is accurate to the
                mol files used in the rest of the experiment (fp classif.

"""

# --------------------------------- IMPORTS-------------------------------------
import argparse, os
import json
import time
import urllib
from datetime import date


import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger  # for muting warnings
from tqdm import tqdm
import numpy as np


# ================================ FUNCTIONS ===================================
# ================================== input =====================================
def cli():
    """command line interface"""
    parser = argparse.ArgumentParser(
        description="npc_ifier: using RDK Chem mols and NPClassifier API, requests \
        classification information with time delay \
        saves the classfications every 1000 molecules \
        at the end, takes the 400 files and combines them"
    )
    parser.add_argument(
        "sdf_file",
        metavar="sdf_file",
        type=str,
        help="sdf file to be used for the experiment",
    )
    parser.add_argument(
        "-s",
        "--savesize",
        type=int,
        default=100,
        help="number of molecules after which the api result is being written",
    )
    parser.add_argument(
        "-c",
        "--chunksize",
        type=int,
        default=5,
        help="number of requests after which the requests 'sleep' for a while",
    )
    parser.add_argument(
        "-t",
        "--sleeptime",
        type=int,
        default=5,
        help="number of seconds to sleep after chunksize requests",
    )
    args = parser.parse_args()
    return args


def npc_resultr(smiles_string: str):
    """requests the classification of a smiles string"""
    url = "https://npclassifier.ucsd.edu"

    r = requests.get(f"{url}/classify?smiles={smiles_string}")
    return r.json()


def mol_to_result(mol: Chem.Mol) -> dict[str]:
    """takes a list of molecules and returns a dictionary of classifications,
    trying again with an encoded query or -1 if there is an error"""
    smiles = Chem.MolToSmiles(mol)  # gives error if there is a # or +
    try:
        npc_result = npc_resultr(smiles)
    except:
        try:
            encoded_smiles = urllib.parse.quote(smiles)
            npc_result = npc_resultr(encoded_smiles)
        except:
            npc_result = -1
    return npc_result


def json_to_classification(js: dict) -> list[str]:
    """takes a json and returns a list of classifications"""
    npc_pathway = ",".join(js["pathway_results"])
    npc_superclass = ",".join(js["superclass_results"])
    npc_class = ",".join(js["class_results"])
    return [npc_pathway, npc_superclass, npc_class]


def json_to_fp(js: dict) -> tuple[np.array]:
    """takes a json and returns a tuple of fp arrays"""
    fp1 = np.array(js["fp1"], dtype=int)
    fp2 = np.array(js["fp1"], dtype=int)
    return fp1, fp2


def handle_npc_result(npc_result: dict, index: int) -> None:
    """writes the npc results to the right files"""
    index_predictions = [str(index)] + json_to_classification(npc_result)
    preds = np.array(index_predictions, dtype=str)

    with open("npc_predictions.tsv", "a") as pr_f:
        np.savetxt(pr_f, preds, fmt="%s", newline="\t")
        pr_f.write("\n")

    try:
        fp1, fp2 = json_to_fp(npc_result)
        fp1 = np.array(([index] + list(fp1)), dtype=int)
        fp2 = np.array(([index] + list(fp2)), dtype=int)
        with open("npc_fp1.tsv", "a") as fp1_f:
            np.savetxt(fp1_f, fp1, fmt="%s", newline="\t")
            fp1_f.write("\n")
        with open("npc_fp2.tsv", "a") as fp2_f:
            np.savetxt(fp2_f, fp2, fmt="%s", newline="\t")
            fp2_f.write("\n")
    except:
        with open("npc_fp_errors.txt", "a") as errors:
            errors.write(f"{index}\n")
    return None


# /////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////
# ++++++++++++++++++++++++++++++++++ main ++++++++++++++++++++++++++++++++++++++
# /////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////


def main():
    args = cli()

    savesize = args.savesize
    chunksize = args.chunksize
    sleep_time = args.sleeptime

    # get absolute path
    sdf_path = os.path.abspath(args.sdf_file)
    sdf_name = sdf_path.split("/")[-1].replace(".sdf", "")

    # make output folder
    start_dir = os.getcwd()
    if not os.path.exists("./npc_ifier/"):
        os.mkdir("./npc_ifier/")
    os.chdir("./npc_ifier/")
    if not os.path.exists(sdf_name):
        os.mkdir(sdf_name)
    os.chdir(sdf_name)

    # print start

    print("\n", 10 * "=", "\n", "NPC_ifier\n", date.today(), "\n", 10 * "=", "\n")

    last_index = -1
    if os.path.exists("npc_results.tsv"):
        indexes = np.loadtxt("npc_results.tsv", dtype=int, delimiter="\t", usecols=(0))
        print("npc_results.tsv exists, starting from last index")
        last_index = indexes[-1]

    suppl = Chem.SDMolSupplier(sdf_path)
    for index, mol in tqdm(enumerate(suppl), total=len(suppl)):
        if index <= last_index:
            continue

        if index % chunksize == 0:
            time.sleep(sleep_time)

        if mol is None:
            with open("mol_errors.txt", "a") as errors:
                errors.write(f"{index}\n")
            continue

        smiles = Chem.MolToSmiles(mol)
        npc_result = mol_to_result(mol)

        if npc_result == -1:
            with open("npc_errors.txt", "a") as errors:
                errors.write(f"{index}\t{smiles}\n")
            continue

        handle_npc_result(npc_result, index)

    print("\n =======\n DONE~~~\n =======\n")
    os.chdir(start_dir)
    exit(0)
    return None


if __name__ == "__main__":
    main()
