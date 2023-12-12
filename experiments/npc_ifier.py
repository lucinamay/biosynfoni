"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: NPC_ifier                    ||
language: python                    ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
student no: 1197142                 ||
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

# /////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////

# |||||||||||||||||||||||||||||||| npc_ifier ||||||||||||||||||||||||||||||||

# /////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////

# --------------------------------- IMPORTS-------------------------------------
import argparse, os
import requests, json
import time
import urllib
from datetime import date


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger  # for muting warnings
from tqdm import tqdm
import numpy as np

# RDLogger.DisableLog('rdApp.*') # mutes warnings

# ------------------------------------------------------------------------------

# ================================= GLOBALS ====================================

# no globals (yet)

# ==============================================================================

# ================================= CLASSES ====================================

# none yet

# ==============================================================================


# /////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////


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


# def mols_getr(sdf_file: str, supplier_only=False):
#     suppl = Chem.SDMolSupplier(sdf_file)
#     print("read sdf, number of entries:", len(suppl))
#     start_time = time.time()
#     if supplier_only == False:
#         mols = [x for x in suppl]
#     end_time = time.time()
#     print("mols_getr took {} seconds to finish".format((end_time - start_time)))
#     for mol in mols:
#         if mol is None:
#             raise Exception(
#                 print(
#                     mol,
#                     "mol is None; check if you are \
#                     using RDKit-based sdf",
#                 )
#             )
#     return mols


# # ================================= handling ===================================


# def save_index_getr(rootname: str, save_size=1000, max_index=408):
#     """checks for temporary file lists, returns (list) indexes(int)
#     still needed in folder
#     filenames are <dateofnpcadd>_<rootname>_npcs_index.tsv
#     """
#     file_list = os.listdir(".")
#     file_roots = [name.strip(".csv").strip(".tsv") for name in file_list]
#     indexes_left = []  # all indexes that are missing

#     for i in range(max_index + 1):  # includes last index
#         # reset per  index
#         found = False

#         # go through list until index found
#         for fi in file_roots:
#             if rootname in fi:
#                 if fi.split("_")[-3:] == [str(save_size), "npcs", str(i)]:
#                     found = True
#                     break
#         if found == True:
#             continue  # skip to next index
#         else:
#             indexes_left.append(int(i))
#     return indexes_left


# def slicer(full_lst: list, slice_ind: int, slice_size: int):
#     # num_slices = (total // slice_size) + 1

#     i = slice_ind
#     slice_start, slice_end = i * slice_size, (i * slice_size + slice_size)

#     return full_lst[slice_start, slice_end]


# --------------------------- get classification -------------------------------

# def npc_apier(mol_list: list, index_base=0, chunk_size=5, sleep_time=5):
#     print("npc_apier running starting at index {}".format(index_base))

#     result = []  # array with rows == molecule
#     rest_count = 0

#     for i in range(len(mol_list)):
#         if rest_count == chunk_size:
#             print("at {}, now sleep".format(i + index_base))
#             time.sleep(sleep_time)
#             rest_count = 0
#         mol = mol_list[i]
#         smiles = Chem.MolToSmiles(mol)  # gives error if there is a # or +
#         try:
#             npc_result = npc_resultr(smiles)
#             print(".", end="")
#         except:
#             print("e{}".format(i), end="")
#             errorfile = outfile_namer("api_errors")
#             with open(errorfile, "a") as errors:
#                 errors.write("\t".join(str(j) for j in [index_base + i, smiles, "\n"]))
#             npc_result = [[], [], []]
#         result.append(npc_result)

#         rest_count += 1

#     return result


def npc_resultr(smiles_string: str):
    """requests the classification of a smiles string"""
    url = "https://npclassifier.ucsd.edu"

    r = requests.get(f"{url}/classify?smiles={smiles_string}")
    # npc_pathway = r.json()["pathway_results"]
    # npc_superclass = r.json()["superclass_results"]
    # npc_class = r.json()["class_results"]
    # return [npc_pathway, npc_superclass, npc_class]
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


# # ==============================================================================

# # ================================= output =====================================


# def outfile_namer(filename_root: str, addition=""):
#     """gives outfile names in the right formats (no extension names)"""
#     today = date.today()
#     my_date = today.strftime("%m%d")
#     if len(addition) > 0:
#         outfile_name = "{}_{}_{}".format(my_date, filename_root, addition)
#     else:
#         outfile_name = "{}_{}".format(my_date, filename_root)

#     return outfile_name


# def nested_list_tsvr(nested_list, filename, sep_large="\t", sep_small=","):
#     totname = filename + ".tsv"
#     with open(totname, "w") as new:
#         for row in nested_list:
#             new.write(sep_large.join(sep_small.join(thing) for thing in row))
#             new.write("\n")
#     print("written all {} lines in file: {}".format(len(nested_list), totname))
#     return None


# # ==============================================================================
# # ============================ handling-output loop ============================


# def npc_over_set(mols: list, rootname: str, savesize=1000, chunksize=5, sleep_time=5):
#     """divides the api requests over 'save chuncks' (i.e. number of molecules
#     after which the api result is being written) (default: 1000) and request
#     chuncks (i.e. number of requests after which the requests 'sleep' for a
#     while
#     """

#     print("starting to loop over full set of {} molecules".format(len(mols)))

#     # move to the right folder (hardcoded) ---------------------------------
#     storage = "./npc_storage"
#     if not os.path.exists(storage):
#         os.mkdir(storage)
#     os.chdir(storage)

#     # save-chunk: every <savesize> molecules (default 1000) ----------------
#     num_of_files = (len(mols) // savesize) + 1
#     print(num_of_files, "should be created in total")

#     indexes_left = save_index_getr(rootname, savesize, num_of_files)
#     # what savechunck
#     num_left = len(indexes_left)
#     diff = num_of_files - num_left
#     print(num_left, "are not yet done")
#     print(
#         "[{}{}]{}%".format((diff * "="), (num_left * " "), (diff / num_of_files * 100))
#     )

#     for i in indexes_left:
#         # get slice to save to one savefile (not necessarily start with 0)
#         slice_start, slice_end = i * savesize, (i * savesize + savesize)
#         savechunk = mols[slice_start:slice_end]

#         # request entire slice smartly - - - - - - - - - - - - - - - - - -
#         npc_result = npc_apier(
#             savechunk,
#             index_base=slice_start,
#             chunk_size=chunksize,
#             sleep_time=sleep_time,
#         )
#         # returns nested list:  [ [[],[],[]], [[],[],[]], [[],[],[]] ]

#         # after completing slice, write  - - - - - - - - - - - - - - - - -
#         name_addition = "{}_npcs_{}".format(savesize, i)
#         newfile_name = outfile_namer(rootname, name_addition)
#         nested_list_tsvr(npc_result, newfile_name)

#     os.chdir("..")  # return
#     return None


# # ----------------------------- afterloop combining ----------------------------
# def output_combinr(rootname, savesize=1000):
#     """combines all the files together"""
#     print("starting output combination ...")

#     storage = "./npc_storage"
#     os.chdir(storage)

#     file_list = os.listdir(".")
#     file_roots = [name.strip(".csv").strip(".tsv") for name in file_list]
#     query = "{}_npcs".format(savesize)

#     new_filename = outfile_namer(rootname, "npcs") + ".tsv"

#     if os.path.exists(new_filename):
#         raise Exception("\n\n{} exists, delete file first".format(new_filename))

#     relevant = []
#     for fi in file_list:
#         if (rootname in fi) and (query in fi):
#             relevant.append(fi)

#     relevant_files = sorted(relevant, key=lambda x: int(x.split("_")[-1]))
#     indexes_done = []

#     for eachfile in relevant_files:
#         index = eachfile.strip(".tsv").split("_")[-1]
#         if index in indexes_done:
#             continue
#         indexes_done.append(index)
#         print("copying {} into {}".format(eachfile, new_filename))
#         temp = []
#         # read in small file
#         with open(eachfile) as readfile:
#             for line in readfile:
#                 temp.append(line.strip())
#         # write to big file
#         with open(new_filename, "a") as total:
#             for line in temp:
#                 total.write(line)
#                 total.write("\n")
#         # write log of copy
#         with open(outfile_namer(rootname, "log"), "a") as log:
#             log.write(
#                 "written {} lines of {} to {}".format(len(temp), fi, new_filename)
#             )
#     return None


# ==============================================================================


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

    # # getmolset
    # mols = mols_getr(sdf_file, supplier_only=False)
    # print("preview of mols", mols[:10], "\n")

    # # output
    # os.chdir("./output")
    # npc_over_set(mols, sdf_name, savesize=savesize, chunksize=5, sleep_time=5)

    # # output combining
    # output_combinr(sdf_name, savesize)

    print("\n =======\n DONE~~~\n =======\n")
    os.chdir(start_dir)
    exit(0)
    return None


if __name__ == "__main__":
    main()
