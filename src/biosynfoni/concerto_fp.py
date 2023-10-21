#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: ConCERTO-fp                  ||
created: 2023-09-22                 ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
____________________________________
 
||||||||||||  ()()()  |||||||||||||||

ConCERTO-fp | CONvert Compound E-Representation TO FingerPrint
description: gets fp for given molecule, also shows which fp version was used
"""
from sys import argv
import argparse
import time

from rdkit import Chem

# my imports
from biosynfoni import def_biosynfoni
from biosynfoni.rdkfnx import get_supplier, get_subsset, save_version
from biosynfoni.inoutput import outfile_namer, csv_writr


# for the
DEFAULT_BIOSYNFONI_VERSION = def_biosynfoni.DEFAULT_BIOSYNFONI_VERSION
# =========================== per-mol functions ===============================


def cli():
    parser = argparse.ArgumentParser()

    # required
    parser.add_argument(
        "input",
        metavar="input: molecule(s) / molecule supplier",
        type=str,
        nargs="+",
        help=(
            "molecule representation to get biosynfoni of,"
            "can be sdf file, smiles string(s), or InChI string(s)."
            "will induce type, but recommended to add representation argument"
        ),
    )

    # optional
    parser.add_argument(
        "-r",
        "--repr",
        type=str,
        required=False,
        action="store",
        help="specify the input type of the molecule(s). If molsupplier, choose sdf",
        choices=["sdf", "smiles", "inchi"],
        default="induce",
    )

    parser.add_argument(
        "-v",
        "--version",
        type=str,
        required=False,
        action="store",
        help="specify the fingerprint version to use. If not specified, will use default",
        choices=def_biosynfoni.FP_VERSIONS.keys(),
        default=def_biosynfoni.DEFAULT_BIOSYNFONI_VERSION,
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        action="store",
        help="specify the output file name. If not specified, will use outfile_namer",
        default=None,
    )

    # flags
    parser.add_argument(
        "-l",
        "--interoverlap",
        "--lessblocking",
        required=False,
        action="store_true",
        help="allows overlap between different substructures",
        default=False,
    )
    parser.add_argument(
        "-n",
        "--intraoverlap",
        "--noblocking",
        "--noblockingstrong",
        required=False,
        action="store_true",
        help="allows all overlap: between different substructures and between same substructures",
        default=False,
    )

    parser.add_argument(
        "-c",
        "--coverage",
        required=False,
        action="store_true",
        help="pass if you want a file with the coverage data",
        default=False,
    )

    parser.add_argument(
        "-p",
        "--printonly",
        required=False,
        action="store_true",
        help=(
            "pass if you want to print the fingerprint to stdout."
            "will only print fingerprint, no extra data"
            "default for more than 1 molecule: False"
            "for 1 molecule will be overwritten to True, unless specified save"
        ),
        default=False,
    )
    parser.add_argument(
        "-s",
        "--save",
        required=False,
        action="store_true",
        help=("pass if you want to overwrite standard printonly for 1 mol"),
        default=None,
    )

    args = vars(parser.parse_args())

    # induce input type
    if args["repr"] == "induce":
        if args["input"][0].endswith(".sdf"):
            args["repr"] = "sdf"
        elif args["input"][0].startswith("InChI="):
            args["repr"] = "inchi"
        elif "." not in str(args["input"][0]):
            args["repr"] = "smiles"
        else:
            args["repr"] = "sdf"

    # overwrite for 1 mol if not specified save
    if (len(args["input"]) == 1) and (not args["save"]) and (not args["repr"] == "sdf"):
        args["printonly"] = True
    return args


def detect_substructures(
    mol: Chem.Mol,
    substructure_set: list,
    blocking_intersub: bool = True,
    blocking_intrasub: bool = True,
) -> list[tuple[tuple[int]]]:
    """for mol, extracts fingerprint & atomnumbers of each
    match

    input: (Chem.Mol) mol -- RDKit Chem Mol object, molecule
           (list) fp_set -- list of RDKit Chem mol files of
                            substructure keys
    output:
        [0] (list[tuples[int]]) atomnumbers of each match per
            substructure

    """
    # error handling
    if not mol:
        return []

    # output init.
    all_unique_matches = []  # empty at start

    # temporary----------------------------------------------------------
    blockedout_atoms = []

    for substructure in substructure_set:
        sub_matches = mol.GetSubstructMatches(substructure)  # all matches
        if not blocking_intersub:
            blockedout_atoms = []
        # start with all matches, then remove non-unique ones:
        unique_sub_matches = list(sub_matches)

        # removing overlapping finds from unique_matches
        for match in sub_matches:  # for each found match for curr subs
            approved_atoms = []  # unique matches for current subs
            if not blocking_intrasub:
                blockedout_atoms = []
            # yields either all atoms for a match or none (when has an overlap)
            for atom in match:
                if atom in blockedout_atoms:
                    # throw away match and all of its atoms, move onto next
                    unique_sub_matches.remove(match)
                    approved_atoms = []
                    break
                elif atom not in blockedout_atoms:
                    approved_atoms.append(atom)

            # if every atom in match is free, approve the match:
            if len(approved_atoms) > 0:  # non-thrown away match
                # valid match, so add to actual fp:
                for unique_atom in approved_atoms:
                    # to keep track of which atoms have been used
                    blockedout_atoms.append(unique_atom)

        all_unique_matches.append(tuple(unique_sub_matches))
    return all_unique_matches


def matches_to_vector(matches: list[tuple[tuple[int]]]) -> list[int]:
    """"""
    counted_vector = []
    for substructure_bit in matches:
        counted_vector.append(len(substructure_bit))
    return counted_vector


def get_biosynfoni(
    mol: Chem.Mol,
    version: str = "",
    substructure_set: list = [],
    return_matches: bool = False,
    blocking_intersub: bool = True,
    blocking_intrasub: bool = True,
) -> list[int]:
    """given a name of the fp version, uses get_subsset to get
    the set of substructures, passing them todetect_substructures
    to get the fingerprint. can also pass a substructure set directly,
    making it faster if repeated"""

    assert (
        version or substructure_set
    ), "please give either the version name or the substructure set"
    if not substructure_set:
        substructure_set = get_subsset(version)

    matches = detect_substructures(
        mol,
        substructure_set,
        blocking_intersub=blocking_intersub,
        blocking_intrasub=blocking_intrasub,
    )
    biosynfoni = matches_to_vector(matches)
    if return_matches:
        return biosynfoni, matches
    else:
        return biosynfoni


def subs_assigned_atoms(matches: list[tuple[tuple[int]]]) -> list[list[int]]:
    """per substructure, gives a list of the atom numbers that have been
    assigned to that substructure"""
    atoms_per_sub = []
    for matches_per_subs in matches:
        if not matches_per_subs:
            continue
        this_subs_atoms = []
        for match_in_subs in matches_per_subs:
            if isinstance(match_in_subs, int):
                # if the tuple tuple decides to unpack itself...
                this_subs_atoms.append(match_in_subs)
            elif len(match_in_subs) == 1:
                this_subs_atoms.append(match_in_subs[0])
            elif len(match_in_subs) > 1:
                for atom in match_in_subs:
                    this_subs_atoms.append(atom)
        # sub_atomcount = list(sum(sub, ()))
        atoms_per_sub.append(this_subs_atoms)
    return atoms_per_sub


def coverage_per_subs(matches: list[tuple[tuple[int]]]) -> list[float]:
    """returns number of atoms covered within molecule by each substructure"""
    coverage = []
    for sub in subs_assigned_atoms(matches):
        coverage.append(len(sub))

    assert len(coverage) == len(
        matches
    ), "error in coverage finding: coverage and matches lengths do not match"
    return coverage


def count_listitems(nested_list: list) -> int:
    count = 0
    for item in nested_list:
        if isinstance(item, list) or isinstance(item, tuple):
            sublist_count = count_listitems(list(item))  # recursive
            count = count + sublist_count
        elif item:
            count += 1
    return count


def get_coverage(mol: Chem.Mol, matches: list[tuple[tuple[int]]]) -> float:
    """gets non-h atom-based coverage of fingerprints over atoms"""
    matched_atoms = count_listitems(matches)
    nonh_mol = Chem.rdmolops.RemoveHs(mol, implicitOnly=False, sanitize=True)
    mol_size = nonh_mol.GetNumAtoms()  # non-H atoms
    coverage = float(matched_atoms) / float(mol_size)
    return coverage


def loop_over_supplier(
    supplier: Chem.SDMolSupplier,  # or list of mols
    fp_version: str = "",
    substructure_set: list = [],
    coverage_info=False,
    blocking_intersub=True,
    blocking_intrasub=True,
) -> list[list[int]]:
    """gets the fingerprint of the entire set"""
    # main functionality
    fingerprint_collection = []
    if coverage_info:
        coverage_collection = []
        for mol in supplier:
            fingerprint, matchlist = get_biosynfoni(
                mol,
                version=fp_version,
                substructure_set=substructure_set,
                return_matches=True,  # ==coverage_info
                blocking_intersub=blocking_intersub,
                blocking_intrasub=blocking_intrasub,
            )
            fingerprint_collection.append(fingerprint)
            coverage_collection.append(get_coverage(mol, matchlist))
        return fingerprint_collection, coverage_collection
    else:
        for mol in supplier:
            fingerprint = get_biosynfoni(
                mol,
                version=fp_version,
                substructure_set=substructure_set,
                return_matches=coverage_info,
                blocking_intersub=blocking_intersub,
                blocking_intrasub=blocking_intrasub,
            )
            fingerprint_collection.append(fingerprint)
        return fingerprint_collection


def handle_outnames(settings: dict) -> str:
    outname, inname_root, added = "", "", ""  # init
    if not settings["output"]:
        if settings["repr"] == "sdf":
            inname_root = settings["input"][0].split("/")[-1].split(".")[0]
        else:
            inname_root = f"{settings['repr']}s"
        if settings["intraoverlap"]:
            added = "_noblock"
        elif settings["interoverlap"]:
            added = "_lessblock"
        else:
            added = ""
        outname_root = f"{inname_root}_{settings['version']}{added}"
        outname = f"{outfile_namer(outname_root)}.bsf"

    else:
        if "." in settings["output"]:
            outname = settings["output"]
        else:
            outname = f"{settings['output']}.bsf"
    return outname


# ==========================  main ============================================


def main():
    settings = cli()
    inputlist = settings["input"]
    input_type = settings["repr"]
    fp_version = settings["version"]
    coverage_info = settings["coverage"]  # default False
    blocking_intersub = not settings["interoverlap"]  # default blocking == True
    blocking_intrasub = not settings["intraoverlap"]  # default blocking == True
    print_fp_only = settings["printonly"]  # default False
    outname = handle_outnames(settings)

    if not print_fp_only:
        print(10 * "=", "\nCONCERTO-FP\n", 10 * "=")
        print(settings)

    if input_type == "smiles":
        unchecked = [Chem.MolFromSmiles(mol) for mol in inputlist]
        inputs = [x for x in unchecked if x]
        errors = [inputlist[i] for i in range(len(unchecked)) if not unchecked[i]]
    elif input_type == "inchi":
        unchecked = [Chem.MolFromInchi(mol) for mol in inputlist]
        inputs = [x for x in unchecked if x]
        errors = [inputlist[i] for i in range(len(unchecked)) if not unchecked[i]]
    elif input_type == "sdf":
        if not print_fp_only:
            print("getting supplier...")
        inputs = get_supplier(
            settings["input"][0]
        )  # no error check due to time constraints
        errors = []
    if not print_fp_only:
        if len(inputs) == 0:
            print("no valid molecules found")
        print("looping over molecules...")

    # to prevent not knowing which result is which, include 'None' type
    if print_fp_only and input_type != "sdf":
        inputs = unchecked

    loopres = loop_over_supplier(
        inputs,
        substructure_set=get_subsset(fp_version),
        coverage_info=coverage_info,
        blocking_intersub=blocking_intersub,
        blocking_intrasub=blocking_intrasub,
    )

    if coverage_info:
        biosynfonies, coverages = loopres
        if not print_fp_only:
            print("writing coverages...")
            csv_writr(
                coverages, f"{outname.replace('.bsf','')}_coverages.tsv", sep="\t"
            )
    else:
        biosynfonies = loopres

    if print_fp_only:
        # print("==============")  # for parsing
        # print(f"{fp_version}")
        # print("--------------")  # for parsing
        for row in biosynfonies:
            if not row:
                print("None")
            else:
                print(*row, sep=",")
    else:
        print(f"writing {len(biosynfonies)} biosynfonies to file...")
        csv_writr(biosynfonies, outname, sep=",")
        if input_type != "sdf":
            with open(f"{outname}_input.txt", "w") as f:
                f.write("\n".join([x for x in inputlist if x not in errors]))
        save_version(fp_version, extra_text="extracted")
        print("done")

    return None


if __name__ == "__main__":
    main()
