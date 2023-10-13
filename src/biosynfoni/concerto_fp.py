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

from rdkit import Chem

# my imports
import def_biosynfoni
from rdkfnx import get_supplier, get_subsset
from inoutput import outfile_namer, csv_writr, save_version


# not used, just for reference
DEFAULT_BIOSYNFONI_VERSION = def_biosynfoni.DEFAULT_BIOSYNFONI_VERSION
# =========================== per-mol functions ===============================


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
) -> list[int]:
    """given a name of the fp version, uses get_subsset to get
    the set of substructures, passing them todetect_substructures
    to get the fingerprint. can also pass a substructure set directly"""

    assert (
        version or substructure_set
    ), "please give either the version name or the substructure set"
    if not substructure_set:
        substructure_set = get_subsset(version)

    matches = detect_substructures(mol, substructure_set)
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
            print(match_in_subs)
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
    supplier: Chem.SDMolSupplier,
    fp_version: str = "",
    substructure_set: list = [],
    coverage_info=False,
):
    """gets the fingerprint of the entire set"""

    print("looping over supplier, coverage info provided:", coverage_info)

    # main functionality
    fingerprint_collection = []
    if coverage_info:
        coverage_collection = []
        for mol in supplier:
            fingerprint, matchlist = get_biosynfoni(
                mol,
                version=fp_version,
                substructure_set=substructure_set,
                return_matches=coverage_info,
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
            )
            fingerprint_collection.append(fingerprint)
        return fingerprint_collection


# ==========================  main ============================================


def main():
    print(10 * "=", "\nCONCERTO-FP\n", 10 * "=")

    coverage_info = False  # default
    blocking = True  # default

    supplier_loc = argv[1]
    assert supplier_loc.split(".")[-1] == "sdf", f"non-sdf file: {supplier_loc}"
    fp_version = argv[2]
    assert (
        fp_version in def_biosynfoni.FP_VERSIONS.keys()
    ), f"invalid fp version: {fp_version}"

    if len(argv) >= 4:
        coverage_info_bool = argv[3]  # if want coverage: write True or T
        if coverage_info_bool in [
            "True",
            "true",
            "T",
            "t",
            "coverage",
            "y",
            "yes",
            "Y",
            "Yes",
        ]:
            coverage_info = True
    if len(argv) >= 5:
        blocking = argv[5]  # if want coverage: write True or T
        if blocking in ["no", "NO", "N", "n", "false", "f", "No", "False", "noblock"]:
            blocking = False
            print(
                "Blocking out set to False, will detect all overlapping",
                "building blocks",
            )

    outname = outfile_namer(supplier_loc, "bsf")
    if not blocking:
        outname = outfile_namer(supplier_loc, "overlap_bsf")

    print("getting supplier...")
    supplier = rdkfnx.get_supplier(supplier_loc)
    print("looping over supplier...")

    if not coverage_info:
        biosynfonies = loop_over_supplier(
            supplier,
            substructure_set=get_subsset(fp_version),
            coverage_info=coverage_info,
        )
    else:
        biosynfonies, coverages = loop_over_supplier(
            supplier,
            substructure_set=get_subsset(fp_version),
            coverage_info=coverage_info,
        )
        csv_writr(coverages, f"{outname}_coverages.tsv", sep="\t")

    print(f"writing {len(biosynfonies)} biosynfonies to file...")
    csv_writr(biosynfonies, f"{outname}.bsf", sep=",")
    save_version(fp_version, extra_text="inbsf")
    print("done")

    return None


if __name__ == "__main__":
    main()
