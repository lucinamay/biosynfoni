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
"""
import argparse
from itertools import chain
from functools import partial

import numpy as np
from rdkit import Chem
from tqdm import tqdm


# my imports
from biosynfoni import def_biosynfoni
from biosynfoni.rdkfnx import (
    supplier,
    get_subs_set,
    save_version,
    nonh_atomcount,
    get_sub_matches,
    smiles_to_mol,
    inchi_to_mol,
)
from biosynfoni.inoutput import outfile_namer, csv_writr, myprint


DEFAULT_BIOSYNFONI_VERSION = def_biosynfoni.DEFAULT_BIOSYNFONI_VERSION


def _induce_input_type(inputlist: list) -> str:
    """induces the input type from the inputlist"""
    if inputlist[0].endswith(".sdf"):
        return "sdf"
    elif inputlist[0].startswith("InChI="):
        return "inchi"
    elif "." not in str(inputlist[0]):
        return "smiles"
    else:
        return "sdf"


def _handle_outnames(args: argparse.Namespace) -> str:
    outname, inname_root, overlap_allowance = "", "", ""  # init
    if args.output is None:
        if args.repr == "sdf":
            inname_root = args.input[0].split("/")[-1].split(".")[0]
        else:
            inname_root = f"{args.repr}s"
        if args.intraoverlap:
            overlap_allowance = "_noblock"
        elif args.interoverlap:
            overlap_allowance = "_lessblock"
        else:
            overlap_allowance = ""
        outname_root = f"{inname_root}_{args.version}{overlap_allowance}"
        outname = f"{outfile_namer(outname_root)}.bsf"

    else:
        if "." in args.output:
            outname = args.output
        else:
            outname = f"{args.output}.bsf"
    return outname


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

    # args = vars(parser.parse_args())
    args = parser.parse_args()

    if args.repr == "induce":
        args.repr = _induce_input_type(args.input)

    # overwrite for 1 mol if not specified save
    if args.save is None:
        if len(args.input) == 1 and not args.repr == "sdf":
            args.printonly = True

    args.output = _handle_outnames(args)

    return args


# ============================= MAIN FUNCTIONS  =================================
# -------------------------- per-mol functions ---------------------------------
def _count_listitems(nested_list: list) -> int:
    count = 0
    for item in nested_list:
        if isinstance(item, list) or isinstance(item, tuple):
            sublist_count = _count_listitems(list(item))  # recursive
            count = count + sublist_count
        elif item:
            count += 1
    return count


class MolFp:
    def __init__(
        self,
        mol: Chem.Mol,
        substructure_set: list = None,
        version_name: str = DEFAULT_BIOSYNFONI_VERSION,
        intersub_overlap: bool = False,
        intrasub_overlap: bool = False,
    ) -> None:
        self.mol = mol
        self.version = version_name
        print(substructure_set)
        self.substructure_set = self.set_substructure_set(substructure_set)
        self.intersub_overlap = intersub_overlap
        self.intrasub_overlap = intrasub_overlap
        self.overlap_allowance = {
            "intersub": self.intersub_overlap,
            "intrasub": self.intrasub_overlap,
        }
        self.matches = self._detect_substructures()
        self.biosynfoni = self.get_biosynfoni()
        self.coverage = self._matches_to_coverage()

    def set_substructure_set(self, substructure_set) -> None:
        if substructure_set is None:
            return get_subs_set(self.version)
        else:
            return substructure_set

    # def set_mol(self, mol: Chem.Mol) -> None:
    #     self.mol = mol

    # def set_substructure_set(self, substructure_set: list) -> None:
    #     self.substructure_set = substructure_set

    # def set_overlap_allowance_intersub(self, intersub_overlap: bool) -> None:
    #     self.overlap_allowance["intersub"] = intersub_overlap

    # def set_overlap_allowance_intersub(self, intrasub_overlap: bool) -> None:
    #     self.overlap_allowance["intrasub"] = intrasub_overlap

    # def set_overlap_allowance(self, overlap_allowance: dict) -> None:
    #     self.overlap_allowance = overlap_allowance

    # def set_matches(self, matches: list[tuple[tuple[int]]]) -> None:
    #     self.matches = matches

    # def set_biosynfoni(self, biosynfoni: list[int]) -> None:
    #     self.biosynfoni = biosynfoni

    # def set_coverage(self, coverage: float) -> None:
    #     self.coverage = coverage

    # def get_mol(self) -> Chem.Mol:
    #     return self.mol

    # def get_substructure_set(self) -> list:
    #     return self.substructure_set

    # def get_overlap_allowance(self) -> dict:
    #     return self.overlap_allowance

    # def get_matches(self) -> list[tuple[tuple[int]]]:
    #     return self.matches

    # def get_biosynfoni(self) -> list[int]:
    #     return self.biosynfoni

    # def get_coverage(self) -> float:
    #     return self.coverage

    def _sub_matches_to_atoms(self, sub_matches: list[list[int]]) -> list[int]:
        # return list(set(chain.from_iterable(sub_matches)))
        return list(chain.from_iterable(sub_matches))

    def __filter_intersub_overlap(
        self, sub_matches: list[list[int]], previous_matches: list[list[int]]
    ) -> list[list[int]]:
        """filters out matches that have overlap with previous matches"""
        filtered_matches = []
        for sub_match in sub_matches:
            overlap = False
            for atom in sub_match:
                if atom in self._sub_matches_to_atoms(previous_matches):
                    overlap = True
                    break
            if not overlap:
                filtered_matches.append(sub_match)
        return filtered_matches

    def __filter_intrasub_overlap(
        self, sub_matches: list[list[int]]
    ) -> list[list[int]]:
        """filters out matches that have overlap within the same substructure"""
        filtered_matches = []
        for sub_match in sub_matches:
            overlap = False
            for atom in sub_match:
                if atom in self._sub_matches_to_atoms(filtered_matches):
                    overlap = True
                    break
            if not overlap:
                filtered_matches.append(sub_match)
        return filtered_matches

    def _detect_substructures(self) -> list[list[list[int]]]:
        """for mol, extracts fingerprint & atomnumbers of each
        match
        """
        # output init.
        all_approved_matches = []  # empty at start
        for substructure in self.substructure_set:
            sub_matches = get_sub_matches(self.mol, substructure)  # any matches
            if not self.intersub_overlap:
                sub_matches = self.__filter_intersub_overlap(
                    sub_matches, all_approved_matches
                )
            if not self.intrasub_overlap:
                sub_matches = self.__filter_intrasub_overlap(sub_matches)
            all_approved_matches.append(sub_matches)
        return all_approved_matches

    def _matches_to_substruct_counts(self) -> list[int]:
        """converts the arrays of atoms of matches to an array with the number of matches per substructure"""
        substruct_counts = []
        for matches_per_sub in self.matches:
            substruct_counts.append(len(matches_per_sub))
        return substruct_counts

    def get_biosynfoni(self) -> list[int]:
        """given a name of the fp version, uses get_subs_set to get
        the set of substructures, passing them todetect_substructures
        to get the fingerprint. can also pass a substructure set directly,
        making it faster if repeated"""

        # error handling
        if self.mol is None:
            return []
        else:
            return self._matches_to_substruct_counts()

    def _matched_atoms_per_substructure(self) -> list[list[int]]:
        """per substructure, gives a list of the atom numbers that have been
        assigned to that substructure"""
        atoms_per_sub = []
        for sub_matches in self.matches:
            sub_matches_atoms = self._sub_matches_to_atoms(sub_matches)
            atoms_per_sub.append(sub_matches_atoms)
        return atoms_per_sub

    def _matches_to_coverage(self) -> float:
        """gets non-h atom-based coverage of fingerprints over atoms"""
        matched_atoms = _count_listitems(self.matches)
        coverage = float(matched_atoms) / float(nonh_atomcount(self.mol))
        return coverage

    def _matches_to_coverage_per_substructure(self) -> list[float]:
        """gets non-h atom-based coverage of fingerprints over atoms"""
        mol_size = nonh_atomcount(self.mol)
        atoms_per_sub = self._matched_atoms_per_substructure()
        coverage_per_sub = [
            float(len(atoms)) / float(mol_size) for atoms in atoms_per_sub
        ]
        return coverage_per_sub


# ------------------------------ multiple mols ----------------------------------


# def loop_over_supplier(
#     supplier: Chem.SDMolSupplier,  # or list of mols
#     fp_version: str = "",
#     substructure_set: list = [],
#     coverage_info=False,
#     blocking_intersub=True,
#     blocking_intrasub=True,
# ) -> list[list[int]]:
#     """gets the fingerprint of the entire set"""
#     # main functionality
#     fingerprint_collection = []
#     if coverage_info:
#         coverage_collection = []
#         for mol in supplier:
#             fingerprint, matchlist = get_biosynfoni(
#                 mol,
#                 version=fp_version,
#                 substructure_set=substructure_set,
#                 return_matches=True,  # ==coverage_info
#                 blocking_intersub=blocking_intersub,
#                 blocking_intrasub=blocking_intrasub,
#             )
#             fingerprint_collection.append(fingerprint)
#             coverage_collection.append(_matches_to_coverage(mol, matchlist))
#         return fingerprint_collection, coverage_collection
#     else:
#         for mol in supplier:
#             fingerprint = get_biosynfoni(
#                 mol,
#                 version=fp_version,
#                 substructure_set=substructure_set,
#                 return_matches=coverage_info,
#                 blocking_intersub=blocking_intersub,
#                 blocking_intrasub=blocking_intrasub,
#             )
#             fingerprint_collection.append(fingerprint)
#         return fingerprint_collection


def _yield_sdf_mols(sdflist: list) -> Chem.Mol:
    """yields mols from supplier (only first one of list)"""
    suppl = supplier(sdflist[0])
    for mol in suppl:
        yield mol


def _yield_smiles_mols(smileslist: list) -> Chem.Mol:
    """yields mols from smileslist"""
    for mol in smileslist:
        yield smiles_to_mol(mol)


def _yield_inchi_mols(inchilist: list) -> Chem.Mol:
    """yields mols from inchilist"""
    for mol in inchilist:
        yield inchi_to_mol(mol)


def _get_yield_mol_function(input_type: str) -> Chem.Mol:
    """yields mols from inputlist"""
    types_functions = {
        "sdf": _yield_sdf_mols,
        "smiles": _yield_smiles_mols,
        "inchi": _yield_inchi_mols,
    }
    return types_functions[input_type]


# ==============================  MAIN =======================================


def main():
    args = cli()
    print_only_fp = args.printonly  # default False

    # set printing on and off depending on print_only_fp
    print_ = partial(myprint, print_only_fp)

    print_(10 * "=", "\nCONCERTO-FP\n", 10 * "=")
    print_(args)

    # main functionality
    # get substructure set to speed up process (and not have to get it each time)
    substructure_set = get_subs_set(args.version)

    total_molecules = (
        len(args.input) if args.repr != "sdf" else len(supplier(args.input[0]))
    )
    yield_mol_function = _get_yield_mol_function(args.repr)

    print_("looping over molecules...")
    biosynfonies, coverages = [], []
    for i, mol in tqdm(
        enumerate(yield_mol_function(args.input)), total=total_molecules
    ):
        if mol is None:
            with open("biosynfoni.err", "a") as ef:
                ef.write(f"{i}\n")

        mol_fp = MolFp(
            mol,
            substructure_set,
            args.version,
            args.interoverlap,
            args.intraoverlap,
        )
        if args.coverage:
            coverages.append(mol_fp.coverage)
        biosynfonies.append(mol_fp.biosynfoni)

    if print_only_fp:
        for row in biosynfonies:
            if not row:
                print_("None")
            else:
                print(*row, sep=",")
        exit(0)
        return None
    if args.coverage:
        print_("writing coverages...")
        csv_writr(coverages, args.output.replace(".bsf", "_coverages.tsv"), sep="\t")

    print_(f"writing {len(biosynfonies)} biosynfonies to {args.output}...")

    biosynfonies_array = np.array(biosynfonies)
    np.savetxt(fname=args.output, X=biosynfonies_array, fmt="%i", delimiter=",")

    if args.repr != "sdf":
        with open(args.output.replace(".bsf", "_input.tsv"), "w") as f:
            f.write("\n".join([x for x in args.input]))
    save_version(args.version, extra_text="extracted")
    print_("done")
    exit(0)
    return None


if __name__ == "__main__":
    main()
