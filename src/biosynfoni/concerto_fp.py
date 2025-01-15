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
import logging
from itertools import chain
from functools import partial


import numpy as np
from rdkit import Chem
from tqdm import tqdm


# my imports
from biosynfoni.cli import cli
from biosynfoni.rdkfnx import (
    supplier,
    save_version,
    nonh_atomcount,
    get_sub_matches,
    smiles_to_mol,
    inchi_to_mol,
)
from biosynfoni.rdkfnx import BiosynfoniVersion
from biosynfoni.inoutput import csv_writr, myprint
from biosynfoni.subkeys import defaultVersion


# =============================  HELPERS  =====================================
# -------------------------- general helpers ----------------------------------
def _count_listitems(nested_list: list) -> int:
    count = 0
    for item in nested_list:
        if isinstance(item, list) or isinstance(item, tuple):
            sublist_count = _count_listitems(list(item))  # recursive
            count = count + sublist_count
        elif item or item == 0:
            count += 1
    return count


def _list_unnest_once(once_nested: list[list[int]]) -> list[int]:
    # return list(set(chain.from_iterable(sub_matches)))
    return list(chain.from_iterable(once_nested))


def _list_unnest_twice(twice_nested: list[list[list[int]]]) -> list[int]:
    return list(chain.from_iterable(chain.from_iterable(twice_nested)))


# ============================= MAIN FUNCTIONS  =================================
# -------------------------- per-mol functions ---------------------------------


class Biosynfoni:
    def __init__(
        self,
        mol: Chem.Mol,
        substructure_set: list = None,
        version_name: str = defaultVersion,
        intersub_overlap: bool = True,
        intrasub_overlap: bool = True,
    ) -> None:
        self.mol = mol
        self.version = version_name
        self.substructure_set = self._set_substructure_set(substructure_set)
        self.intersub_overlap = intersub_overlap
        self.intrasub_overlap = intrasub_overlap
        self.overlap_allowance = {
            "intersub": self.intersub_overlap,
            "intrasub": self.intrasub_overlap,
        }
        self.overlap_matches = self._detect_substructures()
        # self.inter_overlap_matches = self._get_nonoverlaps()
        self.matches = self._filter_matches()
        self.counted_fp = self.get_biosynfoni()
        self.fingerprint = self.counted_fp
        # self.coverage = self._matches_to_coverage()
        self.coverage = None
        self.coverage_per_sub = None

    def _set_substructure_set(self, substructure_set) -> None:
        if substructure_set is None:
            return BiosynfoniVersion(self.version).substructures
        else:
            return substructure_set

    def _detect_substructures(self) -> list[list[list[int]]]:
        """for mol, extracts fingerprint & atomnumbers of each
        match
        """
        # output init.
        all_detected_matches = []  # empty at start
        for substructure in self.substructure_set:
            sub_matches = get_sub_matches(self.mol, substructure)  # any matches
            all_detected_matches.append(sub_matches)
        return all_detected_matches

    def __filter_intrasub_overlap(
        self, sub_matches: list[list[int]]
    ) -> list[list[int]]:
        """
        filters out matches that have overlap within the same substructure
        (i.e. one atom cannot belong to two 'versions' of the same substructure. e.g. for CCCC there can only be two matches for CC)
        """
        # per substructure
        filtered_matches = []
        # sort sub_matches from largest to smallest ("heuristic set coverage")
        sub_matches.sort(key=len, reverse=True)

        for sub_match in sub_matches:
            overlap = False  # reset
            atoms_to_block = _list_unnest_once(filtered_matches)
            for atom in sub_match:
                if atom in atoms_to_block:
                    overlap = True
                    break
            if not overlap:
                filtered_matches.append(sub_match)
        return filtered_matches

    def __filter_intra_smarter(self, sub_matches: list[list[int]]) -> list[list[int]]:
        filtered_matches = []
        all_atoms = _list_unnest_once(sub_matches)
        # sort from small to large index num
        all_atoms.sort()
        for atom in all_atoms:
            # search for the sub_match that contains the atom
            for sub_match in sub_matches:
                overlap = False
                atoms_to_block = _list_unnest_once(filtered_matches)
                if atom in sub_match:
                    for atom_of_match in sub_match:
                        if atom_of_match in atoms_to_block:
                            overlap = True
                            break
                    if not overlap:
                        filtered_matches.append(sub_match)
        return filtered_matches

    def __filter_intersub_overlap(
        self, sub_matches: list[list[int]], previous_matches: list[list[list[int]]]
    ) -> list[list[int]]:
        """filters out matches for one substructure key that have overlap with previous matches"""
        filtered_matches = []
        atoms_to_block = _list_unnest_twice(previous_matches)
        for sub_match in sub_matches:
            overlap = False  # reset
            for atom in sub_match:
                if atom in atoms_to_block:
                    overlap = True
                    break
            if not overlap:
                filtered_matches.append(sub_match)
        return filtered_matches

    def _filter_matches(self) -> list[list[list[int]]]:
        filtered_matches = []
        for i, substructure in enumerate(self.substructure_set):
            all_sub_matches = self.overlap_matches[i]  # all sub_matches
            if not self.intrasub_overlap:
                all_sub_matches = self.__filter_intrasub_overlap(all_sub_matches)
            if not self.intersub_overlap:
                all_sub_matches = self.__filter_intersub_overlap(
                    all_sub_matches, previous_matches=filtered_matches
                )
            filtered_matches.append(all_sub_matches)
        return filtered_matches

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
        total_atoms = nonh_atomcount(self.mol)

        logging.debug(f"matched, total: {matched_atoms}, {total_atoms}")

        return float(matched_atoms) / float(nonh_atomcount(self.mol))

    def _matches_to_coverage_per_substructure(self) -> list[float]:
        """gets non-h atom-based coverage of fingerprints over atoms"""
        mol_size = nonh_atomcount(self.mol)
        atoms_per_sub = self._matched_atoms_per_substructure()
        coverage_per_sub = [
            float(len(atoms)) / float(mol_size) for atoms in atoms_per_sub
        ]
        return coverage_per_sub

    def get_coverage(self) -> float:
        """gets non-h atom-based coverage of fingerprints over atoms"""
        return len(
            set(chain(*[match for match in self.matches for match in match]))
        ) / nonh_atomcount(self.mol)
        # coverage = self._matches_to_coverage()
        # self.coverage = coverage
        return coverage

    def get_coverage_per_substructure(self) -> list[float]:
        """gets non-h atom-based coverage of fingerprints over atoms"""
        coverage_per_sub = self._matches_to_coverage_per_substructure()
        self.coverage_per_sub = coverage_per_sub
        return coverage_per_sub

    def _matches_to_connections(self) -> list[float]:
        """under construction"""
        connections_bonds = []
        connections_atoms = []
        connections_subs = []
        connections = []

        # get all of the substructure matches-atoms
        atoms_per_sub = self._matched_atoms_per_substructure()
        # for each substructure combination, check if there is a bond between them
        for i, one_sub in enumerate(self.matches):
            for j, one_sub2 in enumerate(self.matches):
                for _, sub_matches in enumerate(one_sub):
                    for _, sub_matches2 in enumerate(one_sub2):
                        for atom in sub_matches:
                            for atom2 in sub_matches2:
                                if atom != atom2:
                                    bond = self.mol.GetBondBetweenAtoms(atom, atom2)
                                    if bond:
                                        connection_dict = {}
                                        # get atom indices of atom and atom2
                                        connection_dict["atoms"] = (atom, atom2)
                                        connection_dict["bond"] = bond.GetIdx()
                                        connection_dict["sub_index"] = (i, j)
                                        connections.append(connection_dict)

        sort = sorted(connections, key=lambda x: x["atoms"][0])

        return sort


# ------------------------------ multiple mol operations -------------------------------


class MolsCollection:
    def __init__(
        self, representations: list[Chem.Mol], repr_type: str = "smiles"
    ) -> iter:
        self.representations = representations
        self.repr_type = repr_type
        self.num_mols = self.get_num_mols()

    def __iter__(self) -> iter:
        """yields mols from inputlist"""
        appropriate_yielder = self._get_conversion_yielder()
        for mol in appropriate_yielder():
            yield mol

    def _get_conversion_yielder(self) -> Chem.Mol:
        """yields mols from inputlist"""
        types_functions = {
            "sdf": self._yield_sdf_mols,
            "smiles": self._yield_smiles_mols,
            "inchi": self._yield_inchi_mols,
        }
        return types_functions[self.repr_type]

    def _yield_smiles_mols(self) -> Chem.Mol:
        """yields mols from smileslist"""
        for mol in self.representations:
            yield smiles_to_mol(mol)

    def _yield_inchi_mols(self) -> Chem.Mol:
        """yields mols from inchilist"""
        for mol in self.representations:
            yield inchi_to_mol(mol)

    def _yield_sdf_mols(self) -> Chem.Mol:
        """yields mols from supplier (only first one of list)"""
        suppl = supplier(self.representations[0])
        for mol in suppl:
            yield mol

    def get_num_mols(self) -> int:
        if self.repr_type == "sdf":
            suppl = supplier(self.representations[0])
            return len(suppl)
        else:
            return len(self.representations)


# class Mols_sdf:
#     pass

#     def _yield_sdf_mols(sdflist: list) -> Chem.Mol:
#         """yields mols from supplier (only first one of list)"""
#         suppl = supplier(sdflist[0])
#         for mol in suppl:
#             yield mol

# def _get_yield_mol_function(input_type: str) -> Chem.Mol:
# """yields mols from inputlist"""
# types_functions = {
# "sdf": _yield_sdf_mols,
# "smiles": _yield_smiles_mols,
# "inchi": _yield_inchi_mols,
# }
# return types_functions[input_type]

# =========================== functions to import ==============================


def overlapped_fp(mol: Chem.Mol) -> list[int]:
    """gives the biosynfoni fingerprint for a mol, uses the default version as defined in def_biosynfoni.py
    defaults to full overlap allowance for Random Forest classification purposes(!)"""
    substructure_set = BiosynfoniVersion(defaultVersion).substructures
    mol_fp = Biosynfoni(
        mol=mol,
        substructure_set=substructure_set,
        intersub_overlap=True,
        intrasub_overlap=True,
    )
    return mol_fp.fingerprint


# ==============================  MAIN =======================================


def main():
    args = cli()
    silence_except_fp = args.nosave and not args.verbose  # default False
    if args.verbose:
        logging.getLogger("biosynfoni.concerto_fp").setLevel(logging.INFO)

    logging.info(10 * "=", "\nCONCERTO-FP\n", 10 * "=")
    logging.debug(args)

    # main functionality
    # get substructure set to speed up process (and not have to get it each time)
    substructure_set = BiosynfoniVersion(args.version).substructures

    # total_molecules = (
    # len(args.input) if args.repr != "sdf" else len(supplier(args.input[0]))
    # )
    # yield_mol_function = _get_yield_mol_function(args.repr)

    mol_collection = MolsCollection(args.input, args.repr)

    logging.info("looping over molecules...")
    biosynfonies, coverages = [], []
    for i, mol in tqdm(
        # enumerate(yield_mol_function(args.input)), total=total_molecules
        enumerate(mol_collection.__iter__()),
        total=mol_collection.num_mols,
        disable=silence_except_fp,
    ):
        if mol is None:
            with open("biosynfoni.err", "a") as ef:
                ef.write(f"{i}\t{args.input}\n")

        mol_fp = Biosynfoni(
            mol=mol,
            substructure_set=substructure_set,
            version_name=args.version,
            intersub_overlap=not args.intersub_block,
            intrasub_overlap=not args.intrasub_block,
        )
        if args.coverage:
            coverages.append(mol_fp.get_coverage())
        biosynfonies.append(mol_fp.fingerprint)

    if args.nosave:
        for row in biosynfonies:
            if not row:
                print("None")
            else:
                print(*row, sep=",")
        logging.info("coverages:")
        logging.info("\n".join(coverages))
        exit(0)
        return None

    if args.coverage:
        logging.info("writing coverages...")
        csv_writr(coverages, args.output.replace(".bsf", "_coverages.tsv"), sep="\t")

    logging.info(f"writing {len(biosynfonies)} biosynfonies to {args.output}...")

    biosynfonies_array = np.array(biosynfonies)
    np.savetxt(fname=args.output, X=biosynfonies_array, fmt="%i", delimiter=",")

    if args.repr != "sdf":
        with open(args.output.replace(".bsf", "_input.tsv"), "w") as f:
            f.write("\n".join([x for x in args.input]))
    # save_version(args.version, extra_text="extracted")
    # save version
    try:
        save_version(args.version)
    except:
        logging.warning("could not save version")

    logging.info("done")
    exit(0)
    return None


if __name__ == "__main__":
    main()
