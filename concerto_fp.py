#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
(unfinished: does not output yet)

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: ConCERTO-fp                  ||
created: 2023-09-22                 ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
student no: 1197142                 ||
____________________________________
 
||||||||||||  ()()()  |||||||||||||||

'CONvert Compound E-Representation TO FingerPrint'
description: gets fp for given molecule, also shows which fp version was used
"""

from rdkit import Chem
from sys import argv
import def_biosynfoni
from inoutput import outfile_namer, fp_writr
import rdkfnx
import numpy as np
    
#=========================== per-mol functions ===============================
    
def detect_substructures(
        mol:Chem.Mol, 
        substructure_set:list
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
    #output init.
    all_unique_matches = [] #empty at start

    #temporary----------------------------------------------------------
    blockedout_atoms = []

    for substructure in substructure_set:
        sub_matches = mol.GetSubstructMatches(substructure) #all matches
        
        #start with all matches, then remove non-unique ones:
        unique_sub_matches = list(sub_matches)

        #removing overlapping finds from unique_matches
        for match in sub_matches:        #for each found match for curr subs
            approved_atoms = []       #unique matches for current subs
            
            #yields either all atoms for a match or none (when has an overlap)
            for atom in match:      
                if atom in blockedout_atoms:
                    #throw away match and all of its atoms, move onto next 
                    unique_sub_matches.remove(match)
                    approved_atoms = []
                    break
                elif atom not in blockedout_atoms:
                    approved_atoms.append(atom)

            #if every atom in match is free, approve the match:
            if len(approved_atoms) > 0:    #non-thrown away match
                # valid match, so add to actual fp:
                for unique_atom in approved_atoms:
                    # to keep track of which atoms have been used
                    blockedout_atoms.append(unique_atom)
                
        all_unique_matches.append(tuple(unique_sub_matches)) 
    return all_unique_matches

    
def matches_to_vector(matches:list[tuple[tuple[int]]]) -> list[int]:
    """"""
    vector = []
    for substructure_bit in matches:
        vector.append(len(substructure_bit))
    return vector


def get_biosynfoni(
        mol:Chem.Mol, 
        version:str = '', 
        substructure_set:list = [],
        return_matches:bool = False
        ) -> list[int]:
    """given a name of the fp version, uses get_subsset to get 
    the set of substructures, passing them todetect_substructures 
    to get the fingerprint. can also pass a substructure set directly"""
    
    print("getting fingerprint")
    assert (version or substructure_set),\
        'please give either the version name or the substructure set'
    if not substructure_set:
        substructure_set = def_biosynfoni.get_subsset(version)
    
    matches = detect_substructures(mol, substructure_set)
    biosynfoni = matches_to_vector(matches)
    if return_matches:
        return biosynfoni, matches
    else:
        return biosynfoni


def subs_assigned_atoms(matches:list[tuple[tuple[int]]]) -> list[list[int]]:
    """per substructure, gives a list of the atom numbers that have been
    assigned to that substructure"""
    atoms_per_sub = []
    for sub in matches:
        sub_atomcount = list(sum(sub,()))
        atoms_per_sub.append(sub_atomcount)
    return atoms_per_sub


def coverage_per_subs(matches:list[tuple[tuple[int]]]) -> list[float]:
    coverage = []
    for sub in subs_assigned_atoms(matches):
        coverage.append(len(sub))
    
    assert (len(coverage) == len(matches)),\
        "error in coverage finding: coverage and matches lengths do not match"
    return coverage


def count_listitems(nested_list:list) -> int:
    count = 0
    for item in nested_list:
        if isinstance(item,list) or isinstance(item,tuple):
            sublist_count = count_listitems(list(item)) #recursive
            count = count + sublist_count
        elif item:
            count+=1
    return count 


def get_coverage(
        mol:Chem.Mol,
        matches:list[tuple[tuple[int]]]
        ) -> float:
    """gets non-h atom-based coverage of fingerprints over atoms"""
    matched_atoms = count_listitems(subs_assigned_atoms(matches))
    nonh_mol = Chem.rdmolops.RemoveHs(mol,
                                      implicitOnly=False,
                                      sanitize=True)
    
    mol_size = nonh_mol.GetNumAtoms()   #non-H atoms
    coverage = (float(matched_atoms)/float(mol_size))
    return coverage

#==========================  main ============================================

def main():
    print(10*"=","\nCONCERTO-FP\n",10*"=")
    
    coverage_info = False #default
    supplier_loc = argv[1]
    fp_version= argv[2]
    if len(argv) == 4:
        coverage_info_bool= argv[3] #if want coverage: write True or T
        if coverage_info_bool in ['True','true','T','t',
                                  'coverage','y','yes','Y','Yes']:
            coverage_info=True
    
    supplier = rdkfnx.get_supplier(supplier_loc)
    biosynfonies = rdkfnx.loop_over_supplier(supplier,
                       fp_version = fp_version,
                       coverage_info = coverage_info)

    if coverage_info:
        fp_writr(biosynfonies[0], outfile_namer(supplier_loc, 'biosynfoni'))
        fp_writr(biosynfonies[1], outfile_namer(supplier_loc, 'atoms'))
    else:
        fp_writr(biosynfonies, outfile_namer(supplier_loc, 'biosynfoni'))
    
    fp_writr()

if __name__ == "__main__":
    main()

