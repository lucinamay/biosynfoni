#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
(unfinished: does not output yet)

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: CONCERTO-FP                  ||
language: python                    ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
student no: 1197142                 ||
____________________________________
 
||||||||||||  ()()()  |||||||||||||||

CONvert Compound E-Representation TO FingerPrint
description:    gets fp per molecule, also shows which fp version was used
"""

from rdkit import Chem
from sys import argv
import def_biosynfoni

#================================ prep =======================================
def get_supplier(sdf_file:str, supplier_only:bool = True)->list:
    suppl = Chem.SDMolSupplier(sdf_file)
    
    print('reading sdf, number of entries:', len(suppl))
    
    if supplier_only:
        return suppl
    else:
        mols = [x for x in suppl]
        return mols
    
def get_subsset(
            fp_version_name:str,
            subs_smarts:dict = def_biosynfoni.SUBSTRUCTURES,
            fp_versions:dict[list] = def_biosynfoni.FP_VERSIONS)->list:
    """ gives list of rdkit.Chem.Mols of substructures of choice 
    input:   fp_version_name (str) -- name of the version 
             subs_smarts (dict)
                         (key) substructure names (e.g. 'fp1')
                         (val) substructure RDK molfiles (f/SMARTS)
             fp_versions (dict)
                         (key) version name (e.g. fps_full_2)
                         (val) (list) substructure names (e.g. 'fp1')

    output: (list) rdkit.Chem.rdchem.Mol files for substructure keys
    """
    print("getting substructure set {}".format(fp_version_name))
    if not fp_version_name:
        raise 'No version name provided to select substructure set'
    
    substructures = []
    #getting the list of the chosen version's substructures
    chosen_version = fp_versions[fp_version_name]
    for substructure_name in chosen_version:  
        substructures.append(subs_smarts[substructure_name])
    return substructures

#=========================== per-mol functions ===============================

def detect_substructures(mol, 
                         substructure_set:list, 
                         matchnum:bool=False):
    """for mol, extracts fingerprint & atomnumbers of each match
    input: (Mol) mol-- RDKit Chem Mol object
           (list) fp_set -- list of RDKit Chem mol files of substructure keys
    output: 
        [0] (list[int])fingerprint 
        [1] (list[tuples[int]]) atomnumbers of each match per substructure
            
    """
    #output init.
    fingerprint = []
    all_unique_matches = [] #empty at start

    #temporary----------------------------------------------------------
    blockedout_atoms = []

    for substructure in substructure_set:
        num_approved_matches = 0 #number of non-overlapping matches
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
                num_approved_matches += 1       
                for unique_atom in approved_atoms:
                    # to keep track of which atoms have been used
                    blockedout_atoms.append(unique_atom)
                
        fingerprint.append(num_approved_matches)
        all_unique_matches.append(tuple(unique_sub_matches)) 

    if matchnum:
        return fingerprint, all_unique_matches
    else:
        return fingerprint

def get_biosynfoni(mol, version:str='', substructure_set:list=[],
                   matchnum:bool=False):
    """given a name of the fp version, uses get_subsset to get 
    the set of substructures, passing them todetect_substructures 
    to get the fingerprint. can also pass a substructure set directly"""
    print("getting fingerprint")
    if not substructure_set:
        substructure_set = get_subsset(version)
    biosynfoni = detect_substructures(mol, substructure_set,matchnum=matchnum)
    return biosynfoni

def get_atoms_per_sub(matches)->list[list[int]]:
    """per substructure, gives a list of the atom numbers that have been
    assigned to that substructure"""
    atoms_per_sub = []
    for sub in matches:
        sub_atomcount = list(sum(sub,()))
        atoms_per_sub.append(sub_atomcount)
    return atoms_per_sub

def count_listitems(nested_list:list)->int:
    count = 0
    for item in nested_list:
        if isinstance(item,list) or isinstance(item,tuple):
            sublist_count = count_listitems(list(item))
            count = count + sublist_count
        elif item:
            count+=1
    return count 

def get_coverage(mol, matches):
    matched_atoms = count_listitems(get_atoms_per_sub(matches))
    nonh_mol = Chem.rdmolops.RemoveHs(mol,
                                      implicitOnly=False,
                                      sanitize=True)
    mol_size =  nonh_mol.GetNumAtoms() #non-H atoms
    coverage = (float(matched_atoms)/float(mol_size))
    return coverage

def loop_over_supplier(
        supplier,
        fp_version:str='',
        substructure_set:list=[],
        subs_smarts:dict = def_biosynfoni.SUBSTRUCTURES,
        fp_versions:dict[list] = def_biosynfoni.FP_VERSIONS,
        coverage_info=False):
    """gets the fingerprint of the entire set"""
    
    #pass coverage_info bool to matchnum-returning option in get_biosynfoni
    matchnum=coverage_info
    print("looping over supplier, coverage info provided:", matchnum)
    
    #main functionality
    fingerprint_collection = []
    if matchnum:
        coverage_collection = []
        for mol in supplier:
            fingerprint, matches=get_biosynfoni(
                mol,
                version=fp_version,
                substructure_set=substructure_set,
                matchnum=matchnum)
            fingerprint_collection.append(fingerprint)
            coverage_collection.append(get_coverage(mol, matches))
    else:
        for mol in supplier:
            fingerprint=get_biosynfoni(
                mol,
                version=fp_version,
                substructure_set=substructure_set,
                matchnum=matchnum)
            fingerprint_collection.append(fingerprint)
    return fingerprint_collection


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
    
    supplier = get_supplier(supplier_loc)
    loop_over_supplier(supplier,
                       fp_version = fp_version,
                       coverage_info = coverage_info)

    mol_f_smiles = Chem.MolFromSmiles(
        "O=C(OC(C)CC=1C(OC)=C(O)C=2C(=O)C=C(OC)C3=C4C(OC)=CC(=O)C5=C(O)C(OC)=C(C(=C54)C1C23)CC(O)C)C=6C=CC=CC6")

if __name__ == "__main__":
    main()

