#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:59:28 2023

@author: lucina-may
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger  #for muting warnings
import def_biosynfoni

def sdf_writr(mols:list, outfile:str) -> None:
    """writes sdf of mols, 
    * needs RDKit Chem
    """
    writer = Chem.SDWriter(outfile)
    for mol in mols:
        writer.write(mol)
    return None

def get_supplier(sdf_file:str, supplier_only:bool = True)->list:
    suppl = Chem.SDMolSupplier(sdf_file)
    
    print('reading sdf, number of entries:', len(suppl))
    
    if supplier_only:
        return suppl
    else:
        mols = [x for x in suppl]
        return mols
    
def loop_over_supplier(
        supplier,
        fp_version:str='',
        substructure_set:list=[],
        subs_smarts:dict = def_biosynfoni.SUBSTRUCTURES,
        fp_versions:dict[list] = def_biosynfoni.FP_VERSIONS,
        coverage_info=False):
    """gets the fingerprint of the entire set"""
    
    #pass coverage_info bool to matches-returning option in get_biosynfoni
    matches=coverage_info
    print("looping over supplier, coverage info provided:", matches)
    
    #main functionality
    fingerprint_collection = []
    if matches:
        coverage_collection = []
        for mol in supplier:
            fingerprint, matches=get_biosynfoni(
                mol,
                version=fp_version,
                substructure_set=substructure_set,
                matches=matches)
            fingerprint_collection.append(fingerprint)
            coverage_collection.append(get_coverage(mol, matches))
    else:
        for mol in supplier:
            fingerprint=get_biosynfoni(
                mol,
                version=fp_version,
                substructure_set=substructure_set,
                matches=matches)
            fingerprint_collection.append(fingerprint)
    return fingerprint_collection
