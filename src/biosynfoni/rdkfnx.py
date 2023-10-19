# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:59:28 2023
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: RDKit functions              ||
created: 2023-09-13                 ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
____________________________________
 
||||||||||||  ()()()  |||||||||||||||

contains general reusable functionalities depending on RDKit packages,
mainly supplier-handling
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger  # for muting warnings

from biosynfoni.def_biosynfoni import (
    SUBSTRUCTURES,
    FP_VERSIONS,
    get_smarts,
    DEFAULT_BIOSYNFONI_VERSION,
)
from biosynfoni.inoutput import outfile_namer, open_json
import biosynfoni.leaf_subs_nod12 as leaffile


def sdf_writr(mols: list, outfile: str) -> None:
    """writes sdf of mols,
    * needs RDKit Chem
    """
    writer = Chem.SDWriter(outfile)
    for mol in mols:
        writer.write(mol)
    return None


def get_supplier(sdf_file: str, supplier_only: bool = True) -> list:
    suppl = Chem.SDMolSupplier(sdf_file)

    print("reading sdf, number of entries:", len(suppl))

    if supplier_only:
        return suppl
    else:
        mols = [x for x in suppl]
        return mols


def get_leaf_substructures(dict_list: list[dict[str, str]]) -> list[Chem.Mol]:
    leaf_substructures = []
    names = []
    for subs_dict in dict_list:
        mol = Chem.MolFromSmarts(subs_dict["smarts"])
        if isinstance(mol, Chem.Mol):
            leaf_substructures.append(mol)
            names.append(subs_dict["name"])
        else:
            print(
                f"substructure {subs_dict['name']} could not be converted to mol: skipped"
            )

    print("extracted:")
    print("\n".join(names))
    # with open(outfile_namer('json_saved.txt'),'w') as f:
    #    f.write('\n'.join(names))
    return leaf_substructures


def get_subsset(
    fp_version_name: str,
    subs_smarts: dict = SUBSTRUCTURES,
    fp_versions: dict[list] = FP_VERSIONS,
) -> list[Chem.Mol]:
    """gives list of rdkit.Chem.Mols of substructures of choice
    input:   fp_version_name (str) -- name of the version
             subs_smarts (dict)
                         (key) substructure names (e.g. 'fp1')
                         (val) substructure RDK molfiles (f/SMARTS)
             fp_versions (dict)
                         (key) version name (e.g. fps_full_2)
                         (val) (list) substructure names (e.g. 'fp1')

    output: (list) rdkit.Chem.rdchem.Mol files for substructure keys
    """
    # print(f"getting substructures from set {fp_version_name}")

    # ===========================================================================
    # test for leaf:
    if fp_version_name == "leaf":
        return get_leaf_substructures(leaffile.leaf)

    if not fp_version_name:
        raise "No version name provided to select substructure set"
    # ===========================================================================

    substructures = []
    successful_subs_names = []
    # getting the list of the chosen version's substructures
    chosen_version = fp_versions[fp_version_name]
    for substructure_name in chosen_version:
        if isinstance(subs_smarts[substructure_name], Chem.Mol):
            substructures.append(subs_smarts[substructure_name])
        elif isinstance(subs_smarts[substructure_name], str):
            dirtymol = Chem.MolFromSmarts(subs_smarts[substructure_name])
            if isinstance(dirtymol, Chem.Mol):
                substructures.append(dirtymol)
                successful_subs_names.append(substructure_name)
            else:
                print(
                    f"substructure {substructure_name} could not be converted to mol: skipped"
                )
    # print(f"added {','.join([x for x in successful_subs_names])} to substructure set")
    return substructures


def save_version(
    fp_version: str, window_size=(1000, 1000), extra_text: str = ""
) -> None:
    outfilename = outfile_namer("version", f"{fp_version}_{extra_text}")

    # svg_text = fm.drawfp(fp_version, window_size=window_size)
    # with open(f'{outfilename}.svg', 'w') as f:
    #    f.write(svg_text)

    smarts = get_smarts(fp_version)
    with open(f"{outfilename}.smarts", "w") as f:
        for smart in smarts:
            f.write(f"{smart[0]}\t{smart[1]}\n")
    return None
