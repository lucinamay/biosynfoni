#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: definition biosynfoni        ||
created: 2023-09                    ||
author: Lucina-May Nollen           ||
institute: WUR Bioinformatics       ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

Biosynfoni Definition. Can be converted into an sdf file per 'version'

"""
from rdkit import Chem
from enum import Enum

DEFAULT_BIOSYNFONI_VERSION = 'regular_1008'  # not used here, just for reference

# ============================
SUBSTRUCTURES = {
    # dewick's building blocks (d_), conserved aromaticity, free rest ------
    # shikimate
    'd_indoleC2N_12': 'c1cccc2c1c(~[#6]~[#6]~[#7])cn2',
    'd_phenylC2N_9': 'c1ccccc1[#6][#6][#7]',
    # acetate
    'd_c5n_6': '[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1',
    'd_c4n_5': '[#6]~1~[#6]~[#6]~[#6]~[#7]~1',
    # shikimate
    'd_phenylC3_9_1007': 'c1ccccc1~[#6]~[#6]~[#6]',
    'd_phenylC2_8_1007': 'c1ccccc1~[#6]~[#6]',
    'd_phenylC1_7_1007': 'c1ccccc1~[#6]',
    'd_phenylC3_9': '[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]',
    'd_phenylC2_8': '[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]',
    'd_phenylC1_7': '[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]',
    # mevalonate/MEP
    'd_isoprene_5': '[#6]~[#6](~[#6])~[#6]~[#6]',
    # acetate, aromaticity ok
    'd_ethyl_2': '[#6]~[#6]',
    'd_methyl_1': '[CH3]',  # only end methyls
    # sugar-related --------------------------------------------------------
    's_pyranose_C5O4': '[#6]~1~[#8]~[#6]~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~1',
    's_furanose_C4O3': '[#6]~1~[#8]~[#6]~[#6](~[#8])~[#6](~[#8])~1',
    's_openpyr_C6O6': ('[#6](~[#8])~[#6](~[#8])~[#6]'
                       '(~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])'),
    's_openfur_C5O5': ('[#6](~[#8])~[#6](~[#8])~[#6]'
                       '(~[#8])~[#6](~[#8])~[#6](~[#8])'),
    # additional from dewick -----------------------------------------------
    # acetate
    'd2_acetyl_C2O1': '[#6][#6]~[#8]',  # CC=O
    'd2_methylmalonyl_C3': '[#6]~[#6][CH3]',  # check coverage
    # halogens -------------------------------------------------------------
    'hal_f': '[#9]',
    'hal_cl': '[#17]',
    'hal_br': '[#35]',
    'hal_i': '[#53]',
    # phosphates, sulfonates -----------------------------------------------
    'phosphate_2': 'P~O',  # if at end
    'phosphate_5': 'O~P(~O)(~O)~O',
    'sulfonate_2': 'S~O',
    'sulfonate_5': 'O~S(~O)(~O)~O',
    # additional
    'n_nitrate_1': '[NX3;H2,H1]',
    'o_epoxy_1': '[!H]1-O-[!H]1',  # take only oxygen
    'o_ether_1': '[!H]-O-[!H]',  # check detection #+epoxy?
    'o_hydroxyl_1': '[OH]',  # check detection
    # Coenzymes, catalytic units etc ---------------------------------------
    'co_coa_pubchem': (
        'SCCNC(~O)CCNC(~O)C(C(C)(C)COP(=O)(~O)OP(=O)'
        '(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O'
    ),
    'co_coa': (
        'SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)'
        '(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O'
    ),
    'co_nadh': (
        'C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)'
        '(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O'
    ),
    'co_nadph': (
        'C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)'
        'OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O'
    )
    # maybe do a loop where we do substructure searches intra-fingerprint
    # any matches get a 'distance' of one
}


class Biosynfoni(Enum):
    co_coa = 0
    co_nadh = 1
    co_nadph = 2
    s_pyranose_C5O4 = 3
    s_furanose_C4O3 = 4
    s_openpyr_C6O6 = 5
    s_openfur_C5O5 = 6
    d_indoleC2N_12 = 7
    d_phenylC2N_9 = 8
    d_c5n_6 = 9
    d_c4n_5 = 10
    d_phenylC3_9 = 11
    d_phenylC2_8 = 12
    d_phenylC1_7 = 13
    d_isoprene_5 = 14
    d2_acetyl_C2O1 = 15
    d2_methylmalonyl_C3 = 16
    d_ethyl_2 = 17
    d_methyl_1 = 18
    phosphate_2 = 19
    sulfonate_2 = 20
    hal_f = 21
    hal_cl = 22
    hal_br = 23
    hal_i = 24
    n_nitrate_1 = 25
    o_epoxy_1 = 26
    o_ether_1 = 27
    o_hydroxyl_1 = 28

    def as_fp(self):
        return [Chem.Mol(SUBSTRUCTURES[x]) for x in self.keys]


SUBS_PATHWAYS = {
    'd_indoleC2N_12': ['shikimate'],
    'd_phenylC2N_9': ['shikimate'],
    # acetate
    'd_c5n_6': ['acetate'],
    'd_c4n_5': ['acetate'],
    # shikimate
    'd_phenylC3_9': ['shikimate'],
    'd_phenylC2_8': ['shikimate'],
    'd_phenylC1_7': ['shikimate'],
    # mevalonate/MEP
    'd_isoprene_5': ['mevalonate', 'methylerythritol phosphate'],
    # acetate, aromaticity ok
    'd_ethyl_2': ['acetate'],
    'd_methyl_1': ['acetate'],  # only end methyls
    # sugar-related --------------------------------------------------------
    's_pyranose_C5O4': ['sugar'],
    's_furanose_C4O3': ['sugar'],
    's_openpyr_C6O6': ['sugar'],
    's_openfur_C5O5': ['sugar'],
    # additional from dewick -----------------------------------------------
    # acetate
    'd2_acetyl_C2O1': ['acetate'],
    'd2_methylmalonyl_C3': ['acetate']
}

FP_VERSIONS = {
    'fps_full_0814': ['fp8', 'fp12', 'fp9', 'fp10', 'fp11', 'fp13',
                      'fp1', 'fp2', 'fp3', 'fp4', 'fp5', 'fp6', 'fp7',
                      'hal1', 'hal2', 'hal3', 'hal4',
                      'fp98', 'fp99'],

    'fps_full_2_0814': ['fp8_2', 'fp12_2', 'fp9_2', 'fp10_2', 'fp11_2', 'fp13_2',
                        'fp1_2', 'fp2_2', 'fp3_2', 'fp4_2', 'fp5_2', 'fp6_2', 'fp7_2',
                        'hal1_2', 'hal2_2', 'hal3_2', 'hal4_2',
                        'fp98_2', 'fp99_2'],

    'fps_full_3_0814': ['fp8_3', 'fp12_3', 'fp9_3', 'fp10_3', 'fp11_3', 'fp13_3',
                        'fp1_3', 'fp2_3', 'fp3_3', 'fp4_3', 'fp5_3', 'fp6_3', 'fp7_3',
                        'hal1_2', 'hal2_2', 'hal3_2', 'hal4_2',
                        'fp98_3', 'fp99_3_old'],
    'regular_1007': [
        'co_coa',
        'co_nadh',
        'co_nadph',
        's_pyranose_C5O4',
        's_furanose_C4O3',
        's_openpyr_C6O6',
        's_openfur_C5O5',
        'd_indoleC2N_12',
        'd_phenylC2N_9',
        'd_c5n_6',
        'd_c4n_5',
        'd_phenylC3_9_1007',
        'd_phenylC2_8_1007',
        'd_phenylC1_7_1007',
        'd_isoprene_5',
        'd2_acetyl_C2O1',
        'd2_methylmalonyl_C3',
        'd_ethyl_2',
        'd_methyl_1',
        'phosphate_2',
        'sulfonate_2',
        'hal_f',
        'hal_cl',
        'hal_br',
        'hal_i',
        'n_nitrate_1',
        'o_epoxy_1',
        'o_ether_1',
        'o_hydroxyl_1',
    ],
    'regular_1008': [  # new version with d_phenylC1_7, after isoprene
        'co_coa',
        'co_nadh',
        'co_nadph',
        's_pyranose_C5O4',
        's_furanose_C4O3',
        's_openpyr_C6O6',
        's_openfur_C5O5',
        'd_indoleC2N_12',
        'd_phenylC2N_9',
        'd_c5n_6',
        'd_c4n_5',
        'd_phenylC3_9',
        'd_phenylC2_8',
        'd_phenylC1_7',
        'd_isoprene_5',
        # 'd2_acetyl_C2O1',
        # 'd2_methylmalonyl_C3',
        'd_ethyl_2',
        'd_methyl_1',
        'phosphate_2',
        'sulfonate_2',
        'hal_f',
        'hal_cl',
        'hal_br',
        'hal_i',
        'n_nitrate_1',
        'o_epoxy_1',
        'o_ether_1',
        'o_hydroxyl_1',
    ]
}

# --------------------------- substructure sets ---------------------------
SUBSTRUCTURES_old = {
    'fp1': Chem.MolFromSmarts(
        'c1cccc2c1c(CCN)cn2'),
    'fp1_2ar': Chem.MolFromSmarts(
        'c1cccc2c1c([#6][#6][#7])cn2'),
    'fp1_2': Chem.MolFromSmarts(
        '[#6]1[#6][#6][#6][#6]2[#6]1[#6]([#6][#6][#7])[#6][#7]2'),
    'fp1_3': Chem.MolFromSmarts(
        '[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6](~[#6]~[#6]~[#7])~[#6]~[#7]~2'),\
    # fp2  -- Phe.C2N --  shikimate
    'fp2': Chem.MolFromSmarts('c1ccccc1CCN'),\
    'fp2_2ar': Chem.MolFromSmarts('c1ccccc1[#6][#6][#7]'),\
    'fp2_2': Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1[#6][#6][#7]'),\
    'fp2_3': Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#7]'),\
    # fp3 -- cyclic C5N -- acetylCoA  (Krebs)
    'fp3': Chem.MolFromSmarts('C1CCCCN1'),\
    'fp3_2': Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#7]1'),\
    'fp3_3': Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1'),\
    # fp4 -- cyclic C4N -- acetylCoA (Krebs)
    'fp4': Chem.MolFromSmarts('C1CCCN1'),\
    'fp4_2': Chem.MolFromSmarts('[#6]1[#6][#6][#6][#7]1'),\
    'fp4_3': Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#7]~1'),\
    # fp49-- c6C3 -- shikimate
    'fp49': Chem.MolFromSmarts('c1ccccc1CCC'),\
    'fp49_2ar': Chem.MolFromSmarts('c1ccccc1[#6][#6][#6]'),\
    'fp49_2': Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1[#6][#6][#6]'),\
    'fp49_3': Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]'),\
    # fp5 -- c6C2 -- shikimate
    'fp5': Chem.MolFromSmarts('c1ccccc1CC'),\
    'fp5_2ar': Chem.MolFromSmarts('c1ccccc1[#6][#6]'),\
    'fp5_2': Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1[#6][#6]'),\
    'fp5_3': Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]'),\
    # fp6 -- c6C1 -- shikimate
    'fp6': Chem.MolFromSmarts('c1ccccc1C'),\
    'fp6_2ar': Chem.MolFromSmarts('c1ccccc1[#6]'),\
    'fp6_2': Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1[#6]'),\
    'fp6_3': Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]'),\
    # fp7 -- isoprene -- mevalonic/methylerythriol
    'fp7': Chem.MolFromSmarts('CC(C)CC'),\
    'fp7_2': Chem.MolFromSmarts('[#6][#6]([#6])[#6][#6]'),\
    'fp7_3': Chem.MolFromSmarts('[#6]~[#6](~[#6])~[#6]~[#6]'),\
    'fp7_2_old': Chem.MolFromSmarts('[#6]~[#6]([#6])~[#6][#6]'),\
    # in the stats, biochemically correcter(?why)
    # fp98 -- C2 -- acetylCoA  *'acetyls also end up in aromatic systems', Dewick
    'fp98': Chem.MolFromSmarts('CC'),\
    'fp98_2': Chem.MolFromSmarts('[#6][#6]'),\
    'fp98_3': Chem.MolFromSmarts('[#6]~[#6]'),\
    # fp99 -- C -- acetylCoA * methyl only (matches as 1 atom)
    'fp99': Chem.MolFromSmarts('[CH3]'),\
    'fp99_2': Chem.MolFromSmarts('C'),\
    ''
    'fp99_3_old': Chem.MolFromSmarts('[#6]'),\
    # === personal additions
    # phosphate
    'fp8': Chem.MolFromSmarts('O~P(~O)(~O)~O'),\
    'fp8_2': Chem.MolFromSmarts('O~P(~O)~O'),\
    'fp8_3': Chem.MolFromSmarts('P~O'),\
    # sulfate group
    'fp12': Chem.MolFromSmarts('O~S(~O)(~O)~O'),\
    'fp12_2': Chem.MolFromSmarts('O~S(~O)~O'),\
    'fp12_3': Chem.MolFromSmarts('S~O'),\
    # --- sugars
    # pyranose -- sugar
    'fp9': Chem.MolFromSmarts('C1OCC(O)C(O)C(O)1'),\
    'fp9_2': Chem.MolFromSmarts(
        '[#6]1[#8][#6][#6]([#8])[#6]([#8])[#6]([#8])1'),\
    'fp9_3': Chem.MolFromSmarts(
        '[#6]~1~[#8]~[#6]~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~1'),\
    # furanose -- sugar
    'fp10': Chem.MolFromSmarts('C1OCC(O)C(O)1'),\
    'fp10_2': Chem.MolFromSmarts(
        '[#6]1[#8][#6][#6]([#8])[#6]([#8])1'),\
    'fp10_3': Chem.MolFromSmarts(
        '[#6]~1~[#8]~[#6]~[#6](~[#8])~[#6](~[#8])~1'),
    # open pyranose -- sugar
    'fp11': Chem.MolFromSmarts('C(O)C(O)C(O)C(O)C(O)C(O)'),
    'fp11_2': Chem.MolFromSmarts(
        '[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])'),
    'fp11_3': Chem.MolFromSmarts(
        '[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])'),
    # open furanose -- sugar
    'fp13': Chem.MolFromSmarts('C(O)C(O)C(O)C(O)C(O)'),\
    'fp13_2': Chem.MolFromSmarts(
        '[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])'),\
    'fp13_3': Chem.MolFromSmarts(
        '[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])'),\
    # halogens
    'hal1': Chem.MolFromSmarts('F'),\
    'hal1_2': Chem.MolFromSmarts('[#9]'),\
    'hal2': Chem.MolFromSmarts('Cl'),\
    'hal2_2': Chem.MolFromSmarts('[#17]'),\
    'hal3': Chem.MolFromSmarts('Br'),\
    'hal3_2': Chem.MolFromSmarts('[#35]'),\
    'hal4': Chem.MolFromSmarts('I'),\
    'hal4_2': Chem.MolFromSmarts('[#53]'),\

    # additionals:
    'acetylO': Chem.MolFromSmarts('CC=O'),
    # double and single bond only ?
    'acetylO_2': Chem.MolFromSmarts('[#6][#6]O'),
}


# --------------------------- pathway annotation --------------------------------

SUBS_PATHWAYS_old = {
    'fp1': ['shikimate'],
    'fp2': ['shikimate'],
    'fp3': ['acetate'],
    'fp4': ['acetate'],
    'fp49': ['shikimate'],
    'fp5': ['shikimate'],
    'fp6': ['shikimate'],
    'fp7': ['mevalonate', 'methylerythritol phosphate'],
    'fp98': ['acetate'],
    'fp99': ['acetate'],
    'fp11': ['sugar'],
    'fp13': ['sugar'],
    'fp9': ['sugar'],
    'fp10': ['sugar']
}

# --------------------------- get substructure set ----------------------------


def get_smarts(
        fp_version_name: str,
        subs_smarts: dict = SUBSTRUCTURES,
        fp_versions: dict[list] = FP_VERSIONS) -> list[list[str, str]]:
    """ gives list of smarts of substructures of choice
    input:   fp_version_name (str) -- name of the version
             subs_smarts (dict)
                         (key) substructure names (e.g. 'fp1')
                         (val) substructure RDK molfiles (f/SMARTS)
             fp_versions (dict)
                         (key) version name (e.g. fps_full_2)
                         (val) (list) substructure names (e.g. 'fp1')
    """
    chosen_sub_names = fp_versions[fp_version_name]
    return [[x, subs_smarts[x]] for x in chosen_sub_names]


def get_subsset(
        fp_version_name: str,
        subs_smarts: dict = SUBSTRUCTURES,
        fp_versions: dict[list] = FP_VERSIONS) -> list:
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
    # print(f"getting substructures from set {fp_version_name}")
    if not fp_version_name:
        raise 'No version name provided to select substructure set'

    substructures = []
    # getting the list of the chosen version's substructures
    chosen_version = fp_versions[fp_version_name]
    for substructure_name in chosen_version:
        if isinstance(subs_smarts[substructure_name], Chem.Mol):
            substructures.append(subs_smarts[substructure_name])
        elif isinstance(subs_smarts[substructure_name], str):
            dirtymol = Chem.MolFromSmarts(subs_smarts[substructure_name])
            substructures.append(dirtymol)
    return substructures

# get_smarts(DEFAULT_BIOSYNFONI_VERSION, SUBSTRUCTURES, FP_VERSIONS)
# get_subsset(DEFAULT_BIOSYNFONI_VERSION, SUBSTRUCTURES, FP_VERSIONS)

# get_subsset('regular_1007')
