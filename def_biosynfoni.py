#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Biosynfoni Definition. Can be converted into an sdf file per 'version'

@author: lucina-may
"""
from rdkit import Chem

SUBSTRUCTURES = { 
'fp1' : Chem.MolFromSmarts(
    'c1cccc2c1c(CCN)cn2'),\
'fp1_2ar' : Chem.MolFromSmarts(
    'c1cccc2c1c([#6][#6][#7])cn2'),\
'fp1_2' : Chem.MolFromSmarts(
    '[#6]1[#6][#6][#6][#6]2[#6]1[#6]([#6][#6][#7])[#6][#7]2'),\
'fp1_3' : Chem.MolFromSmarts(
    '[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6](~[#6]~[#6]~[#7])~[#6]~[#7]~2'),\
#fp2  -- Phe.C2N --  shikimate
'fp2' : Chem.MolFromSmarts('c1ccccc1CCN'),\
'fp2_2ar' : Chem.MolFromSmarts('c1ccccc1[#6][#6][#7]'),\
'fp2_2' : Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1[#6][#6][#7]'),\
'fp2_3' : Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#7]'),\
#fp3 -- cyclic C5N -- acetylCoA  (Krebs)
'fp3' : Chem.MolFromSmarts('C1CCCCN1'),\
'fp3_2' : Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#7]1'),\
'fp3_3' : Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1'),\
#fp4 -- cyclic C4N -- acetylCoA (Krebs)
'fp4' : Chem.MolFromSmarts('C1CCCN1'),\
'fp4_2' : Chem.MolFromSmarts('[#6]1[#6][#6][#6][#7]1'),\
'fp4_3' : Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#7]~1'),\
#fp49-- c6C3 -- shikimate
'fp49' : Chem.MolFromSmarts('c1ccccc1CCC'),\
'fp49_2ar' : Chem.MolFromSmarts('c1ccccc1[#6][#6][#6]'),\
'fp49_2' : Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1[#6][#6][#6]'),\
'fp49_3' : Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]'),\
#fp5 -- c6C2 -- shikimate
'fp5' : Chem.MolFromSmarts('c1ccccc1CC'),\
'fp5_2ar' : Chem.MolFromSmarts('c1ccccc1[#6][#6]'),\
'fp5_2' : Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1[#6][#6]'),\
'fp5_3' : Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]'),\
#fp6 -- c6C1 -- shikimate
'fp6' : Chem.MolFromSmarts('c1ccccc1C'),\
'fp6_2ar' : Chem.MolFromSmarts('c1ccccc1[#6]'),\
'fp6_2' : Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1[#6]'),\
'fp6_3' : Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]'),\
#fp7 -- isoprene -- mevalonic/methylerythriol
'fp7' : Chem.MolFromSmarts('CC(C)CC'),\
'fp7_2' : Chem.MolFromSmarts('[#6][#6]([#6])[#6][#6]'),\
'fp7_3' : Chem.MolFromSmarts('[#6]~[#6](~[#6])~[#6]~[#6]'),\
'fp7_2_old' : Chem.MolFromSmarts('[#6]~[#6]([#6])~[#6][#6]'),\
                  #in the stats, biochemically correcter
#fp98 -- C2 -- acetylCoA  *'acetyls also end up in aromatic systems', Dewick
'fp98' : Chem.MolFromSmarts('CC'),\
'fp98_2' : Chem.MolFromSmarts('[#6][#6]'),\
'fp98_3' : Chem.MolFromSmarts('[#6]~[#6]'),\
#fp99 -- C -- acetylCoA * methyl only
'fp99' : Chem.MolFromSmarts('[CH3]'),\
'fp99_2' : Chem.MolFromSmarts('C'),\
'fp99_3' : Chem.MolFromSmarts('[#6]'),\
# === personal additions
#phosphate
'fp8' : Chem.MolFromSmarts('O~P(~O)(~O)~O'),\
'fp8_2' : Chem.MolFromSmarts('O~P(~O)~O'),\
'fp8_3' : Chem.MolFromSmarts('P~O'),\
#sulfate group
'fp12' : Chem.MolFromSmarts('O~S(~O)(~O)~O'),\
'fp12_2' : Chem.MolFromSmarts('O~S(~O)~O'),\
'fp12_3' : Chem.MolFromSmarts('S~O'),\
# --- sugars
#pyranose -- sugar
'fp9' : Chem.MolFromSmarts('C1OCC(O)C(O)C(O)1'),\
'fp9_2' : Chem.MolFromSmarts('[#6]1[#8][#6][#6]([#8])[#6]([#8])[#6]([#8])1'),\
'fp9_3' : Chem.MolFromSmarts('[#6]~1~[#8]~[#6]~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~1'),\
#furanose -- sugar
'fp10' : Chem.MolFromSmarts('C1OCC(O)C(O)1'),\
'fp10_2' : Chem.MolFromSmarts('[#6]1[#8][#6][#6]([#8])[#6]([#8])1'),\
'fp10_3' : Chem.MolFromSmarts('[#6]~1~[#8]~[#6]~[#6](~[#8])~[#6](~[#8])~1'),\
#open pyranose -- sugar
'fp11' : Chem.MolFromSmarts('C(O)C(O)C(O)C(O)C(O)C(O)'),\
'fp11_2' : Chem.MolFromSmarts('[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])'),\
'fp11_3' : Chem.MolFromSmarts('[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])'),\
#open furanose -- sugar
'fp13' : Chem.MolFromSmarts('C(O)C(O)C(O)C(O)C(O)'),\
'fp13_2' : Chem.MolFromSmarts('[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])'),\
'fp13_3' : Chem.MolFromSmarts('[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])'),\
#halogens
'hal1' : Chem.MolFromSmarts('F'),\
'hal1_2' : Chem.MolFromSmarts('[#9]'),\
'hal2' : Chem.MolFromSmarts('Cl'),\
'hal2_2' : Chem.MolFromSmarts('[#17]'),\
'hal3' : Chem.MolFromSmarts('Br'),\
'hal3_2' : Chem.MolFromSmarts('[#35]'),\
'hal4' : Chem.MolFromSmarts('I'),\
'hal4_2' : Chem.MolFromSmarts('[#53]'),\
}

#--------------------------- substructure sets --------------------------------
FP_VERSIONS = {
    'fps_full': ['fp8', 'fp12', 'fp9', 'fp10', 'fp11', 'fp13',\
            'fp1', 'fp2', 'fp3', 'fp4', 'fp5', 'fp6', 'fp7',\
            'hal1', 'hal2', 'hal3', 'hal4',\
            'fp98', 'fp99'],

    'fps_full_2':['fp8_2', 'fp12_2', 'fp9_2', 'fp10_2', 'fp11_2','fp13_2',\
              'fp1_2', 'fp2_2', 'fp3_2', 'fp4_2', 'fp5_2', 'fp6_2', 'fp7_2',\
              'hal1_2', 'hal2_2', 'hal3_2', 'hal4_2',\
              'fp98_2', 'fp99_2'],
    
    'fps_full_3' : ['fp8_3', 'fp12_3', 'fp9_3', 'fp10_3', 'fp11_3', 'fp13_3',\
              'fp1_3', 'fp2_3', 'fp3_3', 'fp4_3', 'fp5_3', 'fp6_3', 'fp7_3',\
              'hal1_2', 'hal2_2', 'hal3_2', 'hal4_2',\
              'fp98_3', 'fp99_3']
}
    
#--------------------------- pathway annotation --------------------------------

SUBS_PATHWAYS = {
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
    'fp10':['sugar']
    }

#probably change this to a file-system

