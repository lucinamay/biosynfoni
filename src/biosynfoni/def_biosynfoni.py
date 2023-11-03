"""
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: biosynfoni arrangement       ||
created: 2023-09                    ||
author: Lucina-May Nollen           ||
institute: WUR Bioinformatics       ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

Biosynfoni Definition. Can be converted into an sdf file per 'version'

"""

from enum import Enum

DEFAULT_BIOSYNFONI_VERSION = "full_1103"  # not used here, just for reference

# ============================
SUBSTRUCTURES = {
    # dewick's building blocks (d_), conserved aromaticity, free rest ------
    # shikimate
    "d_indoleC2N_12": "c1cccc2c1c(~[#6]~[#6]~[#7])cn2",
    "d_phenylC2N_9": "c1ccccc1[#6][#6][#7]",
    # acetate
    "d_c5n_6": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1",
    "d_c4n_5": "[#6]~1~[#6]~[#6]~[#6]~[#7]~1",
    # shikimate
    "d_phenylC3_9_1007": "c1ccccc1~[#6]~[#6]~[#6]",
    "d_phenylC2_8_1007": "c1ccccc1~[#6]~[#6]",
    "d_phenylC1_7_1007": "c1ccccc1~[#6]",
    "d_phenylC3_9": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6;!$([r6])]~[#6;!$([r6])]",  # last two should not be next to ring
    "d_phenylC2_8": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6;!$([r6])]",  # last one should not be next to ring
    "d_phenylC1_7": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]",  # this one is allowed to be ring, to catch all
    "d_phenylC3_9_strict": "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]~[#6;!$([r6])]",  # not allowed to be in fused rings
    "d_phenylC2_8_strict": "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]",
    "d_phenylC1_7_strict": "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]",  # no fused rings at all
    # mevalonate/MEP
    "d_isoprene_5": "[#6]~[#6](~[#6])~[#6]~[#6]",
    # acetate, aromaticity ok
    "d_ethyl_2": "[#6]~[#6]",
    "d_methyl_1": "[C;D1;h3]",  # only end methyls
    # sugar-related --------------------------------------------------------
    "s_pyranose_C5O4": "C~1~[#8]~C~C(~[#8])~C(~[#8])~C(~[#8])~1",
    "s_furanose_C4O3": "C~1~[#8]~C~C(~[#8])~C(~[#8])~1",
    "s_openpyr_C6O6": "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
    "s_openfur_C5O5": "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
    # additional from dewick -----------------------------------------------
    # acetate
    "d2_acetyl_C2O1": "[#6]~[#6]~[#8]",  # CC=O
    "d2_methylmalonyl_C3": "[#6]~[#6][C;D1;h3]",  # check coverage
    # halogens -------------------------------------------------------------
    "hal_f": "[#9]",
    "hal_cl": "[#17]",
    "hal_br": "[#35]",
    "hal_i": "[#53]",
    # phosphates, sulfonates -----------------------------------------------
    "phosphate_2": "P~O",  # if at end
    "phosphate_5": "O~P(~O)(~O)~O",
    "sulfonate_2": "S~O",
    "sulfonate_5": "O~S(~O)(~O)~O",
    # additional
    "n_nitrate_1": "[N;D1]",
    "o_epoxy_1": "[O;x2;r3]",  # counts only the oxygen (2ringbond,size3)
    "o_ether_1": "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S)]",  # not ester,twoconn,noringbond,noH
    "o_hydroxyl_1": "[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(P);!$(S)]",  # OH, O-,  only attached to C/N, no acid hydroxyls, no phosphate/sulfonate,
    # Coenzymes, catalytic units etc ---------------------------------------
    "co_coa_pubchem": (
        "SCCNC(~O)CCNC(~O)C(C(C)(C)COP(=O)(~O)OP(=O)"
        "(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O"
    ),
    "co_coa": (
        "SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)"
        "(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O"
    ),
    "co_nadh": (
        "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)"
        "(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"
    ),
    "co_nadph": (
        "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)"
        "OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O"
    ),
    # amino acids:
    "allstnd_aminos": (
        "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),"
        "$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),"
        "$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N])]"  # from daylight.com A.A. Template for 20 standard a.a.s
    ),
    "nonstnd_aminos": (
        "[$([NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]);"
        "!$([$([$([NX3H,NX4H2+]),"
        "$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),"
        "$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),"
        "$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),"
        "$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),"
        "$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),"
        "$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),"
        "$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),"
        "$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),"
        "$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),"
        "$([CH2X4][CHX4]([CH3X4])[CH3X4]),"
        "$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),"
        "$([CH2X4][CH2X4][SX2][CH3X4]),"
        "$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),"
        "$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),"
        "$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),"
        "$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),"
        "$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]"
        # "Generic amino acid but not a "standard" amino acid
        # ("standard" refers to the 20 normal side chains).
        # Won't hit amino acids that are non-standard due solely to the fact
        # that groups are terminally-appended to the polypeptide chain (N or C term).
        # format is [$(generic a.a.); !$(not a standard one)]
        # Hits single a.a.s and specific residues w/in polypeptides (internal, or terminal)."
    ),
    # maybe do a loop where we do substructure searches intra-fingerprint
    # any matches get a 'distance' of one
    # additional, unsupported:
    "r_c3": "[#6]~1~[#6]~[#6]~1",
    "r_c4": "[#6]~1~[#6]~[#6]~[#6]~1",
    "r_c5": "[#6]~1~[#6]~[#6]~[#6]~[#6]~1",
    "r_c6": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1",
    "r_c7": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1",
    "r_c8": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1",
    "r_c9": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1",
    "r_c10": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1",
    "sterol": "C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",  # make fuzzy
}


# test
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

    def as_smarts(self):
        return [SUBSTRUCTURES[x] for x in self.name]


SUBS_PATHWAYS = {
    "d_indoleC2N_12": ["shikimate"],
    "d_phenylC2N_9": ["shikimate"],
    # acetate
    "d_c5n_6": ["acetate"],
    "d_c4n_5": ["acetate"],
    # shikimate
    "d_phenylC3_9": ["shikimate"],
    "d_phenylC2_8": ["shikimate"],
    "d_phenylC1_7": ["shikimate"],
    # mevalonate/MEP
    "d_isoprene_5": ["mevalonate", "methylerythritol phosphate"],
    # acetate, aromaticity ok
    "d_ethyl_2": ["acetate"],
    "d_methyl_1": ["acetate"],  # only end methyls
    # sugar-related --------------------------------------------------------
    "s_pyranose_C5O4": ["sugar"],
    "s_furanose_C4O3": ["sugar"],
    "s_openpyr_C6O6": ["sugar"],
    "s_openfur_C5O5": ["sugar"],
    # additional from dewick -----------------------------------------------
    # acetate
    "d2_acetyl_C2O1": ["acetate"],
    "d2_methylmalonyl_C3": ["acetate"],
    # amino acids
    "allstnd_aminos": ["amino acids"],  # can be pathway-specific later on
    "nonstnd_aminos": ["amino acids"],  # can be pathway-specific later on
}

FP_VERSIONS = {
    "fps_full_0814": [
        "fp8",
        "fp12",
        "fp9",
        "fp10",
        "fp11",
        "fp13",
        "fp1",
        "fp2",
        "fp3",
        "fp4",
        "fp5",
        "fp6",
        "fp7",
        "hal1",
        "hal2",
        "hal3",
        "hal4",
        "fp98",
        "fp99",
    ],
    "fps_full_2_0814": [
        "fp8_2",
        "fp12_2",
        "fp9_2",
        "fp10_2",
        "fp11_2",
        "fp13_2",
        "fp1_2",
        "fp2_2",
        "fp3_2",
        "fp4_2",
        "fp5_2",
        "fp6_2",
        "fp7_2",
        "hal1_2",
        "hal2_2",
        "hal3_2",
        "hal4_2",
        "fp98_2",
        "fp99_2",
    ],
    "fps_full_3_0814": [
        "fp8_3",
        "fp12_3",
        "fp9_3",
        "fp10_3",
        "fp11_3",
        "fp13_3",
        "fp1_3",
        "fp2_3",
        "fp3_3",
        "fp4_3",
        "fp5_3",
        "fp6_3",
        "fp7_3",
        "hal1_2",
        "hal2_2",
        "hal3_2",
        "hal4_2",
        "fp98_3",
        "fp99_3_old",
    ],
    "regular_1007": [
        "co_coa",
        "co_nadh",
        "co_nadph",
        "s_pyranose_C5O4",
        "s_furanose_C4O3",
        "s_openpyr_C6O6",
        "s_openfur_C5O5",
        "d_indoleC2N_12",
        "d_phenylC2N_9",
        "d_c5n_6",
        "d_c4n_5",
        "d_phenylC3_9_1007",
        "d_phenylC2_8_1007",
        "d_phenylC1_7_1007",
        "d_isoprene_5",
        "d2_acetyl_C2O1",
        "d2_methylmalonyl_C3",
        "d_ethyl_2",
        "d_methyl_1",
        "phosphate_2",
        "sulfonate_2",
        "hal_f",
        "hal_cl",
        "hal_br",
        "hal_i",
        "n_nitrate_1",
        "o_epoxy_1",
        "o_ether_1",
        "o_hydroxyl_1",
    ],
    "regular_1008": [  # new version with d_phenylC1_7, after isoprene
        "co_coa",
        "co_nadh",
        "co_nadph",
        "s_pyranose_C5O4",
        "s_furanose_C4O3",
        "s_openpyr_C6O6",
        "s_openfur_C5O5",
        "d_indoleC2N_12",
        "d_phenylC2N_9",
        "d_c5n_6",
        "d_c4n_5",
        "d_phenylC3_9",
        "d_phenylC2_8",
        "d_phenylC1_7",
        "d_isoprene_5",
        "d2_acetyl_C2O1",
        "d2_methylmalonyl_C3",
        "d_ethyl_2",
        "d_methyl_1",
        "phosphate_2",
        "sulfonate_2",
        "hal_f",
        "hal_cl",
        "hal_br",
        "hal_i",
        "n_nitrate_1",
        "o_epoxy_1",
        "o_ether_1",
        "o_hydroxyl_1",
    ],
    "strictphenyl_1016": [  # stricter version with d_phenylC1_7_strict (1016)
        "co_coa",
        "co_nadh",
        "co_nadph",
        "s_openpyr_C6O6",  # open ones first
        "s_openfur_C5O5",
        "s_pyranose_C5O4",
        "s_furanose_C4O3",
        "d_indoleC2N_12",
        "d_phenylC2N_9",
        "d_c5n_6",
        "d_c4n_5",
        "d_phenylC3_9_strict",
        "d_phenylC2_8_strict",
        "d_phenylC1_7_strict",
        "d_isoprene_5",
        "d2_acetyl_C2O1",
        "d2_methylmalonyl_C3",
        "d_ethyl_2",
        "d_methyl_1",
        "phosphate_2",
        "sulfonate_2",
        "hal_f",
        "hal_cl",
        "hal_br",
        "hal_i",
        "n_nitrate_1",
        "o_epoxy_1",
        "o_ether_1",
        "o_hydroxyl_1",
    ],
    "aa_1018": [  # includes 'full' amino acids
        "co_coa",
        "co_nadh",
        "co_nadph",
        "allstnd_aminos",
        "nonstnd_aminos",
        "s_openpyr_C6O6",  # open ones first
        "s_openfur_C5O5",
        "s_pyranose_C5O4",
        "s_furanose_C4O3",
        "d_indoleC2N_12",
        "d_phenylC2N_9",
        "d_c5n_6",
        "d_c4n_5",
        "d_phenylC3_9_strict",
        "d_phenylC2_8_strict",
        "d_phenylC1_7_strict",
        "d_isoprene_5",
        "d2_acetyl_C2O1",
        "d2_methylmalonyl_C3",
        "d_ethyl_2",
        "d_methyl_1",
        "phosphate_2",
        "sulfonate_2",
        "hal_f",
        "hal_cl",
        "hal_br",
        "hal_i",
        "n_nitrate_1",
        "o_epoxy_1",
        "o_ether_1",
        "o_hydroxyl_1",
    ],
    "full_1103": [  # includes 'full' amino acids
        "co_coa",
        "co_nadh",
        "co_nadph",
        "allstnd_aminos",
        "nonstnd_aminos",
        "s_openpyr_C6O6",  # open ones first
        "s_openfur_C5O5",
        "s_pyranose_C5O4",
        "s_furanose_C4O3",
        "d_indoleC2N_12",
        "d_phenylC2N_9",
        "d_c5n_6",
        "d_c4n_5",
        "d_phenylC3_9_strict",
        "d_phenylC2_8_strict",
        "d_phenylC1_7_strict",
        "d_isoprene_5",
        "d2_acetyl_C2O1",
        "d2_methylmalonyl_C3",
        "d_ethyl_2",
        "d_methyl_1",
        "phosphate_2",
        "sulfonate_2",
        "hal_f",
        "hal_cl",
        "hal_br",
        "hal_i",
        "n_nitrate_1",
        "o_epoxy_1",
        "o_ether_1",
        "o_hydroxyl_1",
        "r_c3",
        "r_c4",
        "r_c5",
        "r_c6",
        "r_c7",
        "r_c8",
        "r_c9",
        "r_c10",
    ],
    "leaf": [],
}

# --------------------------- substructure sets ---------------------------
SUBSTRUCTURES_old = {
    "fp1": "c1cccc2c1c(CCN)cn2",
    "fp1_2ar": "c1cccc2c1c([#6][#6][#7])cn2",
    "fp1_2": "[#6]1[#6][#6][#6][#6]2[#6]1[#6]([#6][#6][#7])[#6][#7]2",
    "fp1_3": "[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6](~[#6]~[#6]~[#7])~[#6]~[#7]~2",  # fp2  -- Phe.C2N --  shikimate
    "fp2": "c1ccccc1CCN",
    "fp2_2ar": "c1ccccc1[#6][#6][#7]",
    "fp2_2": "[#6]1[#6][#6][#6][#6][#6]1[#6][#6][#7]",
    "fp2_3": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#7]",  # fp3 -- cyclic C5N -- acetylCoA  (Krebs)
    "fp3": "C1CCCCN1",
    "fp3_2": "[#6]1[#6][#6][#6][#6][#7]1",
    "fp3_3": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1",  # fp4 -- cyclic C4N -- acetylCoA (Krebs)
    "fp4": "C1CCCN1",
    "fp4_2": "[#6]1[#6][#6][#6][#7]1",
    "fp4_3": "[#6]~1~[#6]~[#6]~[#6]~[#7]~1",  # fp49-- c6C3 -- shikimate
    "fp49": "c1ccccc1CCC",
    "fp49_2ar": "c1ccccc1[#6][#6][#6]",
    "fp49_2": "[#6]1[#6][#6][#6][#6][#6]1[#6][#6][#6]",
    "fp49_3": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]",  # fp5 -- c6C2 -- shikimate
    "fp5": "c1ccccc1CC",
    "fp5_2ar": "c1ccccc1[#6][#6]",
    "fp5_2": "[#6]1[#6][#6][#6][#6][#6]1[#6][#6]",
    "fp5_3": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]",  # fp6 -- c6C1 -- shikimate
    "fp6": "c1ccccc1C",
    "fp6_2ar": "c1ccccc1[#6]",
    "fp6_2": "[#6]1[#6][#6][#6][#6][#6]1[#6]",
    "fp6_3": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]",  # fp7 -- isoprene -- mevalonic/methylerythriol
    "fp7": "CC(C)CC",
    "fp7_2": "[#6][#6]([#6])[#6][#6]",
    "fp7_3": "[#6]~[#6](~[#6])~[#6]~[#6]",
    "fp7_2_old": "[#6]~[#6]([#6])~[#6][#6]",  # in the stats, biochemically correcter(?why)
    # fp98 -- C2 -- acetylCoA  *'acetyls also end up in aromatic systems', Dewick
    "fp98": "CC",
    "fp98_2": "[#6][#6]",
    "fp98_3": "[#6]~[#6]",  # fp99 -- C -- acetylCoA * methyl only (matches as 1 atom)
    "fp99": "[CH3]",
    "fp99_2": "C",
    "" "fp99_3_old": "[#6]",  # === personal additions
    # phosphate
    "fp8": "O~P(~O)(~O)~O",
    "fp8_2": "O~P(~O)~O",
    "fp8_3": "P~O",  # sulfate group
    "fp12": "O~S(~O)(~O)~O",
    "fp12_2": "O~S(~O)~O",
    "fp12_3": "S~O",  # --- sugars
    # pyranose -- sugar
    "fp9": "C1OCC(O)C(O)C(O)1",
    "fp9_2": "[#6]1[#8][#6][#6]([#8])[#6]([#8])[#6]([#8])1",
    "fp9_3": "[#6]~1~[#8]~[#6]~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~1",  # furanose -- sugar
    "fp10": "C1OCC(O)C(O)1",
    "fp10_2": "[#6]1[#8][#6][#6]([#8])[#6]([#8])1",
    "fp10_3": "[#6]~1~[#8]~[#6]~[#6](~[#8])~[#6](~[#8])~1",
    # open pyranose -- sugar
    "fp11": "C(O)C(O)C(O)C(O)C(O)C(O)",
    "fp11_2": "[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])",
    "fp11_3": "[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])",
    # open furanose -- sugar
    "fp13": "C(O)C(O)C(O)C(O)C(O)",
    "fp13_2": "[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])[#6]([#8])",
    "fp13_3": "[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])~[#6](~[#8])",  # halogens
    "hal1": "F",
    "hal1_2": "[#9]",
    "hal2": "Cl",
    "hal2_2": "[#17]",
    "hal3": "Br",
    "hal3_2": "[#35]",
    "hal4": "I",
    "hal4_2": "[#53]",
    # additionals:
    "acetylO": "CC=O",
    # double and single bond only ?
    "acetylO_2": "[#6][#6]O",
}


# --------------------------- pathway annotation --------------------------------

SUBS_PATHWAYS_old = {
    "fp1": ["shikimate"],
    "fp2": ["shikimate"],
    "fp3": ["acetate"],
    "fp4": ["acetate"],
    "fp49": ["shikimate"],
    "fp5": ["shikimate"],
    "fp6": ["shikimate"],
    "fp7": ["mevalonate", "methylerythritol phosphate"],
    "fp98": ["acetate"],
    "fp99": ["acetate"],
    "fp11": ["sugar"],
    "fp13": ["sugar"],
    "fp9": ["sugar"],
    "fp10": ["sugar"],
}

# --------------------------- get substructure set ----------------------------


def get_smarts(
    fp_version_name: str,
    subs_smarts: dict = SUBSTRUCTURES,
    fp_versions: dict[str, list[str]] = FP_VERSIONS,
) -> list[list[str, str]]:
    """gives list of smarts of substructures of choice
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


def main():
    # get_smarts(DEFAULT_BIOSYNFONI_VERSION, SUBSTRUCTURES, FP_VERSIONS)
    # get_subsset(DEFAULT_BIOSYNFONI_VERSION, SUBSTRUCTURES, FP_VERSIONS)

    # get_subsset('regular_1007')
    return


if __name__ == "__main__":
    main()
