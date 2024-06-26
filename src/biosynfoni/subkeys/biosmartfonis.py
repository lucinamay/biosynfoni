substructureSmarts = {
    # dewick's building blocks (d_), conserved aromaticity, free rest ------
    # shikimate
    "d_indoleC2N_12": {
        "name": "indoleC2N",
        "smarts": "c1cccc2c1c(~[#6]~[#6]~[#7])cn2",
        "explanation": "rings have to be aromatic, but not necessarily the rest",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C", "N"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.5,  # 0.0 completely flexible, 0.5 some atoms have defined (non-)aromaticity, 1.0 all atoms have defined (non-)aromaticity
        "ideal_detection_order": 0.0,
    },
    "d_phenylC2N_9": {
        "name": "phenylC2N",
        "smarts": "c1ccccc1[#6][#6][#7]",
        "explanation": "ring has to be aromatic, but not necessarily the rest",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C", "N"],
        "fuzziness": 0.0,
        "correlated_substructures": ["r_c6"],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    # acetate
    "d_c5n_6": {
        "name": "C5N",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1",
        "explanation": "",
        "pathway": ["acetate"],
        "dewick": True,
        "elements": ["C", "N"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "d_c4n_5": {
        "name": "C4N",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#7]~1",
        "explanation": "",
        "pathway": ["acetate"],
        "dewick": True,
        "elements": ["C", "N"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    # shikimate
    "d_phenylC3_9_1007": {
        "name": "phenylC3",
        "smarts": "c1ccccc1~[#6]~[#6]~[#6]",
        "explanation": "",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [
            "r_c6",
            "d_phenylC2_8_1007",
            "d_phenylC1_7_1007",
            "d_phenylC3_9",
            "d_phenylC2_8",
            "d_phenylC1_7",
        ],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "d_phenylC2_8_1007": {
        "name": "phenylC2",
        "smarts": "c1ccccc1~[#6]~[#6]",
        "explanation": "",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [
            "r_c6",
            "d_phenylC1_7_1007",
            "d_phenylC3_9",
            "d_phenylC2_8",
            "d_phenylC1_7",
        ],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "d_phenylC1_7_1007": {
        "name": "phenylC1",
        "smarts": "c1ccccc1~[#6]",
        "explanation": "",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [
            "r_c6",
            "d_phenylC3_9",
            "d_phenylC2_8",
            "d_phenylC1_7",
        ],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "d_phenylC3_9": {
        "name": "phenylC3",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6;!$([r6])]~[#6;!$([r6])]",
        "explanation": "last two C's should not be next to ring",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": ["r_c6", "d_phenylC2_8", "d_phenylC1_7"],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },  # last two should not be next to ring
    "d_phenylC2_8": {
        "name": "phenylC2",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6;!$([r6])]",
        "explanation": "last one should not be next to ring",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": ["r_c6", "d_phenylC1_7"],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },  # last one should not be next to ring
    "d_phenylC1_7": {
        "name": "phenylC1",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6]",
        "explanation": "this one is allowed to be ring, to catch all (cf. d_phenylC3_9, phenylC2_8)",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": ["r_c6", "d_phenylC3_9", "d_phenylC2_8"],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },  # this one is allowed to be ring, to catch all
    "d_phenylC3_9_strict": {
        "name": "phenylC3",
        "smarts": "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]~[#6;!$([r6])]",
        "explanation": "phenyl ring not allowed to be in fused rings to prevent wrong identification of fused isoprenic rings",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [
            "r_c6",
            "d_phenylC2_8_strict",
            "d_phenylC1_7_strict",
        ],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },  # not allowed to be in fused rings
    "d_phenylC2_8_strict": {
        "name": "phenylC2",
        "smarts": "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]",
        "explanation": "phenyl ring not allowed to be in fused rings to prevent wrong identification of fused isoprenic rings",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": ["r_c6", "d_phenylC1_7_strict"],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "d_phenylC1_7_strict": {
        "name": "phenylC1",
        "smarts": "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]",
        "explanation": "no fused rings at all. phenyl ring not allowed to be in fused rings to prevent wrong identification of fused isoprenic rings",
        "pathway": ["shikimate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [
            "r_c6",
            "d_phenylC3_9_strict",
            "d_phenylC2_8_strict",
        ],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },  # no fused rings at all
    # mevalonate/MEP
    "d_isoprene_5": {
        "name": "isoprene",
        "smarts": "[#6]~[#6](~[#6])~[#6]~[#6]",
        "explanation": "does not filter for existance of double bond, to allow detection within fused rings",
        "pathway": ["mevalonate", "methylerythritol_phosphate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    # acetate, aromaticity ok
    "d_ethyl_2": {
        "name": "ethyl",
        "smarts": "[#6]~[#6]",
        "explanation": " detects ethyl's in aromatic rings as well, in accordance to biosynthetic products",
        "pathway": ["acetate"],
        "dewick": True,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": ["d2_methylmalonyl_C3", "d2_acetyl_C2O1"],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "d_methyl_1": {
        "name": "methyl",
        "smarts": "[C;D1;h3]",
        "explanation": "explicitly allows only methyl groups, i.e. a terminal C with 3 H's",
        "pathway": ["acetate"],
        "dewick": True,
        "elements": ["C", "H"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },  # only end methyls
    # sugar-related --------------------------------------------------------
    "s_pyranose_C5O4": {
        "name": "pyranose",
        "smarts": "C~1~[#8]~C~C(~[#8])~C(~[#8])~C(~[#8])~1",
        "explanation": "",
        "pathway": ["sugar"],
        "dewick": False,
        "elements": ["C", "O"],
        "fuzziness": 0.0,
        "correlated_substructures": ["s_openpyr_C6O6"],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "s_furanose_C4O3": {
        "name": "furanose",
        "smarts": "C~1~[#8]~C~C(~[#8])~C(~[#8])~1",
        "explanation": "",
        "pathway": ["sugar"],
        "dewick": False,
        "elements": ["C", "O"],
        "fuzziness": 0.0,
        "correlated_substructures": ["s_openfur_C5O5"],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "s_openpyr_C6O6": {
        "name": "open_pyranose",
        "smarts": "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
        "explanation": "",
        "pathway": ["sugar"],
        "dewick": False,
        "elements": ["C", "O"],
        "fuzziness": 0.0,
        "correlated_substructures": ["s_pyranose_C5O4"],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "s_openfur_C5O5": {
        "name": "open_furanose",
        "smarts": "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
        "explanation": "linear furanose, in open form, oxygen atoms can be either double bonded or single bonded to C's",
        "pathway": ["sugar"],
        "dewick": False,
        "elements": ["C", "O"],
        "fuzziness": 0.0,
        "correlated_substructures": ["s_furanose_C4O3"],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    # additional from dewick -----------------------------------------------
    # acetate
    "d2_acetyl_C2O1": {
        "name": "acetyl",
        "smarts": "[#6]~[#6]~[#8]",
        "explanation": "",
        "pathway": ["acetate"],
        "dewick": False,
        "elements": ["C", "O"],
        "fuzziness": 0.0,
        "correlated_substructures": ["d_ethyl_2"],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },  # CC=O
    "d2_methylmalonyl_C3": {
        "name": "methylmalonyl",
        "smarts": "[#6]~[#6][C;D1;h3]",
        "explanation": "one of the C's has to be a (terminal) methyl group",
        "pathway": ["acetate"],
        "dewick": False,
        "elements": ["C", "H"],
        "fuzziness": 0.0,
        "correlated_substructures": ["d2_acetyl_C2O1", "d_ethyl_2"],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },  # check coverage
    # halogens -------------------------------------------------------------
    "hal_f": {
        "name": "fluorine",
        "smarts": "[#9]",
        "explanation": "any fluorine",
        "pathway": [],
        "dewick": False,
        "elements": ["F"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "hal_cl": {
        "name": "chlorine",
        "smarts": "[#17]",
        "explanation": "a",
        "pathway": [],
        "dewick": False,
        "elements": ["Cl"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "hal_br": {
        "name": "bromine",
        "smarts": "[#35]",
        "explanation": "any bromine",
        "pathway": [],
        "dewick": False,
        "elements": ["Br"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "hal_i": {
        "name": "iodine",
        "smarts": "[#53]",
        "explanation": "any iodine",
        "pathway": [],
        "dewick": False,
        "elements": ["I"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    # phosphates, sulfonates -----------------------------------------------
    "phosphate_2": {
        "name": "phosphate",
        "smarts": "P~O",
        "explanation": "any single or double bonded phosphate with at least one oxygen. does not match additional oxygens",
        "pathway": [],
        "dewick": False,
        "elements": ["P", "O"],
        "fuzziness": 1.0,
        "correlated_substructures": ["phosphate_5"],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },  # if at end
    "phosphate_5": {
        "name": "phosphate",
        "smarts": "O~P(~O)(~O)~O",
        "explanation": "phosphate with 4 oxygens, either single or double bonded",
        "pathway": [],
        "dewick": False,
        "elements": [],
        "fuzziness": 0.0,
        "correlated_substructures": ["phosphate_2"],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },
    "sulfonate_2": {
        "name": "sulfonate",
        "smarts": "S~O",
        "explanation": "sulfur single or double bonded to ('at least') oxygen. does not match additional oxygens",
        "pathway": [],
        "dewick": False,
        "elements": ["S", "O"],
        "fuzziness": 0.0,
        "correlated_substructures": ["sulfonate_5"],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },
    "sulfonate_5": {
        "name": "sulfonate",
        "smarts": "O~S(~O)(~O)~O",
        "explanation": "sulfur with 4 oxygens, either single or double bonded",
        "pathway": [],
        "dewick": False,
        "elements": [],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },
    # additional
    "n_nitrate_1": {
        "name": "nitrate",
        "smarts": "[N;D1]",
        "explanation": "N with one single bond, i.e. terminal N. amount of H's connected does not matter.",
        "pathway": [],
        "dewick": False,
        "elements": ["N"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },
    "o_epoxy_1": {
        "name": "epoxy",
        "smarts": "[O;x2;r3]",
        "explanation": "epoxy oxygen, i.e. oxygen with two separate bonds in a ring of size 3. matches only the oxygen, not the entire ring",
        "pathway": [],
        "dewick": False,
        "elements": ["O"],
        "fuzziness": 0.0,
        "correlated_substructures": ["o_ether_1"],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },  # counts only the oxygen (2ringbond,size3)
    "o_ether_1": {
        "name": "ether",
        "smarts": "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S)]",
        "explanation": "ether oxygen, i.e. oxygen with two separate bonds. does not match esters, does not match C=O's,does not match ethers within rings. O should not have H's",
        "pathway": [],
        "dewick": False,
        "elements": ["O"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },  # not ester,twoconn,noringbond,noH
    "o_hydroxyl_1": {
        "name": "hydroxyl",
        "smarts": "[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(P);!$(S)]",
        "explanation": "hydroxyl group. matches single-bound oxygens with one H or without a h (think dehydrogenated hydroxyl groups), but not if they are part of a C-O, P-O or S-O. does not match =O's",
        "pathway": [],
        "dewick": False,
        "elements": ["O"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.5,  # not sure
        "ideal_detection_order": 0.0,
    },  # OH, O-,  only attached to C/N, no acid hydroxyls, no phosphate/sulfonate,
    # Coenzymes, catalytic units etc ---------------------------------------
    "co_coa_pubchem": {
        "smarts": (
            "SCCNC(~O)CCNC(~O)C(C(C)(C)COP(=O)(~O)OP(=O)"
            "(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O"
        ),
        "name": "coenzyme_a",
        "explanation": "strict coa, no fuzziness, according to definition on pubchem",
        "pathway": [],
        "dewick": False,
        "elements": ["C", "N", "O", "P", "S"],
        "fuzziness": -10.0,
        "correlated_substructures": ["co_coa", "phosphate_5", "phosphate_2"],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "co_coa": {
        "smarts": (
            "SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)"
            "(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O"
        ),
        "name": "coenzyme_a",
        "explanation": "coa with fuzziness, based on definition on pubchem",
        "pathway": [],
        "dewick": False,
        "elements": ["C", "N", "O", "P", "S"],
        "fuzziness": 0.0,
        "correlated_substructures": ["co_coa_pubchem"],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "co_nadh_pubchem": {
        "smarts": (
            "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)"
            "(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"
        ),
        "name": "nadh",
        "explanation": "strict nadh, no fuzziness, according to definition on pubchem",
        "pathway": [],
        "dewick": False,
        "elements": [],
        "fuzziness": -10.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "co_nadh": {
        "smarts": (
            "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)"
            "(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O)~O)~O)~O"
        ),
        "name": "nadh",
        "explanation": "nadh with fuzziness, based on definition on pubchem",
        "pathway": [],
        "dewick": False,
        "elements": [],
        "fuzziness": 4.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "co_nadph_pubchem": {
        "smarts": (
            "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)"
            "OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O"
        ),
        "name": "nadph",
        "explanation": "strict nadph, no fuzziness, according to definition on pubchem",
        "pathway": [],
        "dewick": False,
        "elements": [],
        "fuzziness": -10.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    "co_nadph": {
        "smarts": (
            "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6]"
            "(~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6]"
            "(~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]"
            "~[#7]~[#6]~5~4)~[#7])~O~P(~O)(~O)~O)~O)~O)~O"
        ),
        "name": "nadph",
        "explanation": "nadph with fuzziness, based on definition on pubchem",
        "pathway": [],
        "dewick": False,
        "elements": [],
        "fuzziness": 4.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.5,
        "ideal_detection_order": 0.0,
    },
    # amino acids:
    "allstnd_aminos_wrong": {
        "smarts": (
            "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),"
            "$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),"
            "$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N])]"  # from daylight.com A.A. Template for 20 standard a.a.s
        ),
        "name": "standard_amino_acids",
        "explanation": "all standard amino acids, according to definition on daylight.com's A.A. Template for 20 standard a.a.s, but without replacing the *",
        "pathway": ["amino_acid"],
        "dewick": False,
        "elements": ["C", "N", "O", "S"],
        "fuzziness": 0.0,
        "correlated_substructures": ["nonstnd_aminos"],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },
    "allstnd_aminos": {
        "smarts": (
            "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),"
            "$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),"
            "$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),"
            "$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),"
            "$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),"
            "$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:"
            "[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),"
            "$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),"
            "$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),"
            "$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),"
            "$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),"
            "$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]"  # from daylight.com A.A. Template for 20 standard a.a.s
        ),
        "name": "standard_amino_acids",
        "explanation": "all standard amino acids, according to definition on daylight.com's A.A. Template for 20 standard a.a.s, but without replacing the *",
        "pathway": ["amino_acid"],
        "dewick": False,
        "elements": ["C", "N", "O", "S"],
        "fuzziness": 0.0,
        "correlated_substructures": ["nonstnd_aminos"],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },
    "nonstnd_aminos_wrong": {
        "smarts": (
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
        ),
        "name": "non-standard_amino_acids",
        "explanation": (
            "Wrong version as the * is not replaced. Daylight explanation:"
            "Generic amino acid but not a 'standard' amino acid"
            "('standard' refers to the 20 normal side chains)."
            "Won't hit amino acids that are non-standard due solely to the fact"
            "that groups are terminally-appended to the polypeptide chain (N or C term)."
            "format is [$(generic a.a.); !$(not a standard one)]"
            " Hits single a.a.s and specific residues w/in polypeptides (internal, or terminal)."
        ),
        "pathway": ["amino_acid"],
        "dewick": False,
        "elements": ["C", "N", "O", "S"],
        "fuzziness": 0.0,
        "correlated_substructures": ["allstnd_aminos"],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },
    "nonstnd_aminos": {
        "smarts": (
            "[$([NX3,NX4+][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3])"
            ",$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),"
            "$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:"
            "[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),"
            "$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),"
            "$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),"
            "$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),"
            "$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),"
            "$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[O,N]);"
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
        ),
        "name": "non-standard_amino_acids",
        "explanation": (
            "Generic amino acid but not a 'standard' amino acid"
            "('standard' refers to the 20 normal side chains)."
            "Won't hit amino acids that are non-standard due solely to the fact"
            "that groups are terminally-appended to the polypeptide chain (N or C term)."
            "format is [$(generic a.a.); !$(not a standard one)]"
            " Hits single a.a.s and specific residues w/in polypeptides (internal, or terminal)."
        ),
        "pathway": ["amino_acid"],
        "dewick": False,
        "elements": ["C", "N", "O", "S"],
        "fuzziness": 0.0,
        "correlated_substructures": ["allstnd_aminos"],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },
    # maybe do a loop where we do substructure searches intra-fingerprint
    # any matches get a 'distance' of one
    # additional, unsupported:
    "r_c3": {
        "name": "C3_ring",
        "smarts": "[#6]~1~[#6]~[#6]~1",
        "explanation": "any 3-carbon ring",
        "pathway": [],
        "dewick": False,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "r_c4": {
        "name": "C4_ring",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~1",
        "explanation": "any 4-carbon ring",
        "pathway": [],
        "dewick": False,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "r_c5": {
        "name": "C5_ring",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~1",
        "explanation": "",
        "pathway": [],
        "dewick": False,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "r_c6": {
        "name": "C6_ring",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1",
        "explanation": "",
        "pathway": [],
        "dewick": False,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "r_c7": {
        "name": "C7_ring",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1",
        "explanation": "",
        "pathway": [],
        "dewick": False,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "r_c8": {
        "name": "C8_ring",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1",
        "explanation": "",
        "pathway": [],
        "dewick": False,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "r_c9": {
        "name": "C9_ring",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1",
        "explanation": "",
        "pathway": [],
        "dewick": False,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "r_c10": {
        "name": "C10_ring",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1",
        "explanation": "",
        "pathway": [],
        "dewick": False,
        "elements": ["C"],
        "fuzziness": 0.0,
        "correlated_substructures": [],
        "aromaticity_defined": 0.0,
        "ideal_detection_order": 0.0,
    },
    "sterol": {
        "name": "sterol",
        "smarts": "C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
        "explanation": "unfuzzy sterol, according to definition on pubchem",
        "pathway": [],
        "dewick": False,
        "elements": [],
        "fuzziness": -10.0,
        "correlated_substructures": [""],
        "aromaticity_defined": 1.0,
        "ideal_detection_order": 0.0,
    },  # make fuzzy
}
