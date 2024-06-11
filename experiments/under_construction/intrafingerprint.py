import sys, os

import numpy as np
import pandas as pd
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(sys.path[0], os.pardir, "src")))
from biosynfoni.inoutput import outfile_namer
from biosynfoni.rdkfnx import get_subsset, defaultVersion


def intramatch(substructures: list[Chem.Mol]) -> list[list[int]]:
    intramatches = []
    for sub_target in substructures:
        sub_target_matches = []
        for sub_query in substructures:
            if sub_target == sub_query:
                continue
            else:
                matches = sub_target.GetSubstructMatches(sub_query)
                if matches:
                    sub_target_matches.append(matches)
        intramatches.append(sub_target_matches)
    return intramatches


def intramatch_nonoverlap():
    """uses get_matches of biosynfoni algorithm, with inter- and or only intra- blocking out"""
    return None


def main():
    substructures = get_subsset(defaultVersion)
    intramatches = np.array(intramatch(substructures))
    intramatches.to_csv(f"{outfile_namer(defaultVersion)}.csv")
