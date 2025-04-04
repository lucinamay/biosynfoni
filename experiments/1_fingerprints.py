import sys, logging
from time import perf_counter
from pathlib import Path
from typing import Callable

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

from biosynfoni import Biosynfoni
from biosynfoni.fingerprints import counted_tanimoto_sim
from helper import ChangeDirectory


def maccs(mol: Chem.Mol) -> np.array:
    return AllChem.GetMACCSKeysFingerprint(mol)  # , nBits=167)


def morgan(mol: Chem.Mol) -> list:
    return AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=False, radius=2, nBits=2048
    )


def rdk_fp(mol: Chem.Mol) -> list:
    return Chem.RDKFingerprint(mol, fpSize=2048)


def biosynfoni(mol: Chem.Mol) -> list:
    """returns counted fingerprint list"""
    return Biosynfoni(
        mol,
        intersub_overlap=True,
        intrasub_overlap=True,
    ).fingerprint


def biosynfoni_no_overlap(mol: Chem.Mol) -> list:
    """returns counted fingerprint list"""
    return Biosynfoni(
        mol,
        intersub_overlap=False,
        intrasub_overlap=False,
    ).fingerprint


def biosynfoni_no_intraoverlap(mol: Chem.Mol) -> list:
    """returns counted fingerprint list"""
    return Biosynfoni(
        mol,
        intersub_overlap=True,
        intrasub_overlap=False,
    ).fingerprint


def time_fingerprint(fp_function, mol) -> tuple:
    if not mol:
        return 0, []
    start = perf_counter()
    fp = fp_function(mol)
    end = perf_counter()
    return end - start, fp


def write_coverages(sdf: Path):
    suppl = Chem.SDMolSupplier(sdf)
    coverages = [Biosynfoni(mol).get_coverage() for mol in suppl]
    np.savetxt(f"{sdf.stem}_coverages.csv", coverages, delimiter=",", fmt="%.3f")


def write_fingerprints(sdf: Path):
    suppl = Chem.SDMolSupplier(sdf)

    fp_functions = {
        "maccs": maccs,
        "morgan": morgan,
        "rdk": rdk_fp,
        "bsf": biosynfoni,
    }

    for name, fnx in tqdm(fp_functions.items(), desc="Getting fingerprints"):
        with open(f"{sdf.stem}_{name}.csv", "w") as f:
            pass
        with open(f"{sdf.stem}_{name}_times.csv", "w") as f:
            pass
        for mol in tqdm(suppl, desc=name):
            time, fp = time_fingerprint(fnx, mol)
            with open(f"{sdf.stem}_{name}_times.csv", "a") as times:
                times.write(f"{time}\n")
            with open(f"{sdf.stem}_{name}.csv", "a") as fp_file:
                fp_file.write(f"{','.join(map(str, np.array(fp)))}\n")


def main():
    sdf_folder = Path(sys.argv[1]).resolve(strict=True)
    fp_folder = sdf_folder.parent.parent / "fps"

    for sdf in sdf_folder.glob("*.sdf"):
        with ChangeDirectory(fp_folder):
            write_fingerprints(sdf)
            # write_coverages(sdf)

    return 0


if __name__ == "__main__":
    main()
