import sys, logging
from time import perf_counter
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

from biosynfoni import Biosynfoni


def maccs(mol: Chem.Mol) -> np.array:
    fingerprint = AllChem.GetMACCSKeysFingerprint(mol, nBits=167)
    return np.array(fingerprint)


def morgan(mol: Chem.Mol) -> np.array:
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=False, radius=2, nBits=2048
    )
    return np.array(fingerprint)


def morgan_chiral(mol: Chem.Mol) -> np.array:
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=True, radius=2, nBits=2048
    )
    return np.array(fingerprint)


def rdk_fp(mol: Chem.Mol) -> np.array:
    fingerprint = Chem.RDKFingerprint(mol, fpSize=2048)
    return np.array(fingerprint)


def biosynfoni(mol: Chem.Mol) -> np.array:
    """returns counted fingerprint list"""
    counted_fingerprint = Biosynfoni(
        mol,
        intersub_overlap=True,
        intrasub_overlap=True,
    ).fingerprint
    return np.array(counted_fingerprint)


def biosynfoni_no_overlap(mol: Chem.Mol) -> np.array:
    """returns counted fingerprint list"""
    counted_fingerprint = Biosynfoni(
        mol,
        intersub_overlap=False,
        intrasub_overlap=False,
    ).fingerprint
    return np.array(counted_fingerprint)


def biosynfoni_no_intraoverlap(mol: Chem.Mol) -> np.array:
    """returns counted fingerprint list"""
    counted_fingerprint = Biosynfoni(
        mol,
        intersub_overlap=True,
        intrasub_overlap=False,
    ).fingerprint
    return np.array(counted_fingerprint)


def main():
    sdf = Path(sys.argv[1]).resolve(strict=True)
    suppl = Chem.SDMolSupplier(sdf)

    fp_functions = {
        "maccs": maccs,
        "morgan": morgan,
        "morgan_chiral": morgan_chiral,
        "rdk": rdk_fp,
        "bsf": biosynfoni,
        "bsf_no_overlap": biosynfoni_no_overlap,
        "bsf_no_intraoverlap": biosynfoni_no_intraoverlap,
    }

    logfile = f"{sdf.stem}_times.log"
    with open(logfile, "w") as log:
        log.write()
        log.write(f"fingerprint\ttime for {len(suppl)} mols\n")

    for name, fnx in tqdm(fp_functions.items(), desc="Getting fingerprints"):
        start = perf_counter()
        fps = [fnx(mol) for mol in tqdm(suppl, total=len(suppl), desc=f"{name}")]
        end = perf_counter()

        with open(logfile, "a") as log:
            log.write(f"{name}\t{end - start}\n")

        if len(fps) != len(suppl):
            logging.error(f"fps: {len(fps)} != suppl: {len(suppl)}")

        np.savetxt(f"{sdf.stem}_{name}.csv", np.array(fps), fmt="%i", delimiter=",")

    exit(0)


if __name__ == "__main__":
    main()
