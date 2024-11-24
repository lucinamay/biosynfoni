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


def morgan_chiral(mol: Chem.Mol) -> list:
    return AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=True, radius=2, nBits=2048
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


def write_fingerprints(sdf: Path):
    suppl = Chem.SDMolSupplier(sdf)
    coverages = [Biosynfoni.get_coverage(mol) for mol in suppl]
    sizes = [mol.GetNumHeavyAtoms() for mol in suppl]
    np.savetxt(f"{sdf.stem}_coverages.csv", coverages, delimiter=",", fmt="%.3f")
    np.savetxt(f"{sdf.stem}_sizes.csv", sizes, delimiter=",", fmt="%d")

    fp_functions = {
        "maccs": maccs,
        "morgan": morgan,
        "rdk": rdk_fp,
        "bsf": biosynfoni,
    }

    for name, fnx in tqdm(fp_functions.items(), desc="Getting fingerprints"):
        with open(f"{sdf.stem}_{name}.csv", "w") as f:
            pass
        with open(f"{sdf.stem}_{name}_times.csv", "a") as f:
            pass
        for mol in tqdm(suppl, desc=name):
            time, fp = time_fingerprint(fnx, mol)
            with open(f"{sdf.stem}_{name}_times.csv", "a") as times:
                times.write(f"{time}\n")
            with open(f"{sdf.stem}_{name}.csv", "a") as fp_file:
                fp_file.write(f"{','.join(map(str, np.array(fp)))}\n")

# def write_similarities(fp_folder: np.array):
#     for fp_file in fp_folder.glob("*.csv"):
#         if len(fp_file.stem.split("_")) > 2:
#             continue
#         if not "chebi" in fp_file.stem:
#             continue
#         fps = np.loadtxt(fp_file, delimiter=",", dtype=int)
#         sims = np.empty((fps.shape[0], fps.shape[0]), dtype="f")  # float32

#         # with open(f"{fp_file.stem}_sim.csv", "w") as f:
#         #     pass

#         # with open(f"{fp_file.stem}_sim.csv", "a") as f:

#         # # if enough memory (if 10 humungous servers are available):
#         # minimum = np.minimum(fps[:, None, :], fps[None, :, :]).sum(axis=2)
#         # maximum = np.maximum(fps[:, None, :], fps[None, :, :]).sum(axis=2)
#         # sims = np.where(maximum != 0, minimum / maximum, -1.0)
#         # np.savetxt(f"{fp_file.stem}_sim.csv", sims, delimiter=",", fmt="%.3f")

#         for i, fp1 in tqdm(
#             enumerate(fps), total=fps.shape[0], desc="similarity matrix", unit="fp"
#         ):
#             # i_fps = np.zeros(len(fps), dtype="f")
#             # if all(fp1 == 0):
#             #     continue
#             # i_fps[i] = 1.0
#             # minimum = np.sum(np.minimum(fp1, fps[i:]), axis=1)
#             maximum = np.sum(np.maximum(fp1, fps[i:]), axis=1)
#             sims[i][i:] = np.where(
#                 maximum != 0, np.sum(np.minimum(fp1, fps[i:]), axis=1) / maximum, -1.0
#             )

#         np.savetxt(f"{fp_file.stem}_sim.csv", sims, delimiter=",", fmt="%.3f")

#         # minimum_sum = np.sum(minimum, axis=1)
#         # maximum_sum = np.sum(maximum, axis=1)
#         # i_fps[i:] = minimum_sum / maximum_sum
#         # for j, fp2 in enumerate(fps[i:], start=i):
#         #     if np.sum(fp1 | fp2) == 0:
#         #         return -1.0
#         #     else:
#         #         i_fps[j] = np.sum(np.minimum(fp1, fp2)) / np.sum(
#         #             np.maximum(fp1, fp2)
#         #         )
#         # # append to file
#         # np.savetxt(
#         #     f,
#         #     i_fps,
#         #     delimiter=",",
#         #     fmt="%.3f",
#         # )


def main():
    sdf_folder = Path(sys.argv[1]).resolve(strict=True)
    fp_folder = sdf_folder.parent.parent / "fps"

    for sdf in sdf_folder.glob("*.sdf"):
        with ChangeDirectory(fp_folder):
            write_fingerprints(sdf)

    return 0


if __name__ == "__main__":
    main()
