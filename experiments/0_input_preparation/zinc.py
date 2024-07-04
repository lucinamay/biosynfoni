# -*- coding: utf-8 -*-

import os, subprocess, sys
from pathlib import Path
from functools import partial

from rdkit import Chem

_run = partial(subprocess.run, stdout=subprocess.PIPE, text=True, check=True)


def _get_zid_smiles(zinc_all_path) -> dict:
    cmd1 = ["grep", "-E", "-A3", "ZINC[0-9]{10,}", f"{zinc_all_path}"]
    cmd2 = ["grep", "-E", f"^[^>-].*"]
    print(cmd2)
    grep1 = _run(cmd1)
    grep2 = _run(cmd2, input=grep1.stdout)
    id_smiles = grep2.stdout.splitlines()
    return dict(zip(id_smiles[::2], id_smiles[1::2]))


def _get_zinc_ids(path) -> list:
    assert os.path.exists(path), f"{path} not found"
    cmd = ["grep", "-Eo", "ZINC[0-9]{10,}", f"{path}"]
    zinc_id = _run(cmd).stdout.splitlines()
    return [i for i in zinc_id if i]


def get_synthetics(raw_data_path) -> dict:
    iwd = os.getcwd()
    os.chdir(raw_data_path)
    zid_smiles = _get_zid_smiles("zinc-all-for-sale.sdf")
    nps = _get_zinc_ids("13062018.natural-products.sdf")
    biogens = _get_zinc_ids("12.07.2018.biogenic.smi")
    all_zids = zid_smiles.keys()
    synthetic_zid = set(all_zids) - set(nps + biogens)
    os.chdir(iwd)
    return {k: v for k, v in zid_smiles.items() if k in synthetic_zid}


def write_sdf(synthetics: dict, sdf_out: str):
    writer = Chem.SDWriter(sdf_out)
    for zid, smi in synthetics.items():
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            with open("zinc_sdf.err", "a") as fo:
                fo.write(f"{zid}\t{smi}\n")
            continue
        mol.SetProp("zinc_id", zid)
        writer.write(mol)
    writer.close()
    return None


def main():
    """
    This script is for filtering the synthetic molecules from the ZINC database as provided
    by NaPLeS zenodo files. It reads the synthetic zinc numbers and the zinc numbers and
    smiles from the NaPLeS zenodo files and writes the synthetic smiles to a file.
    """
    raw_data = Path(sys.argv[1]).resolve(strict=True)

    smi_out = "zinc.smi"  # for easier visualisation
    sdf_out = "zinc.sdf"

    synthetics = get_synthetics(raw_data)
    with open(smi_out, "w") as fo:
        for zid, smi in synthetics.items():
            fo.write(f"{zid},{smi}\n")

    write_sdf(synthetics, sdf_out)
    exit(0)


if __name__ == "__main__":
    main()
