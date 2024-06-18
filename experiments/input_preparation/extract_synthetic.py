# -*- coding: utf-8 -*-

import os, subprocess
from functools import partial

from rdkit import Chem

_run = partial(subprocess.run, stdout=subprocess.PIPE, text=True, check=True)


def _get_zid_smiles() -> dict:
    cmd1 = ["grep", "-E", "-A3", "ZINC[0-9]{10,}", "zinc-all-for-sale.sdf"]
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


def get_synthetics() -> dict:
    nps = _get_zinc_ids("13062018.natural-products.sdf")
    biogens = _get_zinc_ids("12.07.2018.biogenic.smi")
    zid_smiles = _get_zid_smiles()
    all_zids = zid_smiles.keys()
    synthetic_zid = set(all_zids) - set(nps + biogens)
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
    # assert directory "raw_data" exists (cwd should be data)
    assert os.path.exists("raw_data"), "raw_data directory not found"

    smi_out = "zinc.smi"  # for easier visualisation
    sdf_out = "zinc.sdf"

    os.chdir("raw_data")
    synthetics = get_synthetics()
    os.chdir("..")
    with open(smi_out, "w") as fo:
        for zid, smi in synthetics.items():
            fo.write(f"{zid},{smi}\n")

    write_sdf(synthetics, sdf_out)
    exit(0)


if __name__ == "__main__":
    main()
