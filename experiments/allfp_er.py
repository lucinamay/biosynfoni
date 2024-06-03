import sys, os, logging

import numpy as np

sys.path.append(os.path.abspath(os.path.join(sys.path[0], os.pardir, "src")))
# for intra-biosynfoni-code running
sys.path.append(
    os.path.abspath(os.path.join(sys.path[0], os.pardir, "src", "biosynfoni"))
)
from biosynfoni.inoutput import outfile_namer
from biosynfoni.rdkfnx import get_supplier
from biosynfoni import fingerprints as fp

FP_FUNCTIONS = {
    # "biosynfoni": fp.biosynfoni_getter,
    # "binosynfoni": fp.binosynfoni_getter,
    "maccs": fp.maccs_getter,
    "morgan": fp.morgan_getter,
    "rdkit": fp.rdk_fp_getter,
    # "maccsynfoni": fp.maccsynfoni_getter,
    # "bino_maccs": fp.bino_maccs_getter,
}

sdf = sys.argv[1]
sdf_name = os.path.basename(sdf).split(".")[0]
suppl = get_supplier(sdf_file=sdf)
fp_dict = {x: [] for x in FP_FUNCTIONS.keys()}

i = 0
j = 0
logging.info(f"starting with the {len(suppl)} molecules...")
for mol in suppl:
    if i % 1000 == 0:
        print(f"{j} k done")
        j += 1
    i += 1
    for fp_name, fp_function in FP_FUNCTIONS.items():
        fingerprint = fp_function(mol)
        if fingerprint.any() != None:
            fp_dict[fp_name].append(fingerprint)
        else:
            if fp_name == "morgan" or fp_name == "rdkit":
                fp_dict[fp_name].append(np.nan)
            else:
                fp_dict[fp_name].append(np.nan)

for fp_name, fp_list in fp_dict.items():
    fp_arr = np.array(fp_list)
    outfile = f"{outfile_namer(sdf_name, fp_name)}.csv"
    np.savetxt(outfile, fp_arr, fmt="%i", delimiter=",")
    logging.info(f"saved {fp_name} to {outfile}")
