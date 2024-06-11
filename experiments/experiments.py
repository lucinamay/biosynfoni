import argparse, os, sys, time, subprocess



def cli():
    parser = argparse.ArgumentParser(
        description="A pipeline for analysis not made with snakemake for fear of zombies on server"
    )
  
    args = parser.parse_args()

    return args



def converture(raw_sdf, sdf_path):
    """Converts raw sdf to rdkit-readable sdf"""
    pass

def get_bsfs(sdf_path, fp_path):
    """Get all fingerprints"""
    
    cmd = ["biosynfoni", "-s", sdf_path, "-o", fp_path, "-c"]
    # cmd_less_overlap = ["biosynfoni", "-s", sdf_path, "-o", fp_path, "-c", "-l"]
    # cmd_overlap = ["biosynfoni", "-s", sdf_path, "-o", fp_path, "-c", "-o"]

    subprocess.run(cmd)
    subprocess.run(cmd_less_overlap)
    subprocess.run(cmd_overlap)

    pass





def run_concerto(sdf_path, fp_path, overlap_flag: str = None):
    fnx_path = os.path.join(BSF_PATH, "concerto.py")
    cmd = ["python3", fnx_path, sdf_path, "-o", fp_path, "-c"]
    if overlap_flag is not None:
        cmd.append(overlap_flag)
    subprocess.run(cmd)
    return sdf_path


def run_fppullup(fp_path):
    fnx_path = os.path.join(EXP_PATH, "fp_avger.py")
    subprocess.run(["python3", fnx_path, fp_path])
    return fp_path


def run_umap(fp_path, class_path):
    fnx_path = os.path.join(EXP_PATH, "create_umap.py")
    sp = subprocess.run(["python3", fnx_path, fp_path, class_path])
    # check if finished successfully through io or return code
    code = sp.returncode
    if code != 0:
        raise Exception(f"UMAP failed with return code {code}")
    return


def run_tsne(fp_path, class_path):
    fnx_path = os.path.join(EXP_PATH, "tsne.py")
    sp = subprocess.run(["python3", fnx_path, fp_path, class_path])
    code = sp.returncode
    if code != 0:
        raise Exception(f"tSNE failed with return code {code}")
    return


def run_clustermap(fp_path, class_path, out_path):
    fnx_path = os.path.join(EXP_PATH, "clustermap.py")
    sp = subprocess.run(["python3", fnx_path, fp_path, class_path, out_path])
    code = sp.returncode
    if code != 0:
        raise Exception(f"Clustermap failed with return code {code}")
    return


def run_rf(fp_path, class_path, out_path):
    fnx_path = os.path.join(EXP_PATH, "rf.py")
    sp = subprocess.run(["python3", fnx_path, fp_path, class_path, out_path])
    code = sp.returncode
    if code != 0:
        raise Exception(f"Random Forest failed with return code {code}")
    return


def run_confusion_heatmap(cm_path):
    fnx_path = os.path.join(EXP_PATH, "confusion_heatmap.py")
    sp = subprocess.run(["python3", fnx_path, cm_path])
    code = sp.returncode
    if code != 0:
        raise Exception(f"Confusion heatmap failed with return code {code}")
    return


def main():
    # get input sdfs ----------------
    # coconut

    # chebi

    # metacyc

    # zinc

    # get all fingerprints ----------------
    # coconut

    # chebi

    # metacyc

    # zinc

    # get pmi ----------------
    # coconut

    # chebi

    # metacyc

    # zinc

    # get biosynthetic distance ----------------

    # get umap ----------------
    # coconut

    # chebi

    # metacyc

    # zinc

    # get tsne ----------------

    # get clustermap ----------------

    # get random forest ----------------

    # get confusion heatmap ----------------

    # get importances histogram ----------------

    # get fingerprint pullup ----------------





def main():
    # get path to experiments directory and src directory
    global EXP_PATH, BSF_PATH, CWD
    EXP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    BSF_PATH = os.path.abspath(os.path.join(EXP_PATH, os.pardir, "src/biosynfoni"))
    CWD = os.getcwd()

    # coconut input sdf, classifications
    sdf_path = run_converture(raw_sdf, sdf_path)  # makes sdf
    fp_path = run_concerto(sdf_path)  # makes fingerprint files
    run_fppullup(fp_path)  # makes plots
    run_umap(fp_path, class_path, out_path)  # makes plots
    run_tsne(fp_path, class_path, out_path)  # makes plots
    run_clustermap(fp_path, class_path, out_path)  # makes plots
    run_rf(
        fp_path, class_path, out_path
    )  # makes model, confusion matrix, feature importances
    run_confusion_heatmap(cm_path)  # makes plot
    run_importances_hist(importances_path)  # makes plot

    # biosynthetic distance
    run_metacyc_extract()  # makes csv with reaction pairs
    run_biosynthetic_distance()  # makes csv with biosynthetic distances and plots
