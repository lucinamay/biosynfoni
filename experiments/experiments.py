import argparse, os, sys, time, subprocess
from pathlib import Path, PosixPath


# --------------------------------- HELPER -------------------------------------
class ChangeDirectory:
    def __init__(self, new_path):
        self.new_path = new_path
        self.original_path = None

    def __enter__(self):
        self.original_path = os.getcwd()
        os.makedirs(self.new_path, exist_ok=True)
        os.chdir(self.new_path)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.original_path is not None:
            os.chdir(self.original_path)


# -------------------------------------------------------------------------------
def run(cmd):
    process = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return process


def input_preparation(script_path: PosixPath, work_dir: PosixPath) -> int:
    raw_data = work_dir / "raw_data"
    scripts = script_path / "0_input_preparation"
    input_path = work_dir / "input"
    with ChangeDirectory(input_path):
        input_types = ["coconut", "chebi", "metacyc", "zinc"]

        sdf_processes = {}
        for ityp in input_types:
            to_sdf = run(["python3", scripts / f"{ityp}.py", raw_data, "&"])
            sdf_processes[ityp] = to_sdf

        print({ityp: process.pid for ityp, process in sdf_processes.items()})
        for process in sdf_processes.values():
            process.wait()

    with ChangeDirectory(work_dir / "fingerprints"):
        fp_processes = {}
        for ityp in input_types:
            to_fp = run(
                ["python3", scripts / "get_fps.py", input_path / f"{ityp}.sdf", "&"]
            )
            fp_processes[ityp] = to_fp

        print({ityp: process.pid for ityp, process in fp_processes.items()})
        for process in fp_processes.values():
            process.wait()
    return fp_processes.values()[-1].returncode

def analyse(script_path: PosixPath, work_dir: PosixPath) -> int:
    


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
    script_path = Path(__file__).resolve(strict=True).parent
    work_dir = Path.cwd().resolve(strict=True)

    rc = input_preparation(script_path, work_dir)
    if rc != 0:
        raise Exception(f"Input preparation failed with return code {rc}")
    
    # run analyses 
    
    rc = analyse(script_path, work_dir)


# def main():
#     # get path to experiments directory and src directory
#     # global EXP_PATH, BSF_PATH, CWD
#     # EXP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
#     # BSF_PATH = os.path.abspath(os.path.join(EXP_PATH, os.pardir, "src/biosynfoni"))

#     # get current working directory


#     # coconut input sdf, classifications
#     sdf_path = run_converture(raw_sdf, sdf_path)  # makes sdf
#     fp_path = run_concerto(sdf_path)  # makes fingerprint files
#     run_fppullup(fp_path)  # makes plots
#     run_umap(fp_path, class_path, out_path)  # makes plots
#     run_tsne(fp_path, class_path, out_path)  # makes plots
#     run_clustermap(fp_path, class_path, out_path)  # makes plots
#     run_rf(
#         fp_path, class_path, out_path
#     )  # makes model, confusion matrix, feature importances
#     run_confusion_heatmap(cm_path)  # makes plot
#     run_importances_hist(importances_path)  # makes plot

#     # biosynthetic distance
#     run_metacyc_extract()  # makes csv with reaction pairs
#     run_biosynthetic_distance()  # makes csv with biosynthetic distances and plots
