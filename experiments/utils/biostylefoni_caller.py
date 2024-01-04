import os


import matplotlib.pyplot as plt


def set_style() -> None:
    # get path of this script
    script_path = os.path.dirname(os.path.realpath(__file__))
    # parent_path = os.path.dirname(script_path)
    # utils_path = os.path.join(parent_path, "utils")
    utils_path = os.path.join(script_path)
    print(utils_path)
    style_path = os.path.join(utils_path, "biostylefoni.mplstyle")
    # set style
    plt.style.use(style_path)
    return None
