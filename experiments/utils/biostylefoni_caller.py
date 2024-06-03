import os, logging


import matplotlib.pyplot as plt


def set_style() -> None:
    """
    Set the style of the plot to biostylefoni
    """
    # get path of this script
    script_path = os.path.dirname(os.path.realpath(__file__))
    # parent_path = os.path.dirname(script_path)
    # utils_path = os.path.join(parent_path, "utils")
    utils_path = os.path.join(script_path)
    logging.debug(f"utils_path: {utils_path}")
    style_path = os.path.join(utils_path, "biostylefoni.mplstyle")
    # set style
    try:
        plt.style.use(style_path)
    except:
        logging.warning("Could not set style to biostylefoni. Prepare for ugly plots")
    return None
