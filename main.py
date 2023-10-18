# from src.concerto_fp import get_biosynfoni
# from src.tests import tsne
import argparse

parser = argparse.ArgumentParser(description="obtain the biosynfoni of molecule(s)")

parser.add_argument(
    "molecule(s)",
    metavar="mol",
    type=str,
    nargs="+",
    help="molecule representation to get biosynfoni of",
)


def main():
    # required arguments -----------------------------------------------------
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="path to molecule supplier, or smiles/InChI string if -s or - is specified.",
    )

    # optionals --------------------------------------------------------------
    # --output--
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        help="Path to output tsv file. If not specified, will use outfile_namer",
    )
    # --representation--
    parser.add_argument(
        "-r",
        "--representation",
        type=str,
        help="specify the input type of the molecule(s). If molsupplier, choose sdf (default)",
        action="store_true",
        default="sdf",
        choices=["sdf", "smiles", "inchi"],  # "sdfgz", "sdfbz2", "sdfzip"],
    )
    # -- coverage --
    parser.add_argument(
        "-c",
        "--coverage",
        type=bool,
        required=False,
        default=False,
        help="specify if you want an output of the coverage of the biosynfoni",
    )

    # flags ------------------------------------------------------------------
    parser.add_argument(
        "--header",
        action="store_true",
        help="Flag to indicate that input file contains header.",
    )


args = parser.parse_args()
print(args)
# print(args.accumulate(args.integers))
