import argparse

from biosynfoni.subkeys import get_version
from biosynfoni.inoutput import outfile_namer


defaultVersion = get_version.defaultVersion
allVersions = get_version.fpVersions.keys()


def _induce_input_type(inputlist: list) -> str:
    """induces the input type from the inputlist"""
    if inputlist[0].endswith(".sdf"):
        return "sdf"
    elif inputlist[0].startswith("InChI="):
        return "inchi"
    elif "." not in str(inputlist[0]):
        return "smiles"
    else:
        return "sdf"


def _handle_outnames(args: argparse.Namespace) -> str:
    outname, inname_root, overlap_allowance = "", "", ""  # init
    if args.output is None:
        if args.repr == "sdf":
            inname_root = args.input[0].split("/")[-1].split(".")[0]
        else:
            inname_root = f"{args.repr}s"
        if args.intrasub_overlap:
            overlap_allowance = "_noblock"
        elif args.intersub_overlap:
            overlap_allowance = "_lessblock"
        else:
            overlap_allowance = ""
        outname_root = f"{inname_root}_{args.version}{overlap_allowance}"
        outname = f"{outfile_namer(outname_root)}.bsf"

    else:
        if "." in args.output:
            outname = args.output
        else:
            outname = f"{args.output}.bsf"
    return outname


def cli():
    parser = argparse.ArgumentParser()

    # required
    parser.add_argument(
        "input",
        metavar="input: molecule(s) / molecule supplier",
        type=str,
        nargs="+",
        help=(
            "molecule representation to get biosynfoni of,"
            "can be sdf file, smiles string(s), or InChI string(s)."
            "will induce type, but recommended to add representation argument"
        ),
    )

    # optional
    parser.add_argument(
        "-r",
        "--repr",
        type=str,
        required=False,
        action="store",
        help="specify the input type of the molecule(s). If molsupplier, choose sdf",
        choices=["sdf", "smiles", "inchi"],
        default="induce",
    )

    parser.add_argument(
        "-v",
        "--version",
        type=str,
        required=False,
        action="store",
        help="specify the fingerprint version to use. If not specified, will use default",
        choices=allVersions,
        default=defaultVersion,
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        action="store",
        help="specify the output file name. If not specified, will use outfile_namer",
        default=None,
    )

    # flags
    parser.add_argument(
        "-l",
        "--intersub_overlap",
        "--lessblocking",
        required=False,
        action="store_true",
        help="allows overlap between different substructures",
        default=False,
    )
    parser.add_argument(
        "-a",
        "--intrasub_overlap",
        "--intrasub",
        required=False,
        action="store_true",
        help="allows overlap between same substructures",
        default=False,
    )

    parser.add_argument(
        "-n",
        "--noblocking",
        "--noblockingstrong",
        required=False,
        action="store_true",
        help="allows all overlap: between different substructures and between same substructures",
        default=False,
    )

    parser.add_argument(
        "-c",
        "--coverage",
        required=False,
        action="store_true",
        help="pass if you want a file with the coverage data",
        default=False,
    )

    parser.add_argument(
        "-p",
        "--nosave",
        "--printonly",
        required=False,
        action="store_true",
        help=(
            "pass if you want to print the fingerprint to stdout."
            "will only print fingerprint, no extra data, unless you pass verbose"
            "default for more than 1 molecule: False"
            "for 1 molecule will be overwritten to True, unless specified save"
        ),
        default=False,
    )

    parser.add_argument(
        "-V",
        "--verbose",
        required=False,
        action="store_true",
        help="pass if you want to get prints in non-saving mode. Defaults to true for saving, false for nosave",
        default=None,
    )

    parser.add_argument(
        "-s",
        "--save",
        required=False,
        action="store_true",
        help=("pass if you want to overwrite standard printonly for 1 mol"),
        default=None,
    )

    # args = vars(parser.parse_args())
    args = parser.parse_args()

    # handle overlap flags
    if args.noblocking:
        args.intersub_overlap = True
        args.intrasub_overlap = True

    if args.repr == "induce":
        args.repr = _induce_input_type(args.input)

    # overwrite for 1 mol if not specified save
    if args.nosave and args.verbose is None:
        args.verbose = False

    if args.save is None:
        if len(args.input) == 1 and not args.repr == "sdf":
            # for a single molecule representation it will only print fp.
            args.save = False
            args.nosave = True
            if args.verbose is None:
                args.verbose = False
        if args.repr == "sdf":
            args.save = True
            args.nosave = False

    if args.save and args.verbose is None:
        args.verbose = True

    args.output = _handle_outnames(args)

    return args
