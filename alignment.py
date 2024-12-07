"""
Aligns a given FASTA file using Clustal Omega.

Extra imports:
    api_funcs
        https://github.com/grtakaha/protein_alignment_tool/blob/main/api_funcs.py
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    parse_args() -> Namespace
    main()

Command-Line Arguments:
    --infile
    --out_directory
    --stype
    --email
    --title
"""

import argparse
import api_funcs as af
from helpers import find_path

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", help="Full path of input file.")
    parser.add_argument("-o", "--out_directory", default="./",
                        help="Full path of output directory. Must end with \"/\".")
    parser.add_argument("-s", "--stype", default="protein",
                        help="Sequence type (\"protein\" or \"dna\").")
    parser.add_argument("-t", "--title", default="alignment",
                        help="Alignment title ([TITLE].clustal, [TITLE].pim).")

    return parser.parse_args()

def main(args):
    """
    Aligns an input FASTA file with Clustal Omega.

        Outputs:
            An alignment (.clustal_num) and percent identity matrix (.pim)
            created from the input FASTA file.
    """

    infile = find_path(args.infile, "r", "f").replace("\\", "/")
    print(f"Processing sequences from {infile} \n", flush=True)

    out_directory = find_path(args.out_directory, "w", "d").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)

    stype = args.stype
    if stype not in ["dna", "protein"]:
        print("Given stype not found. Please specify " +
              "\"-stype dna\" OR \"-stype protein\n", flush=True)
        print("Defaulting to stype=\"protein\"\n", flush=True)
        stype = "protein"

    title = args.title

    # Run the alignment via command-line subprocess
    af.align(infile, stype, out_directory, title)

if __name__ == "__main__":
    args = parse_args()
    main(args)
