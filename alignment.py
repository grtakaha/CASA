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
    --title
    --clustal_options
"""

import argparse
import api_funcs as af
from helpers import find_path, parse_kwargs

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser()

    # Overrides -i, --in, --infile clustalo
    parser.add_argument("-i", "--infile", help="Full path of input file.")
    
    # Used to override output directories in multi runs for the following options:
    # -o, --out, --outfile
    # --distmat-out
    # --guidetree-out
    # --clustering-out
    # --posterior-out
    # -l, --log
    # Ex. --distmat-out f"{out_directory}/{title}.pim
    parser.add_argument("-o", "--out_directory", default="./",
                        help="Full path of output directory. Must end with \"/\".")
    parser.add_argument("-s", "--stype", default="protein",
                        help="Sequence type (\"protein\" or \"dna\").")
    parser.add_argument("-t", "--title", default="alignment",
                        help="Alignment title ([TITLE].clustal, [TITLE].pim).")
    parser.add_argument("-copts", "--clustal_options", default="[--full,--outfmt clu,--force]",
                        help="(optional) Bracketed, comma-separated list of valid clustalo input parameters. " +
                        "Valid arguments can be shown via CLI with \"clustalo -h\"." +
                        "File locations (like --hmm-in) MUST be full file paths (not relative). " +
                        "Example Usage: -copts \"[--residuenumber,--iterations 3]\"" )
    
    # Instead, have no control over where it saves for a multi run (or just in general) always save in the same location as the output directory
    
    #subprocess.run(["clustalo",
                    #"--infile", infile,
                    #"--outfile", f"{out_directory}/{title}.clustal",
                    #"--seqtype", stype_conv[stype],
                    #"--distmat-out", f"{out_directory}/{title}.pim",
                    #"--percent-id", "--full",
                    #"--outfmt", "clu", "--force"], check=True)    

    return parser.parse_args()

# NOTE: Removed PIM output by default. Was too annoying to save it to the correct location, 
# --distmat-out
# --guidetree-out
# --clustering-out
# --posterior-out
# --log

# --wrap

def main(args):
    """
    Aligns an input FASTA file with Clustal Omega.

        Outputs:
            An alignment (.clustal_num)
            created from the input FASTA file.
    """

    infile = find_path(args.infile, "r", "f").replace("\\", "/")
    print(f"Processing sequences from {infile} \n", flush=True)

    out_directory = find_path(args.out_directory, "w", "d").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)

    copts = parse_kwargs(args.clustal_options)

    stype = args.stype
    if stype not in ["dna", "protein"]:
        print("Given stype not found. Please specify " +
              "\"-stype dna\" OR \"-stype protein\n", flush=True)
        print("Defaulting to stype=\"protein\"\n", flush=True)
        stype = "protein"

    title = args.title

    # Run the alignment via command-line subprocess
    af.align(infile, stype, out_directory, title, **copts)

if __name__ == "__main__":
    args = parse_args()
    main(args)
