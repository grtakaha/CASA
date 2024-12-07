"""
Downloads current SwissProt release.
BLASTs FASTA sequences against current SwissProt release.
Saves BLAST results (see main() docstring for more info).

Required installations:
    BLAST+ CLI (tested on BLAST 2.12.0+)
        https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Extra imports:
    api_funcs
        https://github.com/grtakaha/protein_alignment_tool/blob/main/api_funcs.py
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Other requirements:
    ~1 GB of storage for SwissProt download and creation of BLAST database

Functions:
    parse_args() -> Namespace
    main()

Command-Line Arguments:
    --infile
    --out_directory
    --stype
    --num_res
"""

import argparse
import api_funcs as af
from helpers import find_path, fasta_to_df

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser(prog="UniProt BLAST script",
                                     description="BLASTs FASTA sequences against current SwissProt release.")

    parser.add_argument("-i", "--infile", help="Full path of input file.")
    parser.add_argument("-o", "--out_directory", default="./",
                        help="Full path of output directory.")
    parser.add_argument("-s", "--stype", default="protein",
                        help="Sequence type (\"protein\" is currently the only option).")
    parser.add_argument("-nr", "--num_res", default="10", help="Number of results.")

    return parser.parse_args()

def main(args):
    """
    Parses an input FASTA file and saves UniProt BLAST results in separate directories.

        Outputs:
            One directory for each sequence (query), each with the following files:
                Table ([QUERY].tsv) and readable ([QUERY].out) BLAST results for that sequence.
                Individual FASTA files with UniProt sequences for each BLAST hit.
                One FASTA file containing all protein sequences, including the query sequence.
    """

    # TODO: Add in translation feature later maybe...or just remove dna.

    infile = find_path(args.infile, "r", "f").replace("\\", "/")
    print(f"Processing sequences from {infile}\n", flush=True)
    infile_df = fasta_to_df(infile)

    out_directory = find_path(args.out_directory, "w", "d").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)

    # Check for Swiss-Prot files in installation path.
    af.verify_sprot()

    # TODO: Add readable results back in. Right now it only outputs outfmt6.
    for protein in infile_df.index.values:
        print(f"BLASTing {protein}...", flush=True)

        sequence = infile_df.loc[protein]["Sequence"]
        accession = infile_df.loc[protein]["Accession"] # Includes ">".

        prot_directory = find_path(f"{out_directory}/{protein}/", "w", "d")
        out_prefix = f"{prot_directory}/{protein}"

        # Saves a FASTA query file.
        query = f"{prot_directory}/{protein}.fasta"
        with open(query, "w", encoding="utf-8") as q_fasta:
            q_fasta.write(f"{accession}\n{sequence}\n")

        af.blast(query, args.stype, f"{out_prefix}", num_res=args.num_res)

        with open(f"{out_prefix}.tsv", "r", encoding="utf-8") as b_res:
            blast_results = b_res.read()

        # Overwrites output all.fasta if it exists.
        with open(f"{prot_directory}/all.fasta", "w", encoding="utf-8") as all_fasta:
            all_fasta.write(f"{accession[0]}QUERY_{accession[1:]}\n{sequence}\n")

        # Parse blast_results (tsv form).
        # Skip last line (empty).
        # No header in command-line blastp outfmt6.
        for line in blast_results.split("\n")[0:-1]:
            hit = line.split("\t")[1]
            print(f"Found BLAST hit: {hit}", flush=True)

            hit_name = hit.split(" ")[0].split("|")[-1]
            hit_fasta = af.get_fasta(hit_name)
            with open(f"{prot_directory}/{hit_name}.fasta", "w", encoding="utf-8") as fasta:
                fasta.write(hit_fasta)
            with open(f"{prot_directory}/all.fasta", "a", encoding="utf-8") as all_fasta:
                all_fasta.write(hit_fasta)

if __name__ == "__main__":
    args = parse_args()
    main(args)
