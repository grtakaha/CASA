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
    --database
    --blast_options
"""

import os
import argparse
import api_funcs as af
from helpers import find_path, fasta_to_df, parse_kwargs

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser(prog="UniProt BLAST script",
                                     description="BLASTs FASTA sequences against a given BLAST database.")

    parser.add_argument("-i", "--infile", help="Full path of input file.")
    parser.add_argument("-o", "--out_directory", default="./",
                        help="Full path of output directory.")
    parser.add_argument("-s", "--stype", default="protein",
                        help="Sequence type (\"protein\" is currently the only option).")
    parser.add_argument("-nr", "--num_res", default="5", help="(optional) Number of results.")
    parser.add_argument("-db", "--database", default=None,
                        help="(optional) Full file path to a protein FASTA file " +
                        "that can be used as a BLAST database. " +
                        "makeblastdb will be run on this file if no BLAST database exists.")
    parser.add_argument("-bopts", "--blast_options", default="[]",
                        help="(optional) Bracketed, comma-separated list of valid blastp input parameters. " +
                        "Valid arguments can be shown via CLI with \"blastp -h\". " +
                        "File locations (like -import_search_strategy) MUST be full file paths (not relative). " +
                        "Example Usage: -bopts \"[-threshold 0,-sorthits 4,-max_hsps 1]\"" )

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

    bopts = parse_kwargs(args.blast_options)

    if not args.database:
        db = None # This is redundant...but just being sure

    else:
        db = find_path(args.database, "r", "f").replace("\\", "/")

    print("Warning: If both -db DATABASE and -bopts \"[-db DATABASE]\"" +
          " were set, the -bopts database will be used.\n",
          flush=True)

    # If the user has set db in kwargs, save here.
    # This will overwrite your previous db,
    # whether or not it was set to None.
    if bopts.get("-db"):
        db = bopts.get("-db")

    else:
        # If the user has not inputted their own database, use Swiss-Prot.
        # This will ONLY execute if the user hasn't specified at least one of the following:
        # -db DATABASE
        # -bopts [-db DATABASE]
        if not db:          
            # Check for Swiss-Prot files in installation path.
            af.verify_sprot()            
            sprot_path = find_path(f"{os.path.abspath(os.path.dirname(__file__))}/SwissProt/",
                                   "w", "d")
            db = f"{sprot_path}/uniprot_sprot.fasta"

    # If neither of the above are true, then use whatever was set as db.
    db = find_path(db, "r", "f").replace("\\", "/")    

    bopts["-db"] = db

    db_fasta = fasta_to_df(db)


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

        # Values in bopts will override num_res and db.
        # This is because bopts must be more intentionally set.
        af.blast(query, args.stype, f"{out_prefix}", num_res=args.num_res, db=db, **bopts)

        with open(f"{out_prefix}.tsv", "r", encoding="utf-8") as b_res:
            blast_results = b_res.read()

        # Overwrites output all.fasta if it exists.
        with open(f"{prot_directory}/all.fasta", "w", encoding="utf-8") as all_fasta:
            all_fasta.write(f"{accession[0]}QUERY_{accession[1:]}\n{sequence}\n")

        # Parse blast_results (tsv form).
        # Skip last line (empty).
        # No header in command-line blastp outfmt6.
        # Keep track of which hits (form: sp|A2WQ39|ANS1_ORYSI) had multiple
        # alignments in BLAST.
        # This will cause problems downstream. Do not allow multiple hits.
        hit_dict = {}
        for line in blast_results.split("\n")[0:-1]:
            hit = line.split("\t")[1]
            print(f"Found BLAST hit: {hit}", flush=True)
            
            if not hit_dict.get(hit):
                hit_name = hit.split(" ")[0].split("|")[-1]
                hit_fasta = af.get_fasta(hit_name, db_fasta) # I'm pretty sure this is >name\nseq\n format
                #hit_fasta = af.get_fasta(hit_name)
                with open(f"{prot_directory}/{hit_name}.fasta", "w", encoding="utf-8") as fasta:
                    fasta.write(hit_fasta)
                with open(f"{prot_directory}/all.fasta", "a", encoding="utf-8") as all_fasta:
                    all_fasta.write(hit_fasta)
                # Set existence of hit in hit_dict
                hit_dict[hit] = 1
            else:
                print(f"Hit already found. It likely has multiple " +
                      "alignments in BLAST results.\nConsider using: -max_hsps 1\n" +
                      "Continuing...\n\n")

if __name__ == "__main__":
    args = parse_args()
    main(args)
