"""
Retrieves annotations for given UniProt IDs.

Extra imports:
    pandas
        https://pandas.pydata.org/docs/getting_started/install.html
    api_funcs
        https://github.com/grtakaha/protein_alignment_tool/blob/main/api_funcs.py
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    parse_args() -> Namespace
    remove_query(pandas.DataFrame) -> pandas.DataFrame
    main()

Command-Line Arguments:
    --infile
    --out_directory
"""

import argparse
import pandas as pd
import api_funcs as af
from helpers import find_path, fasta_to_df, read_alignment

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser()

    # One or more sequences in FASTA format
    parser.add_argument("-i", "--infile", help="Full path of input file.")
    parser.add_argument("-o", "--out_directory", default="./",
                        help="Full path of output directory.")

    return parser.parse_args()

def main(args):
    """
    Parses an input FASTA file and saves UniProt annotations.

        Outputs:
            A collection of files that includes the following:
                Individual annotation files (.ann), one for each unique sequence in the input file.
                One combined annotation file for these sequences (all.ann).
    """

    #args = parse_args()

    infile = find_path(args.infile, "r", "f").replace("\\", "/")
    print(f"Processing sequences from {infile}\n", flush=True)
    proteins = []
    dupe_dict = {} # Dictionary of duplicates that should be skipped
    if infile.endswith(".fasta"):
        infile_df = fasta_to_df(infile)
        #infile_df = remove_query(infile_df)
        for prot in infile_df.index.values:
            # Skips duplicate accessions
            if dupe_dict.get(prot):
                continue
            elif isinstance(infile_df.loc[prot, "Accession"], pd.Series):
                dupe_dict[prot] = 1
                # Only use the first accession that IS NOT the query.
                for acc in infile_df.loc[prot, "Accession"]:
                    if not acc.startswith(">QUERY_"):
                        acc = acc[1:]
                        break
                print(f"Duplicate accessions for ID: {prot}.\nProceeding with {acc}.\n", flush=True)
            else:
                acc = infile_df.loc[prot, "Accession"][1:]
            proteins.append((prot, acc.split(" ")[0]))
    elif (infile.endswith(".clustal") or 
          infile.endswith(".clustal_num") or 
          infile.endswith(".aln-clustal_num")):
        alignment = read_alignment(infile)
        for whole_prot in alignment.proteins:
            # Set to UniProt format the way FASTA files are read in by default.
            prot = whole_prot.split(" ")[0].split("|")[-1]
            # Note that this will not trigger if QUERY is found first. Decide later if I care.
            if dupe_dict.get(prot):
                print(f"Duplicate accessions for ID: {prot}.\nProceeding with {dupe_dict.get(prot)}.\n", flush=True)
            elif not whole_prot.startswith("QUERY_"):
                proteins.append((prot, whole_prot))
                dupe_dict[prot] = whole_prot
    else:
        print("Infile does not have a \".fasta\", \".clustal\", \".aln-clustal_num\" " +
              "or \".clustal_num\" extension.\nAssuming FASTA file.\n", flush=True)
        infile_df = fasta_to_df(infile)
        #infile_df = remove_query(infile_df)
        for prot in infile_df.index.values:
            # Skips duplicate accessions
            if dupe_dict.get(prot):
                continue            
            elif isinstance(infile_df.loc[prot, "Accession"], pd.Series):
                dupe_dict[prot] = 1
                # Only use the first accession that IS NOT the query.
                for acc in infile_df.loc[prot, "Accession"]:
                    if not acc.startswith("QUERY_"):
                        acc = acc[1:]
                        break
                print(f"Duplicate accessions for ID: {prot}.\nProceeding with {acc}.\n", flush=True)
            else:
                acc = infile_df.loc[prot, "Accession"][1:]
            proteins.append((prot, acc.split(" ")[0]))

    out_directory = find_path(args.out_directory, "w", "d").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)

    features_df_list = []
    for prot_names in proteins:
        prot = prot_names[0]
        whole_prot = prot_names[1]
        metadata = af.get_metadata(prot)

        if metadata is not None:
            features = metadata.get("features")
        else:
            features = None

        if features is not None:
            features_df = pd.json_normalize(features)
        else:
            print(f"Could not retrieve features for {prot}...continuing with empty DataFrame\n",
                  flush=True)
            features_df = pd.DataFrame(columns=["type","description",
                                                "evidences","location.start.value",
                                                "location.start.modifier","location.end.value",
                                                "location.end.modifier"])

        features_df.to_csv(f"{out_directory}/{prot}.ann", sep="\t")

        features_df.insert(0, "prot", prot)
        features_df.insert(1, "whole_prot", whole_prot)
        features_df_list.append(features_df)

    all_features_df = pd.concat(features_df_list, axis=0, join='outer')
    all_features_df.to_csv(f"{out_directory}/all.ann", sep="\t")

if __name__ == "__main__":
    args = parse_args()
    main(args)
