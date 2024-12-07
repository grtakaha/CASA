"""
Manager for protein alignment tools.

Extra imports:
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    parse_args() -> Namespace
    execute_tool(Namespace)
    main()

Command-Line Arguments:
    --infile
    --out_directory
    tool_name
    --order
    --stype
    --email
    --num_res
    --title
    --codes
    --nums
    --uniprot_format
    --annotations
"""

import argparse
import os
import subprocess
from helpers import find_path, fasta_to_df
import clustal_to_svg as cts
import retrieve_annotations as annotate
import alignment as align
import search_proteins as sp

# Using one parser here. Consider using subparsers and set_defaults().
def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser(description='CASA Tool Manager')

    parser.add_argument("-i", "--infile", type=str, help="Full path of input file")
    parser.add_argument("-o", "--out_directory", type=str,
                        help="Full path of output directory.")

    # Will provide a list of 1 or more tools to use with these arguments
    parser.add_argument("-ord", "--order", nargs="+",
                        help="Order of tools to run (blast, annotate, align, svg)." +
                        "Ex. --order align svg")
    parser.add_argument("-s", "--stype", default="protein",
                        help="Sequence type (\"protein\" is currently the only option).")
    #parser.add_argument("-e", "--email", default="")
    parser.add_argument("-nr", "--num_res", default="10",
                        help="Number of results.")
    parser.add_argument("-t", "--title", default="alignment",
                        help="Alignment title ([TITLE].clustal, [TITLE].pim).")
    parser.add_argument("-c", "--codes", default="FALSE")
    parser.add_argument("-n", "--nums", default="FALSE",
                        help="When set to TRUE, includes total residue numbers " +
                        "at the end of each line.")
    parser.add_argument("-u", "--uniprot_format", default="TRUE",
                        help="When set to TRUE, truncates all accessions as if " +
                        "they were UniProt entries.\n" +
                        "Ex. sp|P00784|PAPA1_CARPA -> PAPA1_CARPA")
    # Will annotate if this is provided at all; cannot be nonetype.
    parser.add_argument("-a", "--annotations", default="",
                        help="If an annotation file is provided, it will be " +
                        "used to annotate the resulting SVG files.")
    parser.add_argument("-f", "--features",
                        default="Active site:#0000ff,Disulfide bond:#e27441," +
                        "Propeptide:#9e00f2,Signal:#2b7441",
                        help="A comma-separated list of feature:color pairs to include in SVGs." +
                        "Case sensitive. " +
                        "If features include spaces, the list must be enclosed in quotes." +
                        "If no features should be included, use: -f None" +
                        "The following example is default behavior. " +
                        "Ex. -f \"Active site:#0000ff,Disulfide bond:#e27441," +
                        "Propeptide:#9e00f2,Signal:#2b7441\"")

    return parser.parse_args()

# TODO: Make sure non-UniProt entries don't break annotation.
# TODO: Remove subprocesses in favor of just importing each module.
def execute_tool(args):
    """
    Executes the tool(s) specified in args.

        Parameters:
            args (Namespace): Namespace with command-line arguments.
    """

    if len(args.order) > 1:
        tool_name = "multi"
    else:
        tool_name = args.order[0] # Set to the one tool in order.

    # Execute the selected tool.
    print(f"Running {tool_name}...\n", flush=True)
    if tool_name == "annotate": # annotation command main_tool.py ann ...
        annotate.main(args)
        args.annotations = find_outputs(args)[0] # will change args

    elif tool_name == "blast":
        sp.main(args)

    elif tool_name == "align":
        align.main(args)
        args.infile = find_outputs(args)[0] # .clustal_num file

    elif tool_name == "svg":
        cts.main(args)

    elif tool_name == "multi":
        # List of one or more: blast annotate align svg.
        # Set so that changing args.order doesn't affect things.
        # Every call, full_order will shrink by 1.
        full_order = args.order
        for i, tool in enumerate(full_order):
            print(f"In multi run, executing {tool}...\n", flush=True)
            # This should be ok to run sequentially since args.order is set.
            args.order = [tool] # Run again with just one tool.
            execute_tool(args)
            # Will only continue if there is something after blast
            if tool == "blast" and i+1 < len(full_order):
                # Run everything from here on out using multi and different inputs
                infiles = find_outputs(args)
                print("blast tool executed in a multi run.")
                #print("here2", args.order)
                print(f"Remaining calls ({args.order}) will be executed on {infiles}.\n")
                for in_f in infiles:
                    # Reset to original full_order for each blast result.
                    args.order = full_order[i+1:] # Start after blast call.
                    #print("here3", args.order)
                    #args.tool_name = "multi"
                    args.infile = in_f
                    args.out_directory = "/".join(in_f.split("/")[:-1]) + "/"
                    execute_tool(args)
                # Prevents extra execution of annotate, align, or svg after blast.
                break

# TODO: Consider returning outputs, then printing them separately.
def find_outputs(args):
    """
    Returns output file paths for the tool specified in args.

        Parameters:
            args (Namespace): Namespace with command-line arguments.

        Returns:
            new_inputs (list): List of output file paths for the specified tool.
                               Can be used as infiles for the next tool in a multi run.
    """

    new_inputs = []
    tool_name = args.order[0] # First value in args.order.
    if tool_name == "annotate":
        # It also outputs individual annotation files.
        # This is a TSV of all annotations for a particular FASTA.
        new_inputs.append(f"{args.out_directory}/all.ann".replace("//", "/"))
    elif tool_name == "blast":
        proteins = fasta_to_df(find_path(args.infile, "r", "f"))
        for prot in proteins.index.values:
            # Fastas of blast hits.
            new_inputs.append(f"{args.out_directory}/{prot}/all.fasta".replace("//", "/"))
    elif tool_name == "align":
        # TODO: Either make clustalo output .clustal_num, or change this to .clustal.
        new_inputs.append(f"{args.out_directory}/{args.title}.clustal".replace("//", "/"))
    elif tool_name == "svg":
        pass # Not necessary to return anything because svg is the end of the line.

    return new_inputs

def main(args):
    """
    Executes the given protein alignment tools (blast, annotate, align, and svg).

        Outputs:
            Outputs depend on the tool(s) being executed.
    """

    #args = parse_args()
    execute_tool(args)

if __name__ == "__main__":
    args = parse_args()
    main(args)
