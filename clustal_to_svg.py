"""
Converts a given Clustal alignment to a reformatted Inkscape SVG.

Extra imports:
    pandas
        https://pandas.pydata.org/docs/getting_started/install.html
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    get_max_header(protein_classes.Alignment, str) -> int
    format_alignment(protein_classes.Alignment, str, str) -> list
    get_conserved_coords(protein_classes.Alignment, int, str) -> dict
    get_feature_coords(protein_classes.Alignment, int, str) -> dict
    add_conservation(str, int, int) -> str
    create_feature(protein_classes.Feature, int, int) -> str
    create_tspan(str, str, None) -> str
    split_lines(list, int, int) -> list
    create_svg(protein_classes.Alignment, str, str, str, str)
    parse_args() -> Namespace
    main()

Command-Line Arguments:
    --infile
    --out_directory
    --codes
    --nums
    --uniprot_format
    --annotations
"""

import argparse
import pandas as pd
from helpers import find_path, read_alignment

def get_max_header(alignment, nums="FALSE"):
    """
    Returns the true max header length, given the alignment's maximum and nums status.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            nums (str): If "TRUE", incorporates residue numbers into calculations.

        Returns:
            max_header (int): Maximum length of a protein display name.
    """

    # Stored as the default or manually set max header length.
    max_header = alignment.max_header

    # Accounts for space + numbers up to max digits.
    # len(str(alignment.get_max_length())) = length of the maximum integer residue number.
    # Ex. alignment.max_length = 354; max_header -= 2 + 3
    if nums == "TRUE":
        max_header -= 2 + len(str(alignment.get_max_length()))

    return max_header

def format_alignment(alignment, codes="FALSE", nums="FALSE"):
    """
    Returns a list of lines for an alignment, reformatted for Martin Lab figures.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            codes (str): If "TRUE", includes Clustal identity codes at the bottom of each block.
            nums (str): If "TRUE", includes total residue numbers at the end of each line.

        Returns:
            lines (list): List of lines for an alignment, reformatted for Martin lab figures.
    """

    # TODO: Consider storing this formatted alignment, line by line, in an alignment object

    lines = []
    max_header = get_max_header(alignment, nums=nums)
    counts = {} # count of amino acids for each protein

    # Base sequence length around this first Protein object.
    # Messy, but order of proteins matters for line IDs in create_svg().
    prot1 = list(alignment.proteins.values())[0]
    for i in range(0, len(prot1.sequence), 100):
        for name in alignment.proteins: # Gives full names only.
            lines.append("") # start of a new line
            prot = alignment.proteins[name]

            # Add display name and spaces up to max header length + 2.
            lines[-1] += prot.disp_name[:max_header] + (" "*(max_header-len(prot.disp_name))) + "  "
            seq_100 = prot.sequence[i:i+100]
            for i_10 in range(0, len(seq_100), 10): # add spaces between sequence
                lines[-1] += seq_100[i_10:i_10+10] + " "
            lines[-1] = lines[-1][:-1] # remove last " "
            if nums == "TRUE":
                if counts.get(name) is None:
                    counts[name] = 0
                counts[name] += len(seq_100) - seq_100.count("-")
                lines[-1] += f"  {str(counts[name])}"
        lines.append("") # start of new line
        if codes == "TRUE":
            lines[-1] += " "*(max_header+2) # start of code line
            codes_100 = alignment.codes[i:i+100]
            for i_10 in range(0, len(codes_100), 10):
                lines[-1] += codes_100[i_10:i_10+10] + " "
        lines[-1] = lines[-1][:-1] # remove last " "
        lines.append("") # start of new line between blocks
    return lines

def get_conserved_coords(alignment, lines_per_svg=72, nums="FALSE"):
    """
    Returns a dictionary of SVG numbers and their respective conserved residue coordinates.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            lines_per_svg (int): Max number of lines allowed per SVG.
            nums (str): If "TRUE", includes residue numbers in x_coord calculations.

        Returns:
            conserved_dict (dict): Dictionary of SVG#:conserved_residue_list pairs.
                                   Ex. {0:[(residue, center_x, center_y, i, char_num)]}
    """

    #x_start, y_start = 9.11, 20.274 # in px
    #x_start, y_start = 6.188, 15.206 # In pt
    x_start, y_start = 6.988, 15.000 # In pt, correct    
    #char_width, line_height = 5.5, 11 # In px
    char_width, line_height = 4.125, 8.250 # In pt
    #width_gap = (701.576 - 127*char_width) / (127-1) # in px
    width_gap = (526.182 - 127*char_width) / (127-1) # In pt
    #height_gap = (827.5 - 72*line_height) / (72-1) # in px
    height_gap = (620.625 - 72*line_height) / (72-1) # In pt

    lines_per_block = len(alignment.proteins) + 2

    # int() rounds down, so this will be <= 72.
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block
    max_blocks = int(max_lines / lines_per_block)

    # alignment.conserved_res was set by alignment.set_conserved_res() in main().
    conserved_dict = {} # {svg_num:[(f1, x, y), (f2, x, y), etc.], svg_num:[(f3, x, y), etc.], etc.}
    for res, i in alignment.conserved_res:
        # Note: i, the overall position, is 0-indexed in this.
        # This is different from feature.clust_start, which is 1-indexed.

        # I realize this is all very repetitive...it was the only way I could think about it.
        tot_block_num = int(i / 100) + min(1, i % 100) # 1-indexed

        # Necessary for svg #
        svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks))


        block_num = int(tot_block_num % max_blocks)
        if block_num == 0:
            block_num = max_blocks

        char_num = int(i % 100) # Necessary for x
        if char_num == 0:
            char_num = 100

        # Necessary for x
        # Remove one character, then add it back if this is char_num % 10 is not 0.
        char_num += (int(char_num / 10) - 1 + min(1, char_num % 10) + 2 +
                     get_max_header(alignment, nums=nums))

        #line_num = lines_per_block * (block_num - 1) # necessary for y

        # center_x is the center of the circle.
        # 127 is the number of characters that can fit on each line.
        center_x = x_start + ((char_width / 2) * ((2 * char_num) - 1)) + (width_gap * char_num)

        # center_y is the center of the circle.
        # 72 is the number of lines that can fit in this height.
        center_y = (y_start - line_height/2 - height_gap +
                    ((line_height + height_gap) * lines_per_block*(block_num-1)))

        if not conserved_dict.get(svg_num):
            conserved_dict[svg_num] = []
        # Going to be honest, I can't remember why I returned i and char_num...
        conserved_dict[svg_num].append((res, center_x, center_y, i, char_num))

    return conserved_dict

def get_feature_coords(alignment, lines_per_svg=72, nums="FALSE"):
    """
    Returns a dictionary of SVG numbers and their respective feature coordinates.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            lines_per_svg (int): Max number of lines allowed per SVG.
            nums (str): If "TRUE", includes residue numbers in x_coord calculations.

        Returns:
            feature_dict (dict): Dictionary of SVG#:feature_list pairs.
                                 Ex. {0:[(f1, x, y, prot_name1)]}
    """

    # TODO: Consider redoing x-coord. Rectangle widths are still not exactly aligned.
    #x_start, y_start = 9.11, 20.274 # In px
    #x_start, y_start = 6.188, 15.206 # In pt
    x_start, y_start = 6.988, 15.000 # In pt, correct
    #char_width, line_height = 5.280, 11 # In px
    char_width, line_height = 3.960, 8.250 # In pt
    #width_gap = (701.576 - 127*char_width) / (125 + 2) # In px
    width_gap = (526.182 - 127*char_width) / (125 + 2) # In pt

    #height_gap = (827.5 - 72*line_height) / (72-1) # In px
    height_gap = (620.625 - 72*line_height) / (72-1) # In pt

    lines_per_block = len(alignment.proteins) + 2

    # int() rounds down, so this will be <= 72.
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block
    max_blocks = int(max_lines / lines_per_block)

    # Ex. {svg_num1:[(f1, x, y, prot_name1)], svg_num2:[(f3, x, y, prot_name3)]}
    feature_dict = {}
    # i is the number of this protein in an "ordered" dictionary.
    for i, name in enumerate(alignment.proteins):
        protein = alignment.proteins[name]

        # I realize this is all very repetitive...it was the only way I could think about it.
        for feature in protein.features:
            # If clust_start or clust_end is still set to None, skip feature.
            if not feature.clust_start or not feature.clust_end:
                print(f"{feature.feature_type} in {protein.name} " +
                      "is missing some or all position information.", flush=True)
                print(f"Skipping...\n", flush=True)
                # Skip remaining code for this feature
                continue
            
            # The below is used for Disulfide bonds (2 boxes), propeptides, etc.
            box_list = [] # list of individual positions for features with more than one box

            if feature.clust_start != feature.clust_end:
                if feature.feature_type == "Disulfide bond":
                    box_list.append(feature.clust_start)
                    box_list.append(feature.clust_end)
                    #print(f"{protein.name}: {box_list}, {feature.start}, {feature.end}")
                else:
                    # range() is [)
                    for position in range(feature.clust_start, feature.clust_end+1):
                        box_list.append(position)
            else:
                box_list.append(feature.clust_start)

            for box_position in box_list:
                tot_block_num = int(i / 100) + min(1, i % 100) # 1-indexed
    
                # Necessary for svg #
                svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks))
    
                tot_block_num = int(box_position / 100) + min(1, box_position % 100) # 1-indexed
    
                # Necessary for svg #
                svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks))
                #block_num = int(tot_block_num / (svg_num + 1))
                block_num = tot_block_num - (svg_num * max_blocks)
                char_num = int(box_position % 100) # Necessary for x, 1-indexed
                if char_num == 0:
                    char_num = 100
    
                # Necessary for x
                # Remove one character, then add it back if this is char_num % 10 is not 0.
                # get_max_header() takes nums status into account.
                char_num += (int(char_num / 10) - 1 + min(1, char_num % 10) + 2 +
                             get_max_header(alignment, nums=nums))
    
                line_num = lines_per_block * (block_num - 1) + (i+1) # Necessary for y
    
                x_coord = x_start + width_gap + (char_num - 1)*(char_width + width_gap)
                y_coord = y_start + (line_num - 1)*(line_height + height_gap)
    
                if not feature_dict.get(svg_num):
                    feature_dict[svg_num] = []
                # feature is a protein_classes.Feature object
                # x_coord, y_coord are integers
                # protein.disp_name is a string
                feature_dict[svg_num].append((feature, x_coord, y_coord, protein.disp_name))

    return feature_dict

def add_conservation(res, center_x, center_y):
    """
    Creates a transparent rectangle to annotate the given feature.

        Parameters:
            res (str): Conserved residue, used for color coding.
            center_x (int): x-coordinate for the center of an Inkscape circle.
            center_y (int): y-coordinate for the center of an Inkscape circle.

        Returns:
            circle (str): Circle line to be added to the final SVG.
    """

    color_pos = "3953a4"
    color_neg = "ed1c24"
    color_special = "f7ec13"
    color_hphobe = "3bb54a"
    color_others = "000000"

    color_codes = {"R":color_pos, "H":color_others, "K":color_pos,
                   "D":color_neg, "E":color_neg,
                   "S":color_others, "T":color_others, "N":color_others, "Q":color_others,
                   "C":color_special, "G":color_others, "P":color_others,
                   "A":color_hphobe, "V":color_hphobe, "I":color_hphobe,
                   "L":color_hphobe, "M":color_hphobe, "F":color_hphobe,
                   "Y":color_others, "W":color_hphobe, "default": "000000"}

    #circle = f"<circle cx=\"{center_x}\" cy=\"{center_y}\" r=\"3px\" style=\"fill:#"
    circle = f"<circle cx=\"{center_x}pt\" cy=\"{center_y}pt\" r=\"2.25pt\" style=\"fill:#"
    circle += f"{color_codes.get(res)};stroke-width:0;fill-opacity:1;stroke:rgb(0,0,0)\" />"

    return circle

def create_feature(feature, left_x, top_y, feature_colors):
    """
    Creates a transparent rectangle to annotate the given feature.

        Parameters:
            feature (protein_classes.Feature): Feature object to be anntoated.
                                               Not implemented.
            left_x (int): x-coordinate for the left side of an Inkscape rectangle.
            top_y (int): y-coordinate for the top of an Inkscape rectangle.
            feature_colors (dict): Dictionary of feature colors in hex codes.

        Returns:
            rect (str): Rectangle line to be added to the final SVG.
    """
    #if feature_colors == None:
        #feature_colors = {"Active site":"#0000ff",
                          #"Disulfide bond":"#e27441",
                          #"Propeptide":"#9e00f2",
                          #"Signal":"#2b7441"}
    # TODO: Add in functionality for other features.
    # Create the following for a transparent rectangle.
    # Add it to the svg based on x and y of characters.
    #rectangle = f"<rect x=\"{left_x}\" y=\"{top_y}\" width=\"5.5px\" height=\"11px\" "
    rectangle = f"<rect x=\"{left_x}pt\" y=\"{top_y}pt\" width=\"4.125pt\" height=\"8.250pt\" "
    rectangle += f"style=\"fill:{feature_colors[feature.feature_type]};stroke-width:0;fill-opacity:.2\" />"

    return rectangle

def create_tspan(text, tspan_id, feature=None):
    """
    Creates a tspan SVG line with correct tspan ID, text, and fill.

        Parameters:
            text (str): Text to be added to the SVG tspan.
            tspan_id (str): ID used to identify protein, line, and residue numbers.
            feature (None): If set, will override default coloring with specified feature
                            coloring. Not implemented.

        Returns:
            tspan (str): tspan line to be added to the final SVG.
    """

    # TDOO: Consider making this editable for things like active sites
    # TODO: Maybe give proteins something to override this
    # TODO: Add feature markup (for things like active His) that will override default colors

    color_pos = "3953a4"
    color_neg = "ed1c24"
    color_special = "f7ec13"
    color_hphobe = "3bb54a"
    color_others = "000000"

    color_codes = {"R":color_pos, "H":color_others, "K":color_pos,
                   "D":color_neg, "E":color_neg,
                   "S":color_others, "T":color_others, "N":color_others, "Q":color_others,
                   "C":color_special, "G":color_others, "P":color_others,
                   "A":color_hphobe, "V":color_hphobe, "I":color_hphobe,
                   "L":color_hphobe, "M":color_hphobe, "F":color_hphobe,
                   "Y":color_others, "W":color_hphobe, "default": "000000"}

    code = color_codes.get(text)
    if code is None:
        code = color_codes["default"]

    return f"<tspan id=\"{tspan_id}\" style=\"fill:#{code};\">{text}</tspan>"

def split_lines(lines, lines_per_block, lines_per_svg=72):
    """
    Splits text into groups of 72 (with appropriate breaks - only between blocks).

        Parameters:
            lines (list):           List of lines in an alignment.
                                    Formatted by format_alignment().
            lines_per_block (int):  Total number of proteins in alignment + 2.
            lines_per_svg (int):    Max number of lines allowed per SVG.

        Returns:
            svg_list (list): List of lists of lines that go in each SVG.
    """

    # int() rounds down, so max_lines will be <= 72.
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block

    # List of lists of lines that go in each SVG.
    svg_list = []
    for i in range(0, len(lines), max_lines):
        svg_list.append(lines[i:i+max_lines])

    return svg_list

# TODO: Fix bug with identical display names - ID#s will add for both.
# TODO: Add numbers to "-" and " " so that they have real IDs.
def create_svg(alignment, out_directory, codes="FALSE", nums="FALSE", features=None):
    """
    Creates reformatted alignment Inkscape SVGs from a given Clustal alignment.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            out_directory (str): Directory for storage of reformatted alignments.
            codes (str): If "TRUE", includes Clustal identity codes at the bottom of each block.
            nums (str): If "TRUE", includes total residue numbers at the end of each line.
            features (str): If "TRUE", includes feature annotations.

        Outputs:
            A series of Inkscape SVGs with Martin Lab-formatted alignments.
            SVGs are named 0.svg, 1.svg, etc.
    """

    # Adjust as needed; good enough for now.
    line_height = 1.25
    # In px; I have no idea why this doesn't work with the rectangle y...figure out later
    #x_start, y_start = 9.317, 28.5
    #x_start, y_start = 9.317, 20.000
    x_start, y_start = 6.988, 15.000 # In pt, correct

    # Format lines
    text = format_alignment(alignment, codes=codes, nums=nums) # Returns list of lines.

    lines_per_block = len(alignment.proteins) + 2
    svgs = split_lines(text, lines_per_block) # List of svgs with up to 72 lines each.

    max_header = get_max_header(alignment, nums=nums)
    feature_coords = get_feature_coords(alignment, nums=nums)
    conserved_coords = get_conserved_coords(alignment, nums=nums)

    #relevant_features = ["Active site", "Disulfide bond", "Propeptide", "Signal"] # TODO: Add more valid features.

    aa_counts = {} # Use display names and keep track of aa#.
    pos_counts = {} # Use display names and keep track of overall position counts.
    line_counts = {} # Use display names and keep track of line#.
    codes_count = 0
    space_count = 0
    newline_count = 0
    num_count = 0

    aa_list = ["R", "H", "K", "D", "E", "S", "T",
               "N", "Q", "C", "G", "P", "A", "V",
               "I", "L", "M", "F", "Y", "W"]
    codes_map = {"*":"star", ".":"dot", ":":"colon", " ":"space"}
    num_list = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]

    for i, lines in enumerate(svgs):

        # x should never change, but y will be updated by line.
        x = x_start # In pt or px
        #y = y_start + 8.067333 # In px
        y = y_start + 6.05 # In pt

        # TODO: DECIDE WHETHER OR NOT TO TRY TO CHANGE TSPAN STUFF...LEANING TOWARD REVERTING TO OLDER VERSION
        # Inkscape SVG things.
        svg = "<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n"
        svg += "<svg\n"
        svg += "   height=\"9in\"\n   version=\"1.1\"\n   width=\"7.5in\"\n"
        svg += "   xmlns=\"http://www.w3.org/2000/svg\" "
        svg += "   xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\"\n"
        svg += "   xmlns:svg=\"http://www.w3.org/2000/svg\"\n"
        svg += "   xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\">\n"
        svg += "  <defs />\n"

        # Putting features in first puts them behind the text. Decide later if I like this.
        # Add in features (transparent rectangles).
        unique_features = {}
        if features and feature_coords.get(i):
            for f_info in feature_coords[i]:
                feature = f_info[0]
                disp_name = f_info[3]
                unique_feature_id = f"{disp_name}:{feature.feature_type}:{feature.sequence}"
                if feature.feature_type in features:
                    if not unique_features.get(unique_feature_id):
                        if feature.feature_type == "Disulfide bond":
                            # TODO: Fix this later so that Disulfide bonds are represented by only pairs of cysteines.
                            print(f"{feature.feature_type} found for " +
                                  f"{f_info[3]}, residue {feature.start}{feature.sequence}, " +
                                  f"alignment position {feature.clust_start}. " +
                                  f"Adding to {i}.svg.", flush=True)                        
                        else:
                            print(f"{feature.feature_type} found for " +
                                  f"{f_info[3]}, residue {feature.start}{feature.sequence}, " +
                                  f"alignment position {feature.clust_start}. " +
                                  f"Adding to {i}.svg.", flush=True)
                        unique_features[unique_feature_id] = 1
                    svg += create_feature(f_info[0], f_info[1], f_info[2], features) # feature, x, y
        
        svg += "  <g inkscape:label=\"Layer 1\" inkscape:groupmode=\"layer\" id=\"layer1\">\n"
        svg += "    <rect "
        svg += "style=\"fill:none;stroke-width:1;stroke-linecap:round;stroke-miterlimit:100;"
        svg += "stroke:none;stroke-opacity:1;stroke-dasharray:none\" "
        svg += "id=\"box1\" "
        svg += "width=\"7.5in\" height=\"9in\" "
        svg += f"x=\"{x_start}pt\" y=\"{y_start}pt\" />\n"     
        svg += f"    <text style=\"fill:#000000;line-height:{line_height};"
        #svg += "white-space:pre;font-family:Courier New;font-size:9.2px;font-weight:bold;shape-inside:url(#box1);\" "
        svg += "white-space:pre;font-family:Courier New;font-size:6.9pt;font-weight:bold;shape-inside:url(#box1);\" "
        svg += f"xml:space=\"preserve\">"
        #svg += f"<tspan x=\"{x_start}\">"

        # TODO: Consider removing need to re-split.
        prot_num = 0 # Protein number in list(alignment.proteins)
        prot_names = list(alignment.proteins) # Keys. Full protein names.
        for line in lines:

            # Important to add this back in, or they will not act like lines.
            line += "\n"
            if line == "\n":
                # Add a space if empty so user has something to "see".
                header = " \n"
                current_prot = "line"
            elif line[0] == " ":
                header = line[:max_header+2]
                current_prot = "codes"
            else:
                # Protein name and spaces.
                header = line[:max_header+2]
                #current_prot = header.split(" ")[0]
                current_prot = prot_names[prot_num % len(prot_names)]
                if not aa_counts.get(current_prot):
                    # Start at 0 and increment when valid AA is found.
                    aa_counts[current_prot] = 0
                if not pos_counts.get(current_prot):
                    pos_counts[current_prot] = 0
                prot_num += 1 # Increase number once a protein is found.

            if current_prot not in line_counts:
                line_counts[current_prot] = 0
            line_counts[current_prot] += 1
            
            current_prot_id = current_prot.replace("|", "_")
            #svg += f"<tspan x=\"{x}\" y=\"{y}\" id=\"{current_prot_id}{line_counts[current_prot]}\">" # In px
            svg += f"<tspan x=\"{x}pt\" y=\"{y}pt\" id=\"{current_prot_id}{line_counts[current_prot]}\">" # In pt
            
            svg += f"<tspan>{header}</tspan>"
            for residue in line[max_header+2:]:
                # Most of these elif statements are to make this SVG palatable to Illustrator...
                if residue in aa_list:
                    aa_counts[current_prot] += 1
                    pos_counts[current_prot] += 1
                    res_id = f"{current_prot_id}:{residue}{aa_counts[current_prot]}"
                elif residue == "-":
                    pos_counts[current_prot] += 1
                    res_id = f"{current_prot_id}:dash{pos_counts[current_prot]}"
                elif current_prot == "codes" and codes_map.get(residue):
                    codes_count += 1
                    res_id = f"code{codes_count}:{codes_map.get(residue)}"
                elif residue in num_list:
                    num_count += 1
                    res_id = f"num{num_count}:{residue}"
                elif current_prot != "codes" and residue == " ":
                    space_count += 1
                    res_id = f"space{space_count}"
                elif residue == "\n":
                    newline_count += 1
                    res_id = f"newline{newline_count}"
                else:
                    res_id = residue
                svg += create_tspan(residue, res_id)
            svg += "</tspan>"
            #y += 11.5 # In px
            y += 8.625 # In pt
        #svg += "</tspan>"
        svg += "</text>"
        svg += "</g>"

        # Add in conserved residues (circles).
        if conserved_coords.get(i):
            for cons_info in conserved_coords[i]:
                # residue, center_x, center_y
                svg += add_conservation(cons_info[0], cons_info[1], cons_info[2])

        svg += "</svg>"

        out = f"{out_directory}{i}.svg".replace("\\", "/")
        with open(out, "w", encoding="utf-8") as out_svg:
            print(f"Writing part {i} to {out}.\n", flush=True)
            out_svg.write(svg)

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--out_directory", default="./",
                        help="Full path of output directory.")
    parser.add_argument("-c", "--codes", default="FALSE",
                        help="When set to TRUE, includes Clustal identity " +
                        "codes at the bottom of each block.")
    parser.add_argument("-n", "--nums", default="FALSE",
                        help="When set to TRUE, includes total residue numbers " +
                        "at the end of each line.")
    parser.add_argument("-u", "--uniprot_format", default="FALSE",
                        help="When set to TRUE, truncates all accessions as if " +
                        "they were UniProt entries.\n" +
                        "Ex. sp|P00784|PAPA1_CARPA -> PAPA1_CARPA")
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

def main(args):
    """
    Converts a given Clustal alignment to a reformatted Inkscape SVG.

        Outputs:
            A series of Inkscape SVGs with Martin Lab-formatted alignments.
            SVGs are named 0.svg, 1.svg, etc.
    """

    infile = find_path(args.infile, "r", "f").replace("\\", "/")
    print(f"Processing sequences from {infile} \n", flush=True)

    out_directory = find_path(args.out_directory, "w", "d").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)

    alignment = read_alignment(infile)
    alignment.set_disp_names(uniprot_format=args.uniprot_format)

    annotations = args.annotations
    if annotations != "":
        annotations = pd.read_csv(find_path(annotations, "r", "f"), sep="\t")
    else:
        annotations = pd.DataFrame()
    alignment.add_features(annotations) # Add nothing if annotations aren't set.

    if args.features != "None":
        feature_pairs = args.features.split(",")
        features = {}
        for pair in feature_pairs:
            feature, color = pair.split(":")
            features[feature] = color
    else:
        features = None

    create_svg(alignment, find_path(out_directory, "w", "d"),
               codes=args.codes, nums=args.nums, features=features)

if __name__ == "__main__":
    args = parse_args()
    main(args)
