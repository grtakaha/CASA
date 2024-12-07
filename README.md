# Protein Alignment Tool

A centralized tool manager for four scripts used to process and visualize protein sequences.

**NOTE**: These scripts make use of EMBL-EBI and NCBI resources. References for tools and databases used here include:

UniProt:\
The UniProt Consortium. \
“UniProt: The Universal Protein Knowledgebase in 2023.” \
Nucleic Acids Research 51, no. D1 (January 6, 2023): D523–31. https://doi.org/10.1093/nar/gkac1052.

NCBI:\
Sayers, Eric W, Evan E Bolton, J Rodney Brister, Kathi Canese, Jessica Chan, Donald C Comeau, Ryan Connor, et al. \
“Database Resources of the National Center for Biotechnology Information.” \
Nucleic Acids Research 50, no. D1 (December 1, 2021): D20–26. https://doi.org/10.1093/nar/gkab1112.

Clustal Omega:\
Sievers, Fabian, Andreas Wilm, David Dineen, Toby J Gibson, Kevin Karplus, Weizhong Li, Rodrigo Lopez, et al. \
“Fast, Scalable Generation of High‐quality Protein Multiple Sequence Alignments Using Clustal Omega.” \
Molecular Systems Biology 7, no. 1 (January 2011): 539. https://doi.org/10.1038/msb.2011.75.

Sievers, Fabian, and Desmond G. Higgins. \
“Clustal Omega for Making Accurate Alignments of Many Protein Sequences.” \
Protein Science: A Publication of the Protein Society 27, no. 1 (January 2018): 135–45. https://doi.org/10.1002/pro.3290.

BLAST+:\
Camacho, Christiam, George Coulouris, Vahram Avagyan, Ning Ma, Jason Papadopoulos, Kevin Bealer, and Thomas L. Madden. \
“BLAST+: Architecture and Applications.” \
BMC Bioinformatics 10, no. 1 (December 2009): 1–9. https://doi.org/10.1186/1471-2105-10-421.

https://blast.ncbi.nlm.nih.gov/doc/blast-help/references.html#references

### Prerequisites
General:
* Internet connection (when running search_proteins.py and retrieve_annotations.py)
* Python 3.7+
* ~ 1 GB of storage for Swiss-Prot BLAST database (if running search_proteins.py)

Command-line tools:
* Clustal Omega (http://www.clustal.org/omega/)
* NCBI BLAST+ (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

Python Libraries:
* pandas
* requests

### Installing
Before running, ensure that required command-line tools are on your PATH.
* Clustal Omega is required for alignment.py
* NCBI BLAST+ is required for search_proteins.py

Download and add the protein_alignment_tool directory to PATH and PYTHONPATH.

## SCRIPTS

### pat_main.py

Runs one or more of the below scripts in the order given.

**INPUT**: Depends on which tool(s) are being executed. Should be an acceptable input of the first tool executed. Inputs for runs that start with annotate.py may additionally be limited by what can be passed to downstream tools.  Example: --order annotate svg -> must use .clustal or .clustal_num file as input (.fasta file cannot be used by clustal_to_svg.py)

**OUTPUT**: Depends on which tool(s) are being executed. Each tool will have its own output if it is included in a run. All outputs will be split into separate directories in a run that includes search_proteins.py.

**NOTE**: There are currently limited ways to run multiple tools at once (inputs and outputs will vary depending on start and end):
* blast annotate align svg
* blast annotate align
* blast annotate
* blast align annotate svg
* blast align annotate
* blast align svg
* blast align
* blast
* annotate align svg
* annotate align
* annotate svg
* annotate
* align annotate svg
* align annotate
* align svg
* align
* svg

**Example**: python -m pat_main -i ./alignment.clustal -o ./svg_folder/ -ord annotate svg -u TRUE -c FALSE -nums TRUE

**Usage**:
```
usage: pat_main.py [-h] [-i INFILE] [-o OUT_DIRECTORY]
                   [-ord ORDER [ORDER ...]] [-s STYPE] [-nr NUM_RES]
                   [-t TITLE] [-c CODES] [-n NUMS] [-u UNIPROT_FORMAT]
                   [-a ANNOTATIONS] [-f FEATURES]

Protein Alignment Tool Manager

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Full path of input file
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory.
  -ord ORDER [ORDER ...], --order ORDER [ORDER ...]
                        Order of tools to run (blast, annotate, align,
                        svg).Ex. --order align svg
  -s STYPE, --stype STYPE
                        Sequence type ("protein" is currently the only
                        option).
  -nr NUM_RES, --num_res NUM_RES
                        Number of results.
  -t TITLE, --title TITLE
                        Alignment title ([TITLE].clustal, [TITLE].pim).
  -c CODES, --codes CODES
  -n NUMS, --nums NUMS  When set to TRUE, includes total residue numbers at
                        the end of each line.
  -u UNIPROT_FORMAT, --uniprot_format UNIPROT_FORMAT
                        When set to TRUE, truncates all accessions as if they
                        were UniProt entries. Ex. sp|P00784|PAPA1_CARPA ->
                        PAPA1_CARPA
  -a ANNOTATIONS, --annotations ANNOTATIONS
                        If an annotation file is provided, it will be used to
                        annotate the resulting SVG files.
  -f FEATURES, --features FEATURES
                        A comma-separated list of feature:color pairs to
                        include in SVGs.Case sensitive. If features include
                        spaces, the list must be enclosed in quotes.If no
                        features should be included, use: -f NoneThe following
                        example is default behavior. Ex. -f "Active
                        site:#0000ff,Disulfide
                        bond:#e27441,Propeptide:#9e00f2,Signal:#2b7441"
```

### search_proteins.py

Takes one or more protein sequences (FASTA format) as input and BLASTs them against the current Swiss-Prot release (uniprotkb_refprotswissprot).\

Requires ~1 GB of storage for Swiss-Prot download and creation of BLAST database.

**Swiss-Prot download information**:
Current Swiss-Prot release is verified/downloaded from:
* ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/

Swiss-Prot database and associated files are stored in the same directory as this script.

**INPUT**: FASTA-formatted file with at least one sequence.

**OUTPUT**: A set of directories - one for each sequence in the original input file - that contain the following:
* the BLAST results for that sequence (the query) against the current Swiss-Prot release in table form ([QUERY].tsv)
* individual FASTA files with UniProt sequences for each BLAST hit
* one FASTA file containing all protein sequences, including the query sequence (all.fasta)

**NOTE**: --stype dna is currently not supported in any form. May enter an infinite loop. Please do not use --stype dna until updated.

If used at the beginning of a multi-step run, downstream commands will be run on each resulting collection of outputs.

For example: search_proteins.py will yield multiple all.fasta (one for each query), which can be sent to both retrieve_annotations.py and alignment.py. This is why "blast annotate align" is a valid input for the --order optional argument.

**Example from pat_main.py**:
```
  python -m pat_main -i ./unknown_proteins.fasta -o ./unknown_protein_folder/ -ord blast align -s protein -nr 3
```

**Example standalone**:
```
  python -m search_proteins -i ./unknown_proteins.fasta -o ./ -s protein -nr 5
```

**Usage**:
```
usage: UniProt BLAST script [-h] [-i INFILE] [-o OUT_DIRECTORY] [-s STYPE]
                            [-nr NUM_RES]

BLASTs FASTA sequences against current Swiss-Prot release.

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Full path of input file.
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory.
  -s STYPE, --stype STYPE
                        Sequence type ("protein" is currently the only
                        option).
  -nr NUM_RES, --num_res NUM_RES
                        Number of results.
```

### retrieve_annotations.py

Takes one or more protein sequences (FASTA format) as input and retrieves annotations for sequences whos IDs exist in UniProt. 

**INPUT**: A FASTA-formatted, CLUSTAL_NUM-formatted, or CLUSTAL-formatted file with at least one protein sequence.

**OUTPUT**: A collection of files that includes the following:
* individual annotation files (.ann), one for each unique sequence in the input file
* one combined annotation file that includes all annotations for this collection of sequences (all.ann)

**NOTE**: all.ann can be used as input for clustal_to_svg.py. See **ANNOTATION FORMAT** for help formatting annotations by hand.

**Example from pat_main.py**:
```
python -m pat_main -i ./uniprot_proteins.fasta -o ./annotations_folder/ -ord annotate align svg -s protein -u TRUE -c FALSE -n FALSE
```

**Example standalone**:
```
python -m retrieve_annotations -i ./uniprot_proteins.fasta -o ./annotations_folder/
```

**Usage**:
```
usage: retrieve_annotations.py [-h] [-i INFILE] [-o OUT_DIRECTORY]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Full path of input file.
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory.
```

### alignment.py

Takes at least two protein sequences as input and aligns them using Clustal Omega.

**INPUT**: A FASTA-formatted file with at least two sequences

**OUTPUT**: An alignment (.clustal) and associated percent identity matrix (.pim) of the given FASTA file

**NOTE**: A percent identity matrix will only be generated by Clustal Omega for an input of three or more sequences.

**Example from pat_main.py**:
```
python -m pat_main -i ./three_proteins.fasta -o ./alignments_folder/ -ord annotate align -s protein -title three_proteins
```

**Example standalone**:
```
python -m alignment -i ./three_proteins.fasta -o ./alignments_folder/ -s protein -title three_proteins
```

**Usage**:
```
usage: alignment.py [-h] [-i INFILE] [-o OUT_DIRECTORY] [-s STYPE] [-t TITLE]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Full path of input file.
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory. Must end with "/".
  -s STYPE, --stype STYPE
                        Sequence type ("protein" or "dna").
  -t TITLE, --title TITLE
                        Alignment title ([TITLE].clustal, [TITLE].pim).
```

### clustal_to_svg.py

Reformats a .clustal_num or .clustal alignment into an editable Inkscape SVG. Currently annotates conserved residues (automatic, not optional) and a given list of features (optional).

**INPUT**: A CLUSTAL or CLUSTAL_NUM file

**OUTPUT**: A sequential set of SVGs (.svg), numbered 0, 1, 2, etc., with formatted alignments and associated conserved residues and/or annotations.

**NOTE**: These SVG outputs were designed to be edited in Inkscape. They retain full functionality when opened in Inkscape, including multiline text-box editability. They retain some functionality when opened in other SVG viewers/editors. They can be viewed and edited in Illustrator, but do not retain multiline text-box editability (each letter is its own text-box). They can also be viewed in Chrome.

**Example from pat_main.py**:
```
python -m pat_main -i ./proteins.fasta -o ./SVGs/ -ord align svg -s protein -u TRUE -c FALSE -n FALSE -a annotations.ann -f "Active site:blue,Propeptide:#000000"
```

**Example standalone**:
```
python -m clustal_to_svg -i ./alignment.clustal -o ./SVGs/ -u TRUE -c FALSE -n FALSE -a annotations.ann -f "Active site:blue,Propeptide:#000000"
```

**Usage**:
```
usage: clustal_to_svg.py [-h] [-i INFILE] [-o OUT_DIRECTORY] [-c CODES]
                         [-n NUMS] [-u UNIPROT_FORMAT] [-a ANNOTATIONS]
                         [-f FEATURES]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory.
  -c CODES, --codes CODES
                        When set to TRUE, includes Clustal identity codes at
                        the bottom of each block.
  -n NUMS, --nums NUMS  When set to TRUE, includes total residue numbers at
                        the end of each line.
  -u UNIPROT_FORMAT, --uniprot_format UNIPROT_FORMAT
                        When set to TRUE, truncates all accessions as if they
                        were UniProt entries. Ex. sp|P00784|PAPA1_CARPA ->
                        PAPA1_CARPA
  -a ANNOTATIONS, --annotations ANNOTATIONS
                        If an annotation file is provided, it will be used to
                        annotate the resulting SVG files.
  -f FEATURES, --features FEATURES
                        A comma-separated list of feature:color pairs to
                        include in SVGs. Case sensitive. If features include
                        spaces, the list must be enclosed in quotes. If no
                        features should be included, use: -f None The
                        following example is default behavior. Ex. -f "Active
                        site:#0000ff,Disulfide
                        bond:#e27441,Propeptide:#9e00f2,Signal:#2b7441"
```

## ANNOTATION FORMAT

Annotations can be added to an SVG with the -a or --annotations option in a run that calls clustal_to_svg.py.

A truncated, but real example of a valid annotation file can be found in annotation_example.ann.

**NOTE**: If no annotation is explicity provided, and retrieve_annotations.py is called during the run, clustal_to_svg.py will instead use the "all.ann" annotation file retrieved by retrieve_annotations.py. Annotation files provided via -a or --annotations will override those retrieved by retrieve_annotations.py.

**NOTE**: The -f or --features option was added to clustal_to_svg.py on 2024.10.30. It allows for users to provide a comma-separated list of feature:color pairs that can be used to customize SVG annotations.
Default behavior is identical to including the following option in a run that calls clustal_to_svg.py:
```
-f "Active site:#0000ff,Disulfide bond:#e27441,Propeptide:#9e00f2,Signal:#2b7441"
```

**Format**:

Annotation files that include the following columns and VALUES (tab-delimited) can be used as inputs for clustal_to_svg.py:

		prot	whole_prot	type	location.start.value	location.end.value
	ARBITRARY_INDEX	UNIPROT_FORMAT_ACC	FULL_ACCESSION	ANNOTATION_TYPE	START	END

 Other columns, like "description" may be added for record-keeping, but they will not be used when adding annotations to SVGs.

**Example**:

		prot	whole_prot	type	location.start.value	location.end.value
	0	PAPA1_CARPA	sp|P00784|PAPA1_CARPA	Active site	158	158

