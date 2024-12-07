"""
Classes for storing protein alignment data.

Classes:
    Protein
    Alignment
    Feature
"""

import math

class Protein:
    """
    Stores sequences, names, and features of a protein.

    Attributes:
        name (str): Original name of the protein.
        disp_name (str): Name that should be displayed in Martin Lab alignment.
        sequence (str): Aligned protein sequence (includes "-" for gaps).
        seq_raw (str): Raw protein sequence (no "-").
        features (list): List of UniProt-derived Feature objects for this protein.
    """

    def __init__(self, name, sequence):
        """
        Initiates Protein.

        Parameters:
            name (str): Original name of the protein.
            sequence (str): Aligned protein sequence (includes "-" for gaps).
        """

        self.name = name
        self.disp_name = name # Unless set otherwise.
        # In aligned FASTA format, with "-" for gaps (no >name, just sequence).
        self.sequence = sequence
        self.seq_raw = self.sequence.replace("-", "")
        self.features = []

    def add_feature(self, feature):
        """
        Sets a given feature's Clustal positions and adds it to self.features.

        Parameters:
            feature (Feature): Feature object to be added to this protein.
        """

        clust_start = None # 1-indexed
        clust_end = None # 1-indexed

        raw_count = 0
        clust_count = 0
        # Note: will run through whole sequence if feature.start == feature.end.
        for position in self.sequence:
            clust_count += 1
            # Test this first; cannot end on a "-"
            if position == "-":
                continue
            else:
                raw_count += 1                
            if raw_count == feature.start:
                clust_start = clust_count
                if feature.start == feature.end:
                    clust_end = clust_count
                    break # End the loop if start and end found.
            elif raw_count == feature.end:
                clust_end = clust_count
                break # End the loop once it finds the end of the feature.

        # clust_start and clust_end are 1-indexed.
        if feature.start and feature.end: 
            feature.set_seq(self.seq_raw[feature.start-1:feature.end])
        feature.set_clust_positions(clust_start, clust_end)
        self.features.append(feature)

    def set_disp_name(self, length=1000, uniprot_format="FALSE", query="FALSE"):
        """
        Sets a shorter display name for this protein.

        Parameters:
            length (int): Max length of the final display name.
                          This will be overridden by the true
                          max header of an Alignment object.
            uniprot_format (str): If "TRUE", uses a UniProt-specific truncation.
                                  Ex. sp|P00784|PAPA1_CARPA -> PAPA1_CARPA
        """

        if query == "TRUE":
            self.disp_name = self.disp_name.replace("QUERY_", "")[:length]
        elif uniprot_format == "TRUE":
            self.disp_name = self.disp_name.split("|")[-1][:length]
        else:
            self.disp_name = self.disp_name[:length]

class Alignment:
    """
    Stores Clustal Omega alignment information.

    Attributes:
        proteins (dict): Dictionary of Protein objects in the alignment.
                         Stored as Protein.name:Protein key:value pairs.
        codes (str): String of Clustal Omega conservation codes (*:.).
        max_header (None or int): Maximum header length if one is set.
                                  Header length affects how Protein.disp_name
                                  is displayed in the final alignment.
        conserved_res (list): List of Clustal positions for conserved residues.
    """

    def __init__(self):
        """Initiates Alignment."""

        # Dictionary of Proteins; keys are Protein.name; consider using just a list...
        self.proteins = {}
        self.codes = ""
        self.max_header = None
        self.conserved_res = [] # List of clust_positions for conserved residues.

    def set_conserved_res(self):
        """Sets the Clustal positions for conserved residues in this Alignment."""

        # Sets conserved residues only on an external call.
        i = 1
        prot_list = list(self.proteins.values())
        # Break out when sequences have been entirely gone through.
        while i <= len(prot_list[0].sequence):
            conserved = True
            res_0 = prot_list[0].sequence[i-1]
            for protein in prot_list[1:]:
                res_new = protein.sequence[i-1]
                if res_0 != res_new or res_0 == "-" or res_new == "-":
                    conserved = False
                    break # Breaks out of this inner for loop.

            if conserved:
                self.conserved_res.append((res_0, i))
            i += 1

    def add_protein(self, protein):
        """
        Adds the given Protein object to self.proteins (dict).

        Parameters:
            protein (Protein): Protein object to be stored.
        """

        self.proteins[protein.name] = protein

    def add_codes(self, codes):
        """
        Stores the Clustal conservation codes in this Alignment.

        Parameters:
            codes (str): A string of Clustal conservation codes (*:.).
        """

        self.codes = codes

    def set_max_header(self, max_header):
        """
        Sets the max header length for this Alignment.

        Parameters:
            max_header (int): Maximum length for protein names
                              in Martin Lab SVG format.
        """

        self.max_header = max_header

    def get_max_length(self):
        """
        Returns the maximum length of all proteins in this Alignment.

        Returns:
            length (int): Length of the longest raw protein
                          sequence in this Alignment.
        """

        length = 0
        for name in self.proteins:
            length = max(len(self.proteins[name].seq_raw), length)
        return length

    def add_features(self, annotations):
        """
        Creates and adds Feature objects to Protein objects in this Alignment.

        Parameters:
            annotations (pandas.DataFrame): DataFrame of UniProt
                                            annotations for proteins
                                            in this Alignment.
        """

        for row in annotations.itertuples(index=False):

            feature_type = row[annotations.columns.get_loc("type")]
            start = row[annotations.columns.get_loc("location.start.value")]
            end = row[annotations.columns.get_loc("location.end.value")]

            # If no start or end, set those values to None.
            # If they exist, make them integers.
            if math.isnan(float(start)):
                start = None
            else:
                start = int(start)
            if math.isnan(float(end)):
                end = None
            else:
                end = int(end)

            try:
                feature = Feature(feature_type, start, end, row.description)
            except AttributeError:
                feature = Feature(feature_type, start, end, "")

            prot_temp = self.proteins.get(row.whole_prot)
            if prot_temp:
                self.proteins[row.whole_prot].add_feature(feature)

        # TODO: Consider adding a function that returns feature coordinates here.

    def set_disp_names(self, length=1000, uniprot_format="FALSE"):
        """
        Sets display names for Protein objects in this Alignment.

        Parameters:
            length (int): Max length of the final display name.
                          This will be overridden by the true
                          max header of an Alignment object.
            uniprot_format (str): If "TRUE", uses a UniProt-specific truncation.
                                  Ex. sp|P00784|PAPA1_CARPA -> PAPA1_CARPA
        """

        # Length set to 1000 just so full length is kept.
        for name in self.proteins:
            self.proteins[name].set_disp_name(length=length, uniprot_format=uniprot_format)

class Feature:
    """
    Stores information for Protein object annotations (features).

    Attributes:
        feature_type (str): Type of annotation. Ex. "Active Site".
        start (int): Start residue number in the raw protein sequence.
        end (int): End residue number in the raw protein sequence.
        clust_start (None or int): Start Clustal position for this Feature.
        clust_end (None or int): End Clustal position for this Feature.
        description (str): Description of this Feature.
        
    """

    # Deal with mismatched sequences later (maybe use Carter's 5res before/after).
    # Start and end are 1-indexed, like a normal protein sequence.
    def __init__(self, feature_type, start, end, description):
        """
        Initiates Feature.

        Parameters:
            feature_type (str): Type of annotation. Ex. "Active Site".
            start (int): Start residue number in the raw protein sequence.
            end (int): End residue number in the raw protein sequence.
            description (str): Description of this Feature.
        """

        self.feature_type = feature_type
        self.start = start # Position in raw sequence.
        self.end = end # Position in raw sequence.
        self.clust_start = None # Set when added to Protein.
        self.clust_end = None # Set when added to Protein.
        self.description = description
        self.sequence = None # Set when added to Protein.
        # TODO: Consider adding self.sequence in order to identify later...

    # Clustal positions based on start and end.
    # Called from Protein object.
    def set_clust_positions(self, clust_start, clust_end):
        """
        Sets Clustal positions for this Feature.

        Parameters:
            clust_start (int): Start Clustal position for this Feature.
            clust_end (int): End Clustal position for this Feature.
        """
        
        #print("Set Clustal positions " + str(clust_start) +
              #" " + str(clust_end), flush=True)
        self.clust_start = clust_start
        self.clust_end = clust_end

    def set_seq(self, sequence):
        """
        Sets sequence for Feature.

        Parameters:
            sequence (str): Amino acid sequence (no gaps) of this Feature.
        """

        #print("Set sequence " + str(sequence), flush=True)
        self.sequence = sequence