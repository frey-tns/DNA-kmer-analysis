"""
Read a FASTA file and compute basic statistics on biological sequences.

This module provides a function to load sequences from a FASTA-formatted file
and compute simple descriptive statistics, including:

- the number of sequences
- the total length of all sequences

The sequences are stored in a dictionary mapping sequence IDs to their
corresponding nucleotide sequences.

The function also validates the input file path and displays progress
during parsing.

USAGE
    This module is intended to be imported and used in other scripts :

    from sequences import read_fasta
    sequences, total_length = read_fasta("input.fasta")

INPUT
    file_path : str
    Path to the input FASTA file.

AUTHOR
    Anouk RISCH

CONTACT
    https://github.com/frey-tns

URL
    https://github.com/frey-tns/DNA-kmer-analysis

VERSION
    1.2, 2026-04-24

"""

#################
#   Libraries   #
#################
import os

# progress bar
from tqdm import tqdm

#####################################################
#    Function to read and give stat on Fasta file   #
#####################################################
# Input: Fasta file                                 #
#                                                   #
# Outputs: - Number of sequence in the fasta file   #
#          - Sum of length                          #
#####################################################

## READ FASTA FILE
def read_fasta(file_path):
    """
    Load and read FASTA file from local path.

    Args:
        file_path (str): FASTA file.

    Returns:
        dict: dictionary {ID : sequence}
    """

    # Did file path already exist
    if not os.path.exists(file_path):
        # If not stop the program
        raise FileNotFoundError(f"[ERROR] Input FASTA file not found: {file_path}\n"
                                f"Check working directory or use absolute path.")

    # Initialize dictionary who contain FASTA sequences
    dico_sequences = {}
    # Initialize FASTA ID as an empty string
    current_id = ""

    with open(file_path, "r") as fasta_file:
        # For each line present in FASTA file
        for line in tqdm(fasta_file, desc="Reading FASTA file"):
            # Remove the spaces on the right
            line = line.rstrip()

            # FASTA pipe detection
            if line.startswith(">"):
                # Extract ID
                current_id = line[1:].split()[0]
                # Adds FASTA ID as key in sequences dictionary
                dico_sequences[current_id] = ""

            else:
                # If there is a key in sequences dictionary
                if current_id != "":
                    # Adds DNA sequence as value in sequences dictionary
                    dico_sequences[current_id] += line.upper()

    ## STAT

    # Number of sequence
    seq_number = len(dico_sequences)

    # Sum of lengths
    total_length = sum(len(seq) for seq in dico_sequences.values())

    return dico_sequences, total_length, seq_number
