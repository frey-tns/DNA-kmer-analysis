__authors__ = ("Anouk RISCH")
__contact__ = ("anouk.risch@etu.univ-amu.fr")
__date__ = "2026-04-13"
__version__ = "1.2"

#################
#   Libraries   #
#################
# For HTTP requests
import requests
# To parser for command-line arguments
import argparse
# Interaction with the operating system for file management and manipulation
import os
import datetime
# For HTML result
import pandas as pd

# For generate all possible combinations with repetition
from itertools import product
from collections import Counter # count (cf k-mer count)

#####################################################
#    Function to read and give stat on Fasta file   #
#####################################################
# Input: Fasta file from RSAT website (DEMO)        #
#                                                   #
# Outputs: - Number of sequence in the fasta file   #
#          - Sum of length                          #
#####################################################

## READ FASTA FILE
def read_fasta(url):
    """
    Load and read FASTA file from URL.

    Args:
        url (str): FASTA file URL.

    Returns:
        dict: dictionary {ID : sequence}
    """

    # HTTP request
    response = requests.get(url)
    # Status code HTTP
    if response.status_code != 200:
        raise Exception(f"FASTA download error ({response.status_code})")
    # Create variable who contain the FASTA file
    fasta_file = response.text

    # Initialize dictionary who contain FASTA sequences
    dico_sequences = {}
    # Initialize FASTA ID as an empty string
    current_id = ""

    # For each line present in FASTA file
    for line in fasta_file.splitlines():
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
                dico_sequences[current_id] += line

    ## STAT

    # Number of sequence
    print(f"{len(dico_sequences)} loaded sequences")

    # Sum of lengths
    total_length = sum(len(seq) for seq in dico_sequences.values())
    print(f"Sum of lengths : {total_length}")

    return dico_sequences

#################################################################
#   Function to obtain the complementary reverse of the k-mer   #
#################################################################
def reverse_complementary(kmer_seq):
    # Create complementary reverse dictionary
    dico_complement = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join(dico_complement[base] for base in reversed(kmer_seq))

#########################################
#   Function: nucleotides frequencies   #
#########################################
def nucleotide_frequencies(sequences):
    # Initialization total number of nucleotides
    total = 0
    # Initialization counter (A, T, G, C)
    counts = Counter()

    for seq in sequences.values():
        # Count each letter (+1 for each letter)
        counts.update(seq)
        # Adds the length of the sequence
        total += len(seq)

    # Create dico contain expected frequency for each base
    frequencies = {nuc: counts[nuc] / total for nuc in "ATGC"}
    return frequencies

#########################################
#   Function : Excepted frequencies     #
#########################################
def expected_frequencies(single_kmer, frequencies):
    probability = 1
    for base in single_kmer:
        # Expected number of occurrences for word
        probability *= frequencies[base]

    return probability

###################
#   Main code     #
###################
def main():
    #############
    #   Paths   #
    #############

    ## OUTPUT DIRECTORY FILE

    # Specify which command-line options the program is willing to accept
    parser = argparse.ArgumentParser()
    # Define args used by the user (here output path)
    parser.add_argument("-o", "--output", required=True, help="Path to the output file directory")

    # Reads the command typed in the terminal
    args = parser.parse_args()
    # Define variable to use the value in the script
    output_file = args.output

    ## CONDITION : does the files already exist ?
    # Extract the folder from the full path
    folder = os.path.dirname(output_file)

    # If the folder is not empty and doesn't already exist :
    if folder and not os.path.exists(folder):
        # Warning message
        print("Creating folder {}".format(folder))
        # Create the folder
        os.makedirs(folder)

    #############################################
    #   Defined the address of the fasta file   #
    #############################################

    # Input URL of the FASTA file
    url = input("URL of the sequence : ")

    # REFERENCE URL
    # https://rsat.eead.csic.es/plants/tmp/www-data/2026/04/08/tmp_sequence_2026-04-08.113931_zOkxjo.fasta
    sequences = read_fasta(url)

    ## Access to a specific sequence
    #print(sequences.keys())
    #print(sequences["MET8"])

    ######################
    #   Create k-mer     #
    ######################

    # List of nucleotide
    list_nucl = ['A','C','G','T']
    # Length k-mer in base pair
    k_length = int(input("Enter k-mer length (1-10): "))

    # Generates all combinations to create pattern
    k_mer = ["".join(p) for p in product(list_nucl, repeat = k_length)]

    # print(k_mer)
    # # Length k-mer
    # print("Total number of k-mer :", len(k_mer))

    ###################
    #    Statistics   #
    ###################

    # Nucleotides frequency
    frequencies = nucleotide_frequencies(sequences)
    #print(f"Nucleotides frequencies : {frequencies}")

    # Number of all positions T = L - K + 1
    total_positions = sum(len(seq) - k_length + 1 for seq in sequences.values())

    # List who contain output results
    result_analysis = []

    for single_kmer in k_mer:
        # Expected frequencies
        exp_freq = expected_frequencies(single_kmer, frequencies)
        # Expected occurrences
        exp_occ = total_positions * exp_freq
        # Reverse complementary
        rev_comp = reverse_complementary(single_kmer)

        # Oligomer ID = k-mer + reverse complement
        oligomer_id = f"{single_kmer}|{rev_comp}"

        # Stockage in list oligomers
        result_analysis.append({"seq": single_kmer,
                                "id": oligomer_id,
                                "exp_freq": exp_freq,
                                "exp_occ": exp_occ})

    # Current date
    today = str(datetime.date.today()).replace("-", "_")
    # Create DataFrame
    df = pd.DataFrame(result_analysis, columns=["seq", "id", "exp_freq", "exp_occ"])
    # Convert DataFrame into HTML format
    df_HTML = df.to_html(escape=False, index=False)

    # If the output path is a folder
    if os.path.isdir(output_file):
        # Define output path
        output_path = os.path.join(output_file,f"summary_kmer_analysis_{today}.html")
    else:
        # Force the HTML extension
        if not output_file.endswith(".html"):
            output_file += ".html"
        # If it's a file
        output_path = output_file
    # Write HTML file output
    with open(output_path, "w") as html_file:
        html_file.write(df_HTML)

    print(f"Output written to {output_path}")

#####################
#   Executing code  #
#####################
if __name__ == "__main__":
    main()

