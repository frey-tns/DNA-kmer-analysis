#################
#   Libraries   #
#################
# For HTTP requests
import requests

# For generate all possible combinations with repetition
from itertools import product
from collections import Counter # count (cf k-mer count)

#############
#   Paths   #
#############

## Requested enhancement: https://github.com/frey-tns/DNA-kmer-analysis/issues/2

# Output directory
path_output = r"C:\Users\Anouk\Documents\Master1_bsg\Stage"

# Output file
outputs_file = (f"{path_output}\\kmer-analisys.txt")

#####################################################
#    Function to read and give stat on Fasta file   #
#####################################################
# Input: Fasta file from RSAT website (DEMO)        #
#                                                   #
# Outputs: - Number of sequence in the fasta file   #
#          - Sum of lengths                         #
#####################################################

## READ FASTA FILE
def Read_fasta(url):
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
    print(f"Sum of lenghts : {total_length}")

    return dico_sequences

#################################################################
#   Function to obtain the complementary reverse of the k-mer   #
#################################################################
def Reverse_complemantary(kmer_seq):
    # Creat complementary reverse dictionary
    complement = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
    return "".join(complement[base] for base in reversed(kmer_seq))

#########################################
#   Function: nucleotides frequencies   #
#########################################
def Nucleotide_frequencies(sequences):
    # Initialization total number of nucleotides
    total = 0
    # Initialization counter (A, T, G, C)
    counts = Counter()

    for seq in sequences.values():
        # Count each letter (+1 for each letter)
        counts.update(seq)
        # Adds the length of the sequence
        total += len(seq)

    # Calculate expected frequency
    frequencies = {nuc: counts[nuc] / total for nuc in "ATGC"}
    return frequencies

#########################################
#   Function : Excepted occurrences     #
#########################################

## Excepted occurence
def Exp_occ(single_kmer, sequences, frequencies):
    # Number of all positions T = L - K + 1
    # Sum of length of sequences
    L = sum(len(seq) for seq in sequences.values())
    # Length K-mer
    K = len(single_kmer)
    # Possible positions for a k-letter word in the sequence set
    T = L - K + 1

    ## Expected frequencies
    # Probability of k-mer under background nucleotide frenquencies
    probability = 1
    for base in single_kmer:
        # Expected number of occurrences for word
        probability *= frequencies[base]

    print(f"Expected frequencies {single_kmer} : {probability}")
    return T*probability

#############################################
#   Defined the address of the fasta file   #
#############################################

url = input("URL of the sequence : ")

# Reference URL
# "https://rsat.eead.csic.es/plants/tmp/www-data/2026/04/08/tmp_sequence_2026-04-08.113931_zOkxjo.fasta"
sequences = Read_fasta(url)

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

# Generates all combinaisons to create pattern
k_mer = ["".join(p) for p in product(list_nucl, repeat = k_length)]

print(k_mer)
# Length k-mer
print("Total number of k-mer :", len(k_mer))

#####################
#    Statistiques    #
#####################

# Nucleootides frequency
frequencies = Nucleotide_frequencies(sequences)
print(f"Nucleotides frequencies : {frequencies}")

# Expected occurrences
for single_kmer in k_mer:
    exp_occ = Exp_occ(single_kmer, sequences, frequencies)
    print(f"Exepcted occurrence {single_kmer} : {exp_occ}")

