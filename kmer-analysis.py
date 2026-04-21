__authors__ = ("Anouk RISCH")
__contact__ = ("anouk.risch@etu.univ-amu.fr")
__date__ = "2026-04-21"
__version__ = "1.2"

#################
#   Libraries   #
#################
# To parser for command-line arguments
import argparse
# Interaction with the operating system for file management and manipulation
import os
# For Benchmark
import time
import datetime
import sys
import shlex

from tqdm import tqdm
from collections import Counter # count (cf k-mer count)
# Coloring warning text
from colorama import init, Fore
init(autoreset=True)

#####################################################
#    Function to read and give stat on Fasta file   #
#####################################################
# Input: Fasta file from RSAT website (DEMO)        #
#                                                   #
# Outputs: - Number of sequence in the fasta file   #
#          - Sum of length                          #
#####################################################

## READ FASTA FILE
def read_fasta(file_path):
    """
    Load and read FASTA file from local path.

    Args:
        path file (str): FASTA file.

    Returns:
        dict: dictionary {ID : sequence}
    """

    # Did file path already exist
    if not os.path.exists(file_path):
        # If not stop the program
        raise FileNotFoundError(f"{Fore.RED}[ERROR] Input FASTA file not found: {file_path}\n"
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
    print(f"{Fore.CYAN}{len(dico_sequences)} loaded sequences")

    # Sum of lengths
    total_length = sum(len(seq) for seq in dico_sequences.values())
    print(f"{Fore.CYAN}Sum of lengths : {total_length}")

    return dico_sequences, total_length

####################################################
#   Function: complementary reverse of the k-mer   #
####################################################
def reverse_complementary(kmer_seq):
    """
    Reverse complement of a DNA k-mer sequence.

    The input sequence is reversed, then each nucleotide is replaced by its complementary base:
    A <-> T
    C <-> G

    Args:
        kmer_seq (str): DNA k-mer sequence (A, C, G, T).
    Returns:
        str: Reverse complement of DNA k-mer sequence.
    Raises:
        KeyError: If the sequence contains unsupported characters.
    """

    # Create complementary reverse dictionary
    dico_complement = {"A":"T",
                       "T":"A",
                       "C":"G",
                       "G":"C"}

    # Trace the Kmer backwards and replace each base with its complement
    return "".join(dico_complement[base] for base in reversed(kmer_seq))

#####################################
#      Function: Canonical Kmer     #
#  Merge kmer + reverse complement  #
#####################################
def canonic_kmer(kmer):
    """
    Return  the canonical representation of a DNA k-mer.

    The canonical k-mer is defined as the lexicagraphically smallest sequence between the input k-mer and
    its reverse complement sequence. This allows merging canonical occurrences from both DNA strand.

    Args:
        kmer (str): DNA k-mer sequence (A, C, G, T).
    Returns:
        str: Canonical k-mer sequence.
    Raises:
        KeyError: If the sequence contains unsupported characters.
    Examples:
        >>> canonic_kmer("ATGC")
        'ATCG'
        >>> canonic_kmer("GCAT")
        'ATCG'

    """
    # Kmer reverse complement
    reverse_kmer = reverse_complementary(kmer)

    # Keep the smallest alphabetically
    return min(kmer, reverse_kmer)

#########################################
#   Function: nucleotides frequencies   #
#########################################
def nucleotide_frequencies(sequences):
    """
    Calculate nucleotide frequencies from a list of DNA k-mer sequences.

    The function counts occurrences of valid nucleotides (A, C, G, T)
    across all DNA sequences and ignores unknown bases such as N.

    Args:
        sequences (dict): Dictionary of DNA sequences {sequence_id : sequence}.
    Returns:
        dict: Relative frequencies of nucleotide  {"A" : freq_A,
                                                   "C" : freq_C,
                                                   "G" : freq_G,
                                                   "T" : freq_T,}
    Raises:
         ZeroDivisionError: If no valid nucleotide is found.
    """
    # Initialization total number of nucleotides
    total = 0
    # Initialization a dictionary counter where key is base and value the number of times it was counted
    # {"A":0, "T":0, "G":0, "C":0}
    counts = Counter()

    for seq in sequences.values():
        for base in seq:
            # Ignore N nucleotide
            if base in "ACTG":
                # Count each letter (+1 for each letter)
                counts[base] += 1
                # Total number of nucleotides
                total += 1

    # Create dico contain expected frequency for each base
    frequencies = {base: counts[base] / total for base in "ATGC"}

    return frequencies

#########################################
#   Function : Excepted frequencies     #
#########################################
def expected_frequencies(single_kmer, frequencies):
    """
    Calculate expected frequencies from a list of DNA k-mer sequences.

    The function counts occurrences of valid nucleotide (A, C, G, T).
    The excepted frequency is computed as the product of the individual nucleotide frequencies,
    assuming independence between positions.

    For example:
        P(ATG) = P(A) × P(T) × P(G)

    Args:
        single_kmer (str): DNA k-mer sequence (A, C, G, T).
        frequencies (dict): Dictionary of nucleotide frequencies {"A" : freq_A,
                                                                  "C" : freq_C,
                                                                  "G" : freq_G,
                                                                  "T" : freq_T,}.
    Returns:
        float: Expected frequency of the k-mer.
    Raises:
        KeyError: If the k-mer contains unsupported characters.
     Example:
        >>> freqs = {"A": 0.3, "T": 0.2, "G": 0.3, "C": 0.2}
        >>> expected_frequencies("ATG", freqs)
        0.018
    """
    probability = 1

    for base in single_kmer:
        # Expected number of occurrences for each base
        probability *= frequencies[base]

    return probability

#####################################
#   Function : Count observed kmer  #
#####################################
def counts_kmer(k_length, sequence, strand_mode):
    """
    Count observed k-mers in a set of DNA sequences.

    The function scans each sequence using a sliding window of size k and counts all valid k-mers.

    K-mers containing unknown nucleotides (N) are ignored.

    Depending on the strand mode:
        - "single": count k-mers as they appear.
        - "both": merge each k-mer with its reverse complement using the canonical representation.

    Args:
        k_length (int): Length of the k-mer.
        sequence (dict): Dictionary of DNA sequences in the form
            {sequence_id: sequence}.
        strand_mode (str): Strand counting mode:
            "single" or "both".
    Returns:
        Counter(dict): Dictionary-like object containing k-mer counts {kmer: occurrences}
    """

    # Initialize a dictionary counter where key is kmer and value is the number of times it was counted
    kmer_count = Counter()

    # Progress bar
    for seq_fasta in tqdm(sequence.values(), desc="Counting k-mers"):
       # Explore all possible positions
       for i_position in range(len(seq_fasta) - k_length + 1):
           # Extract kmer
           kmer = seq_fasta[i_position: i_position + k_length]

           # Exclude unknow nucleotides
           if "N" in kmer:
               continue

           if strand_mode == "both":
               key = canonic_kmer(kmer)

           else:
               # single strand (forward)
               key = kmer

           # Kmer found at this position, after reverse complement fusion
           # +1 every time we see a pattern
           kmer_count[key] += 1

    return kmer_count

###################
#   Main code     #
###################
def main():

    # Time tracking (Benchmark)
    start_time = time.perf_counter()
    # Job started
    start_time_date = datetime.datetime.now()
    ############################
    #   Command line options   #
    ############################

    ## OUTPUT DIRECTORY FILE

    # Specify which command-line options the program is willing to accept
    parser = argparse.ArgumentParser(description="k-mer analysis")
    # Define args used by the user (here output path)
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Path to the input fasta file")

    parser.add_argument("-k", "--kmer-length",
                        required=True,
                        type=int,
                        help="Length of k-mer sequence (1-10)")

    parser.add_argument("-o", "--output",
                        required=True,
                        help="Path to the output file directory")

    parser.add_argument("-s", "--strand",
                        choices=["single", "both"],
                        default="single",
                        help="Sequence strand(s) to take into account for k-mer counting. Accepted values: single, only count the occurrences on forward strand; both: count on both forward and reverse strands.)")

    # Reads the command typed in the terminal
    args = parser.parse_args()

    # Define variable to use the value in the script
    input_file = args.input
    kmer_length = args.kmer_length
    output_file = args.output
    strand_mode = args.strand

    ## CONDITION : does the files already exist ?

    # Extract the folder from the full path
    folder = os.path.dirname(output_file)

    # If the folder is not empty and doesn't already exist :
    if folder and not os.path.exists(folder):
        # Warning message
        print(f"{Fore.YELLOW}Creating folder {folder}")
        # Create the folder
        os.makedirs(folder)

    #############################################
    #   Defined the address of the fasta file   #
    #############################################
    # Input URL of the FASTA file
    fasta_file = input_file

    sequences, total_length = read_fasta(fasta_file)

    ###################
    #    Statistics   #
    ###################

    # Nucleotides frequency
    frequencies = nucleotide_frequencies(sequences)

    observed_kmer_count = counts_kmer(kmer_length, sequences, strand_mode)
    # Number of all positions T = L - K + 1
    total_positions = sum(len(seq) - kmer_length + 1 for seq in sequences.values())

    # List who contain output results
    result_analysis = []

    # Treat only once the same canonical kmer
    dico_canon_kmer = set()

    # Browses the dictionary of observed k-mers.
    for canon_kmer, occ in observed_kmer_count.items():

        if strand_mode == "reverse_complement":
            canon_kmer = canonic_kmer(canon_kmer)
            # Kmer inverse complement
            reverse_kmer = reverse_complementary(canon_kmer)
            # Define ID
            kmer_id = f"{canon_kmer}|{reverse_kmer}"

        else:
            canon_kmer = canon_kmer
            # Define ID
            kmer_id = canon_kmer

        # If the canon kmer has already been dealt with
        if canon_kmer in dico_canon_kmer:
            # Ignore it
            continue

        dico_canon_kmer.add(canon_kmer)

        threshold = 1
        if occ < threshold:
            continue

        # Expected frequencies
        exp_freq = expected_frequencies(canon_kmer, frequencies)
        # Expected occurrences
        exp_occ = total_positions * exp_freq
        # Observed frequencies
        obs_freq = occ / total_positions

        # Stockage in list oligomers
        result_analysis.append({"seq": canon_kmer,
                                "id": kmer_id,
                                "exp_freq": exp_freq,
                                "exp_occ": exp_occ,
                                "occ": occ,
                                "obs_freq": obs_freq})

    #####################
    #   Output file     #
    #####################

    # Current date
    today = str(datetime.date.today()).replace("-", "_")
    # # Create DataFrame
    # df = pd.DataFrame(result_analysis, columns=["seq", "id", "exp_freq", "exp_occ"])
    # # Convert DataFrame into HTML format
    # df_HTML = df.to_html(escape=False, index=False)

    # If the output path is a folder
    if os.path.isdir(output_file):
        # Define output path
        output_path = os.path.join(output_file,f"summary_kmer_analysis_{today}.tsv")
    else:
        # Force the HTML extension
        if not output_file.endswith(".tsv"):
            output_file += ".tsv"
        # If it's a file
        output_path = output_file
    # Write HTML file output
    with open(output_path, "w") as tsv_file:

        ## Parameter
        # Command line
        tsv_file.write(f"; Command\t{' '.join(shlex.quote(arg) for arg in sys.argv)}\n;\n")
        # URL in input
        tsv_file.write(f"; Fasta URL\t{fasta_file}\n"
                       f"; Input file\t{input_file}\n"
                       f"; Input format\tFasta\n"
                       f"; Output file\t{output_path}\n"
                       f"; Strand mode\t{strand_mode}\n"
                       f"; Oligomer length\t{kmer_length}\n;\n")

        # Summary
        tsv_file.write(f"; Sequence type\tDNA\n"
                       f"; Nb of sequence\t{len(sequences)}\n"
                       f"; Sum of sequence lengths\t{total_length}\n"
                       f"; Sequences:\n")
        for seq_id, seq_fasta in sequences.items():
            tsv_file.write(f"#\t{seq_id}\t{len(seq_fasta)}\n")

        tsv_file.write(f";\n;\n")

        # Column headers
        tsv_file.write(f"; column headers\n")
        tsv_file.write(f"; seq\t oligomer sequence\n"
                       f"; id\t oligomer identifier\n"
                       f"; exp_freq\t expected relative frequency\n"
                       f"; obs_freq\t observed relative frequency\n"
                       f"; occ\t observed occurrences\n"
                       f"; exp_occ\t expected occurrences\n;\n")

        # Kmer analysis table headers
        tsv_file.write(f"# seq\tid\texp_freq\tobs_freq\tocc\texp_occ\n")
        for row in result_analysis:
            # kmer analysis table
            tsv_file.write(f"{row['seq']}\t{row['id']}\t{row['exp_freq']}\t{row['obs_freq']}\t{row['occ']}\t{row['exp_occ']:.2f}\n")

        # Double return to the line
        tsv_file.write(f";\n;\n")

        # End time
        end_time = time.perf_counter()
        # Job ending
        end_time_date = datetime.datetime.now()
        duration = end_time - start_time

        tsv_file.write(f"; Job started\t{start_time_date}\n"
                       f"; Job done\t{end_time_date}\n"
                       f"; Job duration\t{duration:.3f} seconds\n")

    print(f"{Fore.GREEN}Output written to {output_path}")
    print(f"{Fore.CYAN}Duration : {duration:.3f} seconds\n")

#####################
#   Executing code  #
#####################
if __name__ == "__main__":
    main()


