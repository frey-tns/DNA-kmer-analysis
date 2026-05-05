"""
Compute sequence probabilities from a Markov background model.

SYNOPSIS USAGE
    Print usage line:
        seq-proba

    Print help message:
        seq-proba --help

    Usage:
        seq-proba [-h] --input FASTA_FILE
           --markov MARKOV_MATRIX_FILE --output OUTPUT_FILE

DESCRIPTION

This program computes the probability of each sequence (provided as a FASTA-formatted file) under a Markov background model.

The input consists of:
- a FASTA file containing DNA sequences
- a Markov transition matrix in RSAT tabular format
(typically produced by markov-from-seq or markov-from-kmers)

For each sequence, the program:
- reads the initial prefix probability P(prefix)
- multiplies it by the transition probabilities of each successive residue

If a sequence contains a prefix or context absent from the model, its probability it set to 0.

OPTIONS
    -h, --help
        Display this help message and exit.

    -i, --input FASTA_FILE
        Path to the input FASTA file.

    -m, --matrix MATRIX_MARKOV_FILE
        Path to the Markov transition matrix in tabular format.

    -o, --output OUTPUT_FILE
        Output TSV file.

OUTPUT

A tab-separated file with one row per sequence.

Column contents :
- id
    Sequence identifier from the FASTA file.
- length
    Sequence length of the FASTA file.
- proba_b
    Probability of the sequence under the Markov background model.
- log_p
    Log10 of the probability of the sequence under the Markov background model.

The results are written to a tab-separated value file (extension .tsv).

EXAMPLES

    seq-proba -i data/yeast_MET_upstream.fasta  -m data/markov_transitions_m2_2026_05_04.tsv -o ../result

AUTHOR / CREDITS
    Anouk RISCH
    supervised and revised by Jacques van Helden

VERSION
    0.1, 2026-05-04

CONTACT / URL
    https://github.com/frey-tns
    https://github.com/frey-tns/DNA-kmer-analysis

"""
version = 0.1
#################
#   Libraries   #
#################
import argparse
import time
import datetime
import os
import sys
import math

# Coloring warning text
from colorama import init, Fore

############################
#   Internal libraries     #
############################

# Read FASTA file
import bioseq_kmers.sequences as seq
import bioseq_kmers.utils as utils

################################################################
## FUNCTIONS
################################################################

###################################
#   Function: load markov model   #
###################################
def load_markov_matrix(path):
    """
    Read a markov transition matrix from a tabular RSAT-like file.

    The function reads a Markov background model produced by tolls such as `markov-from-seq` or
    `markov-from-kmers`. It parses a tab separated file containing prefix, transition probabilities
    ofr each nucleotide, and prefix probabilities.

    The prefix columns can be named either `pr\\su` or 'pr\\suf`.
    The Markov order is inferred from the length of the observed prefixes.

    Args:
        path (str): Path to te Markov transition matrix file (.tsv).
    Returns:
        dict: A dictionary containing:
        - dict_matrix (dict): Transition probabilities such taht:
        matrix[prefix][base] = P(base | prefix)
        - dict_prefix_prob (dict): Probability of each prefix P(prefix)
        - dict_order (int): Markov order (length of prefix)
    Raises:
        ValueError: If the file is malformed, missing required columns,
        or contains no valid matrix data.
    Notes:
        - Prefix are converted to uppercase internally
        - '.' is interpreted as an empty prefix (order 0 mode)
    """
    # Dict transition probabilities
    dict_matrix = {}
    # Dict where key is base and value is the initial probability per base
    dict_prefixes_prob = {}
    # column index
    col_index = None

    # Read file
    with open(path, 'r') as f:
        for line in f:
            # Remove line breaks and spaces
            line = line.strip()

            # Skip empty lines and RSAT comments
            if not line or line.startswith(";"):
                continue

            # Read header
            if line.startswith("#"):
                header = line.lstrip("#").strip().split()
                col_index = {name: i for i, name in enumerate(header)}
                # The header line should not be treated as a data line.
                continue

            # If header doesnt exist
            if col_index is None:
                # Error message
                raise ValueError("Missing header line starting with '#'.")

            # Columns are separated by tabs
            parts = line.split("\t")
            # If TSV malformed
            if len(parts) == 1:
                # Automatically separates on any space
                parts = line.split()

            # Ignore the comment lines
            if parts[0].startswith(";"):
                continue

            # Prefix column selection.
            if "pr\\su" in col_index:
                prefix_col = "pr\\su"
            elif "pr\\suf" in col_index:
                prefix_col = "pr\\suf"
            else:
                # Error message if prefix is missing
                raise ValueError("Missing prefix column ('pr\\suf' or 'pr\\su').")

            # Required columns
            for col in ["a", "c", "g", "t", "P_prefix"]:
                if col not in col_index:
                    # Error message if not required columns
                    raise ValueError(f"Missing column '{col}'.")

            # Read prefix
            prefix = parts[col_index[prefix_col]]
            if prefix == ".":
                prefix = ""
            else:
                prefix = prefix.upper()

            # Reading transition probabilities
            dict_matrix[prefix] = {"A": float(parts[col_index["a"]]),
                              "C": float(parts[col_index["c"]]),
                              "G": float(parts[col_index["g"]]),
                              "T": float(parts[col_index["t"]])}

            # Reading prefix probabilities
            dict_prefixes_prob[prefix] = float(parts[col_index["P_prefix"]])

    if not dict_matrix:
        raise ValueError("No markov transition matrix data found.")

    # Extract first prefix
    first_prefix = next(iter(dict_matrix))
    # Deducing the order of the model
    order = len(first_prefix)

    return {"matrix": dict_matrix,
            "prefixes_prob": dict_prefixes_prob,
            "order": order}

######################################
#   Function: Sequence probability   #
######################################
def sequence_probability(sequence, model):
    """
    Compute sequence probability of a biological sequence given a Markov background model.

    The probability is computed using an order k Markov chain, where k is defined by the model order.
    For order 0, bases are assumed independent. For higher orders, the probability is computed as:

        P(s) = P(prefix) × Π P(base_i | context_i)

    where context_i is the kmer preceding each base.

    Args:
        sequence (str): biological sequence (automatically converted to uppercase).
        model (dict): Markov model containing:
            - "matrix" (dict): Transition probabilities such that
            matrix[context][base] = P(base_i | context_i)
            - "prefixes_prob" (dict): Probability of each prefix P(prefix)
            -"order" (int): Markov order (length of prefix(kmer))
    Returns:
        float: Probability of the sequence under the Markov model.
        Returns 0.0 if the sequence is too short, if a required prefix
        is missing, or if an unknow context/base is encountered.

    Notes:
        - If order == 0, bases are assumed independent (Bernoulli).
        - All sequences are normalized to uppercase internally

    """

    # Sequence normalization
    sequence = sequence.upper()

    # Extract transition
    matrix = model["matrix"]
    # Extract probabilities
    prefixes_prob = model["prefixes_prob"]
    # Extract kmer length (=order)
    order = model["order"]

    ## Special case: order 0
    if order == 0:
        prob = 1.0
        for base in sequence:
            # Independently multiplies each base
            prob *= matrix[""].get(base.upper(), 0.0)
        return prob

    # If the sequence is too short
    if len(sequence) < order:
        return 0.0

    # Extract first prefix
    prefix = sequence[:order]

    # If no prefix
    if prefix not in prefixes_prob:
        return 0.0

    # Prefix initialize
    prob = prefixes_prob.get(prefix, None)
    if prob is None:
        return 0.0

    # Navigate the sequence starting from order
    for i_index in range(order, len(sequence)):
        # Extract prefix
        context = sequence[i_index - order:i_index]
        # Next base
        base = sequence[i_index]

        # If no prefix
        if context not in matrix:
            return 0.0

        # Markov multiplication P(base∣context)
        prob *= matrix[context].get(base, 0.0)

    return prob

#################
#   Main code   #
#################
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
    parser = argparse.ArgumentParser(description="Sequence probability from background model")

    # Define args used by the user (input, output, Markov order)
    parser.add_argument("-i", "--input",
                        required=True,
                        help="input FASTA file")

    parser.add_argument("-m", "--matrix",
                        required=True,
                        help="Markov matrix TSV file")

    parser.add_argument("-o", "--output",
                        required=False,
                        help="output TSV file")

    # Reads the command typed in the terminal
    args = parser.parse_args()

    ### Define variable to use the value in the script
    # Load sequence
    input_file = args.input
    model_file = args.matrix
    output_file = args.output

    if output_file is None:
        output_file = "."

    ############################
    #   Load data
    ############################

    # sequences
    sequences, _, _ = seq.read_fasta(input_file)

    # Markov model
    model = load_markov_matrix(model_file)

    #####################
    #   Output file     #
    #####################

    # Current date
    today = str(datetime.date.today()).replace("-", "_")

    # If the output path is a folder
    if os.path.isdir(output_file):
        # Define output path
        output_path = os.path.join(output_file, f"seq_proba_{today}.tsv")
    else:
        # Force the TSV extension
        if not output_file.endswith(".tsv"):
            output_file += ".tsv"
        # If it's a file
        output_path = output_file

    # Write TSV file output
    with open(output_path, "w") as tsv_file:
        ## Parameter
        # Command line
        command_line = utils.format_command_line(sys.argv)
        # Write command line
        tsv_file.write(f"; seq-proba\t{command_line}\n;\n")
        # URL in input
        tsv_file.write(f"; Program version\t{version}\n"
                       f"; Input file\t{os.path.relpath(input_file)}\n")
        # Header
        tsv_file.write("#id\tlength\tproba_b\tlog_proba\n")

        for seq_id, sequence in sequences.items():
            p = sequence_probability(sequence, model)
            length_seq = len(sequence)

            # Log probabilities
            if p == 0:
                log_p = float("-inf")
            else:
                log_p = math.log10(p)

            tsv_file.write(f"{seq_id}\t{length_seq}\t{p:.6e}\t{log_p:.2f}\n")


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