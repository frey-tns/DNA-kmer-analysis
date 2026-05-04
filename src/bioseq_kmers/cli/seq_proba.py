"""
Estimate the parameters of a Markov background model from a kmer counting table and output a transition matrix.

SYNOPSIS USAGE
    Print usage line:
        seq-proba

    Print help message:
        seq-proba --help

    Usage:
        seq-proba [-h] --input FASTA_FILE
           --markov MARKOV_MATRIX_FILE --output OUTPUT_FILE

DESCRIPTION


OPTIONS
    -h, --help
        Display this help message and exit.

    -i, --input FASTA_FILE
        Path to the input FASTA file.

    -m, --matrix MATRIX_MARKOV_FILE
        Input Markov matrix.

    -o, --output OUTPUT_FILE
        Output TSV file.

OUTPUT

A tab-separated file with one row per sequence.

Column contents :
- seq ID (sequence ID)
- seq proba (sequence probability)

The results are written to a tab-separated value file (extension .tsv).

EXAMPLES

    seq-proba -i data/3nt_counts.tsv -m markov_transitions_m2_2026_05_04.tsv -o results/

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
    Read a markov transition matrix

    Args:

    Returns:
    """
    matrix = {}
    prefixes_prob = {}

    with open(path, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith(";") or line.startswith('#'):
                continue

            parts = line.split("\t")

            prefix = parts[0]
            prefix = "" if prefix == "." else prefix.upper()

            matrix[prefix] = {"a": float(parts[1]),
                              "c": float(parts[2]),
                              "g": float(parts[3]),
                              "t": float(parts[4])}

            prefixes_prob[prefix] = float(parts[6])

    first_prefix = next(iter(matrix))
    order = len(first_prefix)

    return {"matrix": matrix,
            "prefixes_prob": prefixes_prob,
            "order": order}

######################################
#   Function: Sequence probability   #
######################################
def sequence_probability(sequence, model):
    """
    Compute sequence probability from a Markov model.

    Args:
        sequence (str): sequence
        model(dict): Markov model
    Returns:
        float
    """

    sequence = sequence.upper()

    matrix = model["matrix"]
    prefixes_prob = model["prefixes_prob"]
    order = model["order"]

    if order == 0:
        prob = 1.0
        for base in sequence:
            prob *= matrix[""].get(base.upper(), 0.0)
        return prob

    if len(sequence) < order:
        return 0.0

    prefix = sequence[:order]

    if prefix in prefixes_prob:
        return 0.0

    prob = prefixes_prob.get(prefix, None)
    if prob is None:
        return 0.0

    for i_index in range(order, len(sequence)):
        context = sequence[i_index - order:i_index]
        base = sequence[i_index]

        if context not in matrix:
            return 0.0

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
        # Force the HTML extension
        if not output_file.endswith(".tsv"):
            output_file += ".tsv"
        # If it's a file
        output_path = output_file

    # Write HTML file output
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
        tsv_file.write("seq_id\tseq_proba\n")

        for seq_id, sequence in sequences.items():
            p = sequence_probability(sequence, model)

            tsv_file.write(f"{seq_id}\t{p:.6e}\n")


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

    print("Nb contexts dans modèle:", len(model["matrix"]))
    print("Exemple contexts:", list(model["matrix"].keys())[:10])

#####################
#   Executing code  #
#####################
if __name__ == "__main__":
    main()