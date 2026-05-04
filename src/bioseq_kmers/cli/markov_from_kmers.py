"""
Estimate the parameters of a Markov background model from a kmer counting table and output a transition matrix.

SYNOPSIS USAGE
    Print usage line:
        markov-from-kmers

    Print help message:
        markov-from-kmers --help

    Usage:
        markov-from-kmers [-h] --input KMER_TABLE
           --markov MARKOV_ORDER --output OUTPUT_FILE

DESCRIPTION

This program reads a kmer counting table (tab-separated file), and builds a Markov model
(transition matrix) with a user-specified order. The order of the Markov model corresponds to the size of the
prefix (oligonucleotide) to compute the conditional probabilities P(base|prefix).

The input file is expected to contain at least:
- a header line starting with "#"
- a column describing k-mer sequences (e.g. "seq" or "sequence")
- a column describing k-mer occurrences (e.g. "occ" or "occurrence")

For a Markov model of order m:
- prefixes have length m
- input k-mers must have length m + 1

OPTIONS
    -h, --help
        Display this help message and exit.

    -i, --input KMER_TABLE
        Input kmer counting table.

    -m, --markov MARKOV_ORDER
        Markov model order.

    -o, --output OUTPUT_FILE
        Output TSV file.

OUTPUT

One row per prefix of the Markov model

Column contents :
- prefix
- transition probabilities P(base|prefix)
- prefix probabilities P(prefix)
- row sums (for validation, should equal 1)
- global nucleotide frequencies (P_res)

The model is estimated from  k-mer counts (k = m + 1) observed in the input sequences,

The results are written to a tab-separated value file (extension .tsv).

EXAMPLES

    markov-from-kmers -i data/3nt_counts.tsv -m 2 -o results/
    markov-from-kmers -i data/1nt_counts.tsv -m 0 -o results/

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
# For shell command string manipulation
import shlex
import sys

from collections import defaultdict
# Coloring warning text
from colorama import init, Fore

############################
#   Internal libraries     #
############################

import bioseq_kmers.background_models as markov_bg
import bioseq_kmers.utils as utils

################################################################
## FUNCTIONS
################################################################

#######################################
#   Function: reads a table of kmer   #
#######################################
def load_kmer_counting_table(path):
    """
    Read a k-mer counting table in RSAT tabular format and extract occurrence counts.

    This function parses a tab-separated file containing k-mer statistics produced by
    kmer analysis tools (e.g. RSAT oligo-analysis or kmer-analysis).

    The input file must contain:
    - a header line starting with "#"
    - a tab-separated file containing k-mer column (e.g; 'seq' or 'sequence')
    - a tab-separated file containing occurrence column (e.g. 'occ' or 'occurrence')

    Args:
        path (str): Path to the input k-mer counting table file
    Returns:
        dict : Dictionary mapping k-er sequence to occurrence counts.
            {kmer(str): occ(int)}.
    Raises:
        ValueError: If the input file does not contain a header line starting with "#" or header line
                    missing required columns ('seq', 'occ').
    """

    # Final result
    dico_counts = {}
    # Dict where key is column name and values is the index
    col_index = None

    # Read file
    with open(path,'r') as f:
        for line in f:
            # Remove the line breaks
            line = line.strip()

            # Ignore empty lines or line who start with ; (RSAT comment)
            if not line or line.startswith(";"):
                continue

            if line.startswith("#"):
                # Create list who contain the table header
                header = line.lstrip("#").strip().split()
                # Associates the column name with its index
                col_index = {name: i for i, name in enumerate(header)}
                continue

            # If missing header
            if col_index is None:
                raise  ValueError("Missing header line start with #")

            # Line parsing
            parts = line.split("\t")

            # Fold over if spaces are used instead of tabs
            if len(parts) == 1:
                parts = line.split()

            # Extract kmer
            if "seq" in col_index:
                kmer = parts[col_index["seq"]]
            elif "sequence" in col_index:
                kmer = parts[col_index["sequence"]]
            # Extract occurrence
            if "occ" in col_index:
                occ = int(parts[col_index["occ"]])
            elif "occurrence" in col_index:
                occ = int(parts[col_index["occurrence"]])
            else:
                # If no column find → error
                raise ValueError("No occurrence column found (occ)")

            dico_counts[kmer] = occ

    return dico_counts

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
    parser = argparse.ArgumentParser(description="Markov model")

    # Define args used by the user (input, output, Markov order)
    parser.add_argument("-i", "--input",
                        required=True,
                        help="input fasta file")

    parser.add_argument("-m", "--markov",
                        required=True,
                        help="Order of the markov model (1-7)")

    parser.add_argument("-o", "--output",
                        required=False,
                        help="output TSV file")

    # Reads the command typed in the terminal
    args = parser.parse_args()

    ### Define variable to use the value in the script
    # Load sequence
    output_file = args.output
    input_file = args.input
    order = int(args.markov)

    # Build Markov model
    kmer_counts = load_kmer_counting_table(input_file)

    # Check that all kmer have the same length
    kmer_length = {len(kmer) for kmer in kmer_counts}

    if len(kmer_length) != 1:
        raise ValueError("Input file must contain kmers of identical length.")

    kmer_length = kmer_length.pop()

    # Check compatility between kmer length and Markov order
    expected_order = kmer_length - 1

    if expected_order != order:
        raise ValueError(f"Incompatible order: input kmers have length {kmer_length},"
                         f"so the compatible Markov order is {expected_order} but got {order}.")

    # Markov model from kmers
    matrix, total_all, context_counts = markov_bg.markov_from_kmers(kmer_counts, order)

    #####################
    #   Output file     #
    #####################

    # Current date
    today = str(datetime.date.today()).replace("-", "_")

    # If the output path is a folder
    if os.path.isdir(output_file):
        # Define output path
        output_path = os.path.join(output_file, f"markov_transitions_m{order}_{today}.tsv")
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
        tsv_file.write(f"; markov-from-kmers\t{command_line}\n;\n")
        # URL in input
        tsv_file.write(f"; Program version\t{version}\n"
                       f"; Input file\t{os.path.relpath(input_file)}\n")

        # Header
        tsv_file.write("#pr\\su \ta\tc\tg\tt\tSum\tP_prefix\n")

        base = ["A", "C", "G", "T"]

        # Stock sums for final stats
        sum_bases = {b: 0.0 for b in base}

        for prefix in sorted(matrix.keys()):
            prob = matrix[prefix]

            row_prob = []
            row_sum = 0.0

            display_prefix = prefix if prefix != "" else "."

            # P(prefix)
            P_prefix = context_counts[prefix] / total_all

            for b in base:
                p = prob.get(b, 0.0)
                row_prob.append(p)
                sum_bases[b] += p
                row_sum += p

            tsv_file.write(
                f"{display_prefix.lower()}\t{'\t'.join(f'{p:.5f}' for p in row_prob)}\t{row_sum:.0f}\t{P_prefix:.4f}\n")

            # Global sums
        tsv_file.write(f"; Sum\t"
                       f"{'\t'.join(f'{sum_bases[b]:.5f}' for b in base)}\t"
                       f"{len(matrix)}\n")

        # Residue frequencies
        total_sum = sum(sum_bases.values())

        tsv_file.write(f"; P_res\t"
                       f"{'\t'.join(f'{sum_bases[b] / total_sum:.5f}' for b in base)}\n")


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