__version__ = 0.1
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

###########################
#   Function: load kmer   #
###########################
def load_kmer_counts(path):

    counts = {}
    with open(path,'r') as f:
        for line in f:
            kmer, count = line.split()
            counts[kmer] = int(count)

    return counts

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
    order = args.markov

    # Build Markov model
    kmer_counts = load_kmer_counts(input_file)

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
        tsv_file.write(f"; Program version\t{__version__}\n"
                       f"; Input file\t{os.path.relpath(input_file)}\n")

        # Header
        tsv_file.write("#pr\\suf \ta\tc\tg\tt\tSum\tP_prefix\n")

        base = ["A", "C", "G", "T"]

        # Stock sums for final stats
        sum_bases = {b: 0.0 for b in base}

        for prefix in sorted(matrix.keys()):
            prob = matrix[prefix]

            row_prob = []
            row_sum = 0.0

            # P(prefix)
            P_prefix = context_counts[prefix] / total_all

            for b in base:
                p = prob.get(b, 0.0)
                row_prob.append(p)
                sum_bases[b] += p
                row_sum += p

            tsv_file.write(
                f"{prefix.lower()}\t{'\t'.join(f'{p:.5f}' for p in row_prob)}\t{row_sum:.0f}\t{P_prefix:.4f}\n")

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