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

import bioseq_kmers.sequences as seq
import bioseq_kmers.kmer_stats as kmers

################################################################
## FUNCTIONS
################################################################

###################################
#   Function: transition matrix   #
###################################

def markov_model(sequences, order):
    """
    Build a markov transition matrix of order m from DNA sequences.

    Args:
        sequences (dict): {id: sequence}.
        Order of the markov model :
            order(int) = m = k–1 = 1
    Example:
        >>> sequences = {"s1":"ATCGT"}
        >>> markov_model(sequences, order=1)
        {'A': {'T':1.0}, 'T': {'C':1.0}, 'C': {'G':1.0}, 'G': {'T':1.0}}
    """
    length_kmer = order + 1

    # Create a two-level dictionary (default value = 0.0)
    transition_matrix = defaultdict(lambda: defaultdict(float))

    # Number of occurrences of each k-mer in the sequences (k = m+1)
    kmer_counts = kmers.counts_kmer(sequences, length_kmer, strand_mode="single")

    # The dictionary stores the total number of occurrences per context
    context_counts = defaultdict(int)

    for kmer, count in kmer_counts.items():
        # Extract prefix of Markov model
        prefix = kmer[:-1]
        # Next base extraction
        suffix = kmer[-1]

        # Order occurrences
        context_counts[prefix] += count
        # Increments a matrix cell
        transition_matrix[prefix][suffix] += count

    # Standardization (transforms counts into probabilities)
    for prefix in transition_matrix :
        # Total for this order
        total = context_counts[prefix]
        # Scans all the databases observed for this context
        for base in transition_matrix[prefix]:
            # Transforms an occurrence into a probability
            transition_matrix[prefix][base] /= total

    total_all = sum(context_counts.values())

    return dict(transition_matrix), total_all, context_counts

#################################################
#   Function: Defined the format command line   #
#################################################
def format_command_line(argv):
    """
     Format command line by changing the absolute path to a relative path.

     Args:
         argv (list): Command-line arguments.

     Returns:
         str: Reconstructed command line.
     """
    # Retrieves the folder from which the script is executed
    cwd = os.getcwd()
    # Contains the rebuilt command
    list_cleaned_command = []

    # The next argument is a file path
    skip_next = False

    # Iterates through each element of the order
    for arg in argv:

        # If argument associated with -i or -o
        if skip_next:
            # Converts absolute path to relative path
            rel = os.path.relpath(arg, cwd)
            list_cleaned_command.append(shlex.quote(rel))
            # Returns to the initial state
            skip_next = False

        elif arg in ["-i", "--input", "-o", "--output"]:
            # Keep the current flag
            list_cleaned_command.append(arg)
            # The next argument is a path
            skip_next = True

        else:
            # Normal argument
            list_cleaned_command.append(shlex.quote(arg))

    return " ".join(list_cleaned_command)


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
    dico_sequences, _, _ = seq.read_fasta(args.input)
    output_file = args.output
    input_file = args.input
    order = args.markov


    # Build Markov model
    matrix, total_all, context_counts = markov_model(dico_sequences, int(args.markov))

    #####################
    #   Output file     #
    #####################

    # Current date
    today = str(datetime.date.today()).replace("-", "_")

    # If the output path is a folder
    if os.path.isdir(output_file):
        # Define output path
        output_path = os.path.join(output_file, f" bg_model_m{order}_{today}.tsv")
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
        command_line = format_command_line(sys.argv)
        # Write command line
        tsv_file.write(f"; markov-from-seq\t{command_line}\n;\n")
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


            tsv_file.write(f"{prefix.lower()}\t{'\t'.join(f'{p:.5f}' for p in row_prob)}\t{row_sum:.0f}\t{P_prefix:.4f}\n")

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

