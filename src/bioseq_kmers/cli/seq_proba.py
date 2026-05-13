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

    seq-proba -i data/seq/yeast_MET_upstream.fasta \
        -m data/bg-models/yeast_all-upstream-noorf_Markov_m2.tsv \
        -o results/yeast_MET_upstream_proba.tsv

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
import bioseq_kmers.background_models as bg

################################################################
## FUNCTIONS
################################################################


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

    parser.add_argument("-m", "--model",
                        required=True,
                        help="Markov model, formated as a transition matrix in a TSV file")

    parser.add_argument("-o", "--output",
                        required=False,
                        help="output TSV file")

    # Reads the command typed in the terminal
    args = parser.parse_args()

    ### Define variable to use the value in the script
    # Load sequence
    input_file = args.input
    model_file = args.model
    output_file = args.output

    if output_file is None:
        output_file = "."

    ############################
    #   Load data
    ############################

    # sequences
    sequences, _, _ = seq.read_fasta(input_file)

    # Markov model
    model = bg.load_markov_matrix(model_file)

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


        dict_proba_seq, dict_log_proba = bg.sequence_probability(sequences, model)

        for seq_id in dict_proba_seq.keys():
            sequence = sequences[seq_id]
            length_seq = len(sequence)
            log_p = dict_log_proba[seq_id]
            p = utils.engineer_mode(log_p) # enables printing values smaller than 1e-324

            tsv_file.write(f"{seq_id}\t{length_seq}\t{p}\t{log_p:.2f}\n")

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