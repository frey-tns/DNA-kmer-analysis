"""
Compute a transition matrix of a Markov background model from biological sequences.

This program reads a set of input sequences (provided as a FASTA-formatted file),
and builds a Markov transition matrix (transition matrix) of a specified order:

- transition probabilities P(base | prefix)
- prefix probabilities P(prefix)
- orw sums (should equal 1)
- global nucleotide frequencies (P_res)

The model is estimated from observed k-mer counts (k = m + 1), where m is the
order of the Markov model. For each prefix (context), the program computes the
conditional probabilities of observing each nucleotide (A, C, G, T).

The results are written to a TSV file.

USAGE
    python markov-from-seq -i data/yeast_MET_upstream.fasta -m 2 -o results/
    python markov-model -i data/yeast_MET_upstream.fasta -m 2 -o results/


OPTIONS
    -i : input FASTA file
    -o : output TSV file
    -m : markov model order

AUTHOR
    Anouk RISCH

CONTACT
    https://github.com/frey-tns

URL
    https://github.com/frey-tns/DNA-kmer-analysis

VERSION
    0.1, 2026-04-30

INSTALLATION

    pip3 install -e .

USAGE AND OPTIONS

    For command-line usage, run:
        markov-model --help
        markov-from-seq --help

EXAMPLES

    markov-from-seq -i data/yeast_MET_upstream.fasta -m 2 -o results/

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
from colorama import Fore

############################
#   Internal libraries     #
############################

# Read FASTA file
import bioseq_kmers.sequences as seq
# Markov model background
import bioseq_kmers.background_models as markov_bg
# Format command line
import bioseq_kmers.utils as utils


#################
#   Main code   #
#################

def main():
    """
    Run the complete workflow to build a Markov background model from biological sequences.

    This function parses command-line arguments, loads sequences from a FASTA file,
    computes a Markov model transition matrix of a given order, and write the resulting model
    to a TSV output file.

    The generated output file contains:
        - command line used
        - program version and input file information
        - transition matrix:
            * prefix
            * conditional probabilities P(base | prefix) for A, C, G, T
            * row sum (should be equal to 1)
            * prefix probability P(prefix)
        - global statistics:
            * sum of probabilities per base
            * estimated residue frequencies (P_res)
        - execution timestamps and runtime


    Command-line arguments:
        -i, --input:
            Path to the input FASTA file.

        -m, --markov:
            Order of the Markov model.
            The model uses k-mers of length k = m + 1

        -o, --output:
            Path to the output TSV file path or directory.


    Raises:
        FileNotFoundError: If the input FASTA file does not exist.
        ValueError: If the input FASTA file is malformed.

    Returns:
        None

    Notes:
        - The transition probabilities are estimated from observed k-mer counts.
        - The output format mimics RSAT background model files (e.g. create-background-model).
        - Residue frequencies (P_res) are derived from aggregated transition probabilities

    """

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
    matrix, total_all, context_counts = markov_bg.markov_model(dico_sequences, int(args.markov))

    #####################
    #   Output file     #
    #####################

    # Current date
    today = str(datetime.date.today()).replace("-", "_")

    # If the output path is a folder
    if os.path.isdir(output_file):
        # Define output path
        output_path = os.path.join(output_file, f"bg_model_m{order}_{today}.tsv")
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
        tsv_file.write(f"; markov-from-seq\t{command_line}\n;\n")
        # URL in input
        tsv_file.write(f"; Program version\t{version}\n"
                       f"; Input file\t{os.path.relpath(input_file)}\n")

        # Header
        tsv_file.write("#pr\\suf \ta\tc\tg\tt\tSum\tP_prefix\n")
        
        base = ["A", "C", "G", "T"]

        # Stock sums for final stats
        sum_bases = {b: 0.0 for b in base}

        # Brows the prefixes of the Markov model
        for prefix in sorted(matrix.keys()):
            # Extract conditional probabilities
            prob = matrix[prefix]

            # Initialize values (A, C, G, T)
            row_prob = []
            # Initialize sum of line (=1)
            row_sum = 0.0

            # P(prefix)
            P_prefix = context_counts[prefix] / total_all

            for b in base:
                # Transition probability extract (P(b∣prefix)) (if empty = 0.0)
                p = prob.get(b, 0.0)
                row_prob.append(p)
                # Overall sum (probabilities by base)
                sum_bases[b] += p
                # Row sum (sumP(b∣prefix)=1)
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

