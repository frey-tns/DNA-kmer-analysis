"""
Compute statistics about k-mer occurrences in biological sequences.

SYNOPSIS USAGE
    Print usage line:
        kmer-analysis

    Print help message:
        kmer-analysis --help

    usage:
        kmer-analysis [-h] --input FASTA_FILE
         --kmer-length LENGTH_KMER --output OUTPUT_FILE

DESCRIPTION

This program computes the occurrences of each k-mer in a set of input sequences (provided as a fasta-formatted file),
and derives different statistics, depending on the user-specified options:

- observed k-mer occurrences (counts) and relative frequencies
- expected occurrences and frequencies based on either a Bernoulli or a Markov model
- over- or under-representation statistics (binomial P-value, E-value, significance)

For DNA sequences, k-mers can optionally be grouped by pairs of reverse complements.
This applies to the observed counts/frequencies as well as to the derived statistics.

The results are written to a TSV file.


OPTIONS
    -h, --help
        Display this help message and exit.

    -i, --input FASTA_FILE
        Input sequence file in fasta format.

    -k, --kmer KMER_LENGTH
        Length of the k-mer sequences to analyze (integer ≥ 1).

    -s, --strand {single,both}
        Strand mode to compute counts for.
        Handling of the strands for computing the occurrences and derived statistics.

        single: count occurrences on the forward strand only
        both: regroup k-mers by pairs of reverse complements
            This actually makes the counts and derived statistics strand-insensitive

    -o, --output OUTPUT_FILE
        Output file in tab-separated values, with recommended extension .tsv

    -r, --return FIELD[,FIELD...]
        Comma-separated list of output fields.

-r, --return FIELD[,FIELD...]
    Comma-separated list of output fields (one or more).

    Accepted values:
        occ       : observed occurrences
        obs_freq  : observed relative frequency
        exp_occ   : expected occurrences
        exp_freq  : expected frequency

        proba     : statistical significance measures, including:
            p_value : binomial P-value according to the specified background model
            e_value : E-value, i.e. P-value × T, where T is the number of k-mers tested
            sig     : significance score, defined as -log10(E-value)

OUTPUT

A tab-separated value file (extension  .tsv) with ne row per observed k-mer in the input sequences,
and one column per k-mer statistics.

Column contents :
- k-mer sequence
- observed occurrence count
- observed frequency
- expected occurrence count (optional)
- expected frequency (optional)

EXAMPLES

    kmer-analysis -i data/yeast_MET_upstream.fasta  -k 6 -s both -o results/yeast_MET_upstream_6nt_2str.tsv

    kmer-analysis -i data/yeast_MET_upstream.fasta  -k 6 -s both \\
        --return occ,obs_freq,exp_occ,exp_freq \\
        -o results/yeast_MET_upstream_6nt_2str.tsv

AUTHOR / CREDITS
    Anouk RISCH
    supervised and revised by Jacques van Helden

VERSION
    0.2, 2026-05-05

CONTACT / URL
    https://github.com/frey-tns
    https://github.com/frey-tns/DNA-kmer-analysis


"""

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

# Coloring warning text
from colorama import init, Fore
# To reset color
init(autoreset=True)

############################
#   Internal libraries     #
############################

import bioseq_kmers.sequences as seq
import bioseq_kmers.kmer_stats as kmers
import bioseq_kmers.background_models as bg
import bioseq_kmers.utils as utils

################################################################
## CONSTANTS
################################################################
# Default field definitions
_DEFAULT_FIELDS = ["occ", "obs_freq"]

# Dependencies between dictionaries
_DICT_DEPENDENCIES = {"obs_freq" : ["occ"],
                     "exp_occ" : ["exp_freq"],
                     "exp_freq" : []}

################################################################
## FUNCTIONS
################################################################

#####################################
#   Function: Parse return option   #
#####################################

# User parameter for --return
def parse_return_option(return_str):
    """
        Parse the `--return` command-line option into a list of requested output fields.

        The function splits the input string on commas and removes any surrounding
        whitespace from each field. If no option is provided, a default list of fields
        is returned.

        Args:
            return_str (str or None): Comma-separated string of requested fields
                (e.g. "occ,exp_freq,obs_freq"). If None, default fields are used.
        Returns:
            list of str: Cleaned list of requested output fields.
        Examples:
            >>> parse_return_option("occ, exp_freq, obs_freq")
            ['occ', 'exp_freq', 'obs_freq']
            >>> parse_return_option(None)
            ['occ', 'obs_freq']
        """
    # If user give nothing
    if return_str is None:
        # Return occ and obs freq (default field)
        return _DEFAULT_FIELDS

    # For each element f in the list resulting from .split(","),
    # clean with .strip() and add to a new list
    return [f.strip() for f in return_str.split(",")]


##############################
#   Function: Dependencies   #
##############################

# Automatic addition of dependencies
def resolve_dependencies(fields):
    """
      Resolve dependencies between output fields.

      Some output fields require others to be computed beforehand. This function
      expands the list of requested fields by automatically adding all required
      dependencies, based on a predefined dependency dictionary.

      The resolution is performed iteratively until no new dependencies need to be added.

      Args:
          fields (list of str): List of requested output fields.
      Returns:
          list of str: Extended list including all required dependencies.
      Raises:
          KeyError: If a field is not defined in the dependency dictionary.
      Examples:
          >>> resolve_dependencies(['obs_freq'])
          ['obs_freq', 'occ']
          >>> resolve_dependencies(['exp_occ'])
          ['exp_freq', 'exp_occ']
      """
    # Avoids duplicates
    resolved = set(fields)

    change = True
    while change:
        change = False
        # Browse the requested fields
        for field in list(resolved):
            # For each field, retrieve its dependencies (if none → empty list)
            for dep in _DICT_DEPENDENCIES.get(field, []):
                # If no dependency
                if dep not in resolved:
                    # Adding the missing dependency
                    resolved.add(dep)
                    change = True

    return list(resolved)

#################
#   Main code   #
#################
def main():
    """
    Run the complete k-mer analysis workflow.

    The function parses command-line arguments, loads the input FASTA file , computes nucleotide frequencies,
    counts observed k-mer sequences according to the selected strand mode, estimates expected frequencies and
    occurrences, then writes results to a TSV output file.

    The generated report contains:
        - command line used to launch the program
        - input/output parameters
        - sequence summary statistics
        - k-mer analysis table
        - execution timestamps and runtime

    Command-line arguments:
        -i, --input:
            Path to the input FASTA file.

        -k, --kmer:
            Length of the K-mer sequences to analyse.

        -o, --output:
            Output TSV file path or output directory.

        -s, --strand:
            Strand mode to compute counts for.

        -r, --return:
            Output TSV file path or output directory.

    Raises:
        FileNotFoundError: If input FASTA file does not exist.
        ValueError: If input FASTA file is empty.
        OSError: If the output directory cannot be created.

    Returns:
        None
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

    parser.add_argument("-r","--return",
                        dest = "return_fields",
                        type = str,
                        default = None,
                        help = "Comma-separated list of fields to return. Default: occ,freq. Supported values: occ,obs_freq,exp_occ,exp_freq")
    # Reads the command typed in the terminal
    args = parser.parse_args()

    # Define variable to use the value in the script
    input_file = args.input
    kmer_length = args.kmer_length
    output_file = args.output
    strand_mode = args.strand
    requested_fields = parse_return_option(args.return_fields)

    fields_compute = resolve_dependencies(requested_fields)


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

    sequences, total_length, seq_number = seq.read_fasta(fasta_file)

    ###################
    #    Statistics   #
    ###################

    # Nucleotides frequency
    frequencies = kmers.nucleotide_frequencies(sequences)

    observed_kmer_count = kmers.counts_kmer(sequences, kmer_length, strand_mode)
    # Number of all positions T = L - K + 1
    total_positions = sum(len(seq) - kmer_length + 1 for seq in sequences.values())

    # List who contain output results
    result_analysis = []

    # Treat only once the same canonical kmer
    dico_canon_kmer = set()

    # Browses the dictionary of observed k-mers.
    for canon_kmer, occ in observed_kmer_count.items():

        if strand_mode == "both":
            canon_kmer = kmers.canonic_kmer(canon_kmer)
            # Kmer inverse complement
            reverse_kmer = kmers.reverse_complementary(canon_kmer)
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

        # Occurrence threshold
        threshold = 1
        if occ < threshold:
            continue

        row = {"seq" : canon_kmer,
                "id" : kmer_id}

        # Occurrence
        if "occ" in fields_compute:
            row["occ"] = occ

        # Expected frequencies
        if "exp_freq" in fields_compute:
            exp_freq = bg.expected_frequencies(canon_kmer, frequencies)
            row["exp_freq"] = exp_freq

        # Expected occurrences
        if "exp_occ" in fields_compute:
            # Checks if exp_freq has already been calculated
            if "exp_freq" not in row:
                # exp_freq does not yet exist
                exp_freq = bg.expected_frequencies(canon_kmer, frequencies)
            else:
                # exp_freq already exists
                exp_freq = row["exp_freq"]

            row["exp_occ"] = total_positions * exp_freq

        # Observed frequencies
        if "obs_freq" in fields_compute:
            row["obs_freq"] = occ / total_positions

        # Stockage in list oligomers
        result_analysis.append(row)

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
        command_line = utils.format_command_line(sys.argv)
        # Write command line
        tsv_file.write(f"; Command\t{command_line}\n;\n")
        # URL in input
        tsv_file.write(f"; Fasta URL\t{fasta_file}\n"
                       f"; Input file\t{os.path.relpath(input_file)}\n"
                       f"; Input format\tFasta\n"
                       f"; Output file\t{os.path.relpath(output_path)}\n"
                       f"; Strand mode\t{strand_mode}\n"
                       f"; Oligomer length\t{kmer_length}\n;\n")

        # Summary
        tsv_file.write(f"; Sequence type\tDNA\n"
                       f"; Nb of sequence\t{len(sequences)}\n"
                       f"; Sum of sequence lengths\t{total_length}\n"
                       f"; Sequences:\n")
        for seq_id, seq_fasta in sequences.items():
            tsv_file.write(f";\t{seq_id}\t{len(seq_fasta)}\n")

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
        column = ["seq","id"] + requested_fields
        tsv_file.write(f"# {'\t'.join(column)}\n")

        for row in result_analysis:
            # kmer analysis table
            # Initializes an empty list to contain all the values of the row
            list_values = []
            # Browses the columns
            for col in column:
                # Retrieves the value associated with the column, if the key doesn't exist → empty string
                value = row.get(col, "")

                if col == "exp_occ" and isinstance(value, float):
                    # Rounded exp occurrence
                    value = f"{value:.2f}"
                else:
                    # Else convertion to string
                    value = str(value)

                list_values.append(value)

            # Assemblies of values with tabs
            tsv_file.write(f"{'\t'.join(list_values)}\n")

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


