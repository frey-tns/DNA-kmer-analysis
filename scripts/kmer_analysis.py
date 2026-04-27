"""
Compute statistics about k-mer occurrences in biological sequences.

This program counts the occurrences of all k-mers in a set of input sequences (provided as a fasta-formatted file),
and derives different optional statistics :

- observed k-mer occurrence counts and relative frequencies
- expected occurrences and frequencies based on either a Bernoulli or a Markov model
- over- or under-representation statistics (P-value, E-value, significance)

For DNA sequences, k-mers can optionally be grouped by pairs of reverse-complements. This applies to the observed
counts/frequencies as well as to the derived statistics.

The results are written to a TSV file.

USAGE
    python kmer-analysis.py -i input.fasta -k 6 -o output.tsv

OPTIONS
    -i : input FASTA file
    -o : output TSV file
    -k : k-mer size
    -s : strand mode (single or both)
    -r : parse return option ("occ, exp_occ, exp_freq, obs_freq")

AUTHOR
    Anouk RISCH

CONTACT
    https://github.com/frey-tns

URL
    https://github.com/frey-tns/DNA-kmer-analysis

VERSION
    1.2, 2026-04-24

EXAMPLES
    python3 scripts/kmer_analysis.py -i data/yeast_MET_upstream.fasta  -k 6 -s both -o results/ueast_MET_upstream_6nt_2str.tsv

USAGE AND OPTIONS

    For command-line usage, run:
        python3 scripts/kmer_analysis.py --help


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
# For shell command string manipulation
import shlex

from tqdm import tqdm
from collections import Counter # count (cf k-mer count)
# Coloring warning text
from colorama import init
from colorama import Fore
# from colorama import Fore as _Fore
# To reset color
init(autoreset=True)

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
          >>> resolve_dependencies(['exp_freq'])
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
        file_path (str): FASTA file.

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

    The canonical k-mer is defined as the lexicographically smallest sequence between the input k-mer and
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
    The expected frequency is computed as the product of the individual nucleotide frequencies,
    assuming independence between positions.

    For example:
        P(ATG) = P(A) * P(T) * P(G)

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
                        help = "Comma-separated list of fields to return (default: None)")
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

        if strand_mode == "both":
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
            exp_freq = expected_frequencies(canon_kmer, frequencies)
            row["exp_freq"] = exp_freq

        # Expected occurrences
        if "exp_occ" in fields_compute:
            # Checks if exp_freq has already been calculated
            if "exp_freq" not in row:
                # exp_freq does not yet exist
                exp_freq = expected_frequencies(canon_kmer, frequencies)
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
        command_line = format_command_line(sys.argv)
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


