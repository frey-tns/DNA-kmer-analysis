"""
K-mer and nucleotide statistics utilities for DNA sequence analysis.

This module provides core functions to manipulate and analyze k-mers from biological DNA sequences.
It is designed to be used as a reusable library for downstream analyses such as k-mer statistics, background model
estimation (Bernoulli/Markov models), and sequence probability computation
(e.g. kmer_analysis.py, markov_model.py, and sequence_prob.py).

Main functionalities include :

    - Reverse complement computation of DNA k-mers
    - Canonical k-mer representation (strand-invariant k-mers)
    - K-mer counting in DNA sequences with optional strand handling
    - Nucleotide frequency estimation from sequence datasets
    - Sequence composition (nucleotide frequencies)

These functions are used as building blocks for higher-level tools such as :

    - k-mer frequency analysis
    - Markov model parameter estimation
    - sequence probability computation under background models (e.g. Bernoulli, Markov)

All functions are designed to be modular, independent of Input/Output,
and suitable for integration into bioinformatics workflows.

The canonical k-mer is a standardized representation of a DNA k-mer that accounts for the double-stranded nature of DNA.

USAGE
    This module is intended to be imported and used in other scripts:

    from kmer_counts import reverse_complementary, canonic_kmer, counts_kmer, nucleotide_frequencies

    # Load sequences
    from sequences import read_fasta
    sequences, total_length, seq_number = read_fasta("yeast_MET_upstream.fasta")

    # Reverse complement
    reverse_kmer = reverse_complementary("ATGC")

    # Canonical k-mer
    canon_kmer = canonic_kmer("ATGC")

    # K-mer counting
    observed_kmer_count = counts_kmer(kmer_length,
                                      sequences,
                                      strand_mode)

    # Compute nucleotide frequencies
    frequencies = nucleotide_frequencies(sequences)


AUTHOR
    Anouk RISCH

CONTACT
    https://github.com/frey-tns

URL
    https://github.com/frey-tns/DNA-kmer-analysis

VERSION
    1.2, 2026-04-24

"""

#################
#   Libraries   #
#################
# Progress bar
from tqdm import tqdm
# Count (cf k-mer count)
from collections import Counter

################################################################
## CONSTANTS
################################################################

# Create complementary reverse dictionary
_DICT_COMPLEMENT = {"A": "T",
                    "T": "A",
                    "C": "G",
                    "G": "C"}

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

    # Trace the Kmer backwards and replace each base with its complement
    return "".join(_DICT_COMPLEMENT[base] for base in reversed(kmer_seq))


#####################################
#      Function: Canonical Kmer     #
#  Merge kmer + reverse complement  #
#####################################
def canonic_kmer(kmer):
    """
    Return  the canonical representation of a DNA k-mer.

    For any k-mer, its reverse complement represents the same biological signal on the opposite strand.
    To avoid double counting, the canonical k-mer is defined as:

    the shortest lexicographic sequence between:
        - the k-mer itself
        - its reverse complement

    This ensures that both strands contribute to a single unified representation,
    which is essential for unbiased k-mer counting and downstream statistical analyses.

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


#####################################
#   Function : Count observed kmer  #
#####################################

def counts_kmer(sequence, k_length, strand_mode):
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