"""
Nucleotide frequency analysis utilities for DNA sequence datasets.

This module provides functions to compute nucleotide composition
statistics from DNA sequences. It is designed as a reusable library
for downstream analyses such as k-mer statistics, background model
estimation (Bernoulli/Markov models), and sequence probability computation
(e.g. kmer_analysis.py, markov_model.py, and sequence_prob.py).

Main functionality:

- Computation of nucleotide frequencies (A, T, G, C)
- Handling of ambiguous nucleotides (ignored, e.g. N)

The frequencies are computed globally across all input sequences
by counting valid nucleotide (A, C, T, G) occurrences and normalizing by the total
number of observed valid bases.

This module is typically used as a preprocessing step for:
- k-mer expected frequency estimation
- Markov model parameter inference
- sequence probability estimation

USAGE

    from sequences import read_fasta
    from nucleotide_frequencies import nucleotide_frequencies

    sequences, _ = read_fasta("input.fasta")

    frequencies = nucleotide_frequencies(sequences)

    print(frequencies)

    Example output:
    >>> {'A': 0.31, 'T': 0.29, 'G': 0.20, 'C': 0.20}

INPUT
    sequences : dict
        Dictionary of DNA sequences in the form:
        {sequence_id: sequence_string}

OUTPUT
    dict
        Dictionary of nucleotide relative frequencies:
        {"A": float, "T": float, "G": float, "C": float}

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

# Nucleotide count
from collections import Counter

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