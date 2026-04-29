"""
This library contains functions to handle background models for DNA sequences.

- expected_frequencies()  estimates k-mer probabilities under a Bernoulli model, which assumes
indenpendence between nucleotides along the sequence

In this model, the probability of a k-mer is computed as the product of
individual nucleotide probabilities:

                  k
        P(kmer) = ∏ P * (ri)
                 i=1

This corresponds to a Bernoulli background model (Markov model of order 0),
where no dependency between adjacent nucleotides is assumed.

This function is used in:
- k-mer enrichment / depletion analysis
- statistical significance estimation
- background model comparison (observed vs expected k-mer frequencies)

USAGE

    from exp_freq import expected_frequencies

     exp_freq = {"A": 0.3, "T": 0.2, "G": 0.3, "C": 0.2}

     exp_freq = expected_frequencies(canon_kmer, frequencies)

    print(exp_freq)
    # 0.3 * 0.2 * 0.3 = 0.018

INPUT
    single_kmer : str
        DNA k-mer sequence (A, T, G, C)

    frequencies : dict
        Nucleotide frequency dictionary:
        {"A": float, "T": float, "G": float, "C": float}

OUTPUT
    float
        Expected probability of the k-mer under a Bernoulli model

ERRORS
    KeyError
        If the k-mer contains invalid or unknown nucleotides

AUTHOR
    Anouk RISCH

CONTACT
    https://github.com/frey-tns

URL
    https://github.com/frey-tns/DNA-kmer-analysis

VERSION
    1.2, 2026-04-24
"""

import bioseq_kmers.kmer_stats as kmers
from collections import defaultdict


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

###############################
#   Function : Markov model   #
###############################

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
