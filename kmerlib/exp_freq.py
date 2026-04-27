"""
Expected k-mer frequency computation under a Bernoulli (0-order Markov) model.

This module provides a function to estimate the expected frequency (probability)
of a DNA k-mer under the assumption that nucleotides occur independently.

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
