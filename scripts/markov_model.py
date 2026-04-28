
#################
#   Libraries   #
#################
from collections import defaultdict

############################
#   Internal libraries     #
############################

from kmerlib.kmer_stats import (counts_kmer)

###################################
#   Function: transition matrix   #
###################################

def markov_model(sequences, order):
    """
    Build a markov transition matrix of order m from DNA sequences.

    Args:
        sequences (dict): {id: sequence}.
        order(int) = m = order of the markov model
    Example:
        >>> sequences = {"s1":"ATCGT"}
        >>> markov_model(sequences, order=1)
        {'A': {'T':1.0}, 'T': {'C':1.0}, 'C': {'G':1.0}, 'G': {'T':1.0}}
    """
    length_kmer = order + 1

    # Create a two-level dictionary (default value = 0.0)
    transition_matrix = defaultdict(lambda: defaultdict(float))

    # Number of occurrences of each k-mer in the sequences (k = m+1)
    kmer_counts = counts_kmer(sequences, length_kmer)

    # The dictionary stores the total number of occurrences per context
    order_counts = defaultdict(int)

    for kmer, count in kmer_counts.items():
        # Extract order of Markov model
        order = kmer[:-1]
        # Next base extraction
        next_base = kmer[-1]

        # Order occurrences
        order_counts[order] += count
        # Increments a matrix cell
        transition_matrix[order][next_base] += count

    # Standardization (transforms counts into probabilities)
    for context in transition_matrix :
        # Total for this order
        total = order_counts[context]
        # Scans all the databases observed for this context
        for base in transition_matrix[context]:
            # Transforms an occurrence into a probability
            transition_matrix[context][base] /= total

    return dict(transition_matrix)
