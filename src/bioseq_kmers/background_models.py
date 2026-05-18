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

#################
#   Libraries   #
#################
import math
from tqdm import tqdm
from collections import defaultdict

############################
#   Internal libraries     #
############################
import bioseq_kmers.kmer_stats as kmers
import bioseq_kmers.utils as utils


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

#########################################
#   Function : Markov model from kmers  #
#########################################
def markov_from_kmers(kmer_counts, order):
    """
    Build a Markov transition matrix from precomputed k-mer counts.

    Args:
        kmer_counts (dict): dictionary {kmer: count}
        order (int): Markov order (m), where k = m + 1

    Returns:
        dict: transition matrix {prefix: {base: probability}}
        int: total counts
        dict: context counts
    """
    # Create a two-level dictionary (default value = 0.0)
    transition_matrix = defaultdict(lambda: defaultdict(float))

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
    for prefix in transition_matrix:
        # Total for this order
        total = context_counts[prefix]

        if total == 0:
            continue

        # Scans all the databases observed for this context
        for base in transition_matrix[prefix]:
            # Transforms an occurrence into a probability
            transition_matrix[prefix][base] /= total

    total_all = sum(context_counts.values())

    return dict(transition_matrix), total_all, context_counts

############################################
#   Function : Markov model from sequence  #
############################################
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

    # Number of occurrences of each k-mer in the sequences (k = m+1)
    kmer_counts = kmers.counts_kmer(sequences, length_kmer,
                                    strand_mode="single",
                                    desc=f"Building Markov model (order {order})")


    return markov_from_kmers(kmer_counts, order)


##################################################
#   Function: load markov model from a TSV file  #
##################################################
def load_markov_matrix(path):
    """
    Read a markov transition matrix from a tabular RSAT-like file.

    The function reads a Markov background model produced by tolls such as `markov-from-seq` or
    `markov-from-kmers`. It parses a tab separated file containing prefix, transition probabilities
    ofr each nucleotide, and prefix probabilities.

    The prefix columns can be named either `pr\\su` or 'pr\\suf`.
    The Markov order is inferred from the length of the observed prefixes.

    Args:
        path (str): Path to te Markov transition matrix file (.tsv).
    Returns:
        dict: A dictionary containing:
        - dict_matrix (dict): Transition probabilities such taht:
        matrix[prefix][base] = P(base | prefix)
        - dict_prefix_prob (dict): Probability of each prefix P(prefix)
        - dict_order (int): Markov order (length of prefix)
    Raises:
        ValueError: If the file is malformed, missing required columns,
        or contains no valid matrix data.
    Notes:
        - Prefix are converted to uppercase internally
        - '.' is interpreted as an empty prefix (order 0 mode)
    """
    # Dict transition probabilities
    dict_matrix = {}
    # Dict where key is base and value is the initial probability per base
    dict_prefixes_prob = {}
    # column index
    col_index = None

    # Read file
    with open(path, 'r') as f:
        for line in f:
            # Remove line breaks and spaces
            line = line.strip()

            # Skip empty lines and RSAT comments
            if not line or line.startswith(";"):
                continue

            # Read header
            if line.startswith("#"):
                header = line.lstrip("#").strip().split()
                col_index = {name: i for i, name in enumerate(header)}
                # The header line should not be treated as a data line.
                continue

            # If header doesnt exist
            if col_index is None:
                # Error message
                raise ValueError("Missing header line starting with '#'.")

            # Columns are separated by tabs
            parts = line.split("\t")
            # If TSV malformed
            if len(parts) == 1:
                # Automatically separates on any space
                parts = line.split()

            # Ignore the comment lines
            if parts[0].startswith(";"):
                continue

            # Prefix column selection.
            if "pr\\su" in col_index:
                prefix_col = "pr\\su"
            elif "pr\\suf" in col_index:
                prefix_col = "pr\\suf"
            else:
                # Error message if prefix is missing
                raise ValueError("Missing prefix column ('pr\\suf' or 'pr\\su').")

            # Required columns
            for col in ["a", "c", "g", "t", "P_prefix"]:
                if col not in col_index:
                    # Error message if not required columns
                    raise ValueError(f"Missing column '{col}'.")

            # Read prefix
            prefix = parts[col_index[prefix_col]]
            if prefix == ".":
                prefix = ""
                dict_prefixes_prob[prefix] = 1

            else:
                prefix = prefix.upper()

                # Reading prefix probabilities
                dict_prefixes_prob[prefix] = float(parts[col_index["P_prefix"]])

            # Reading transition probabilities
            dict_matrix[prefix] = {"A": float(parts[col_index["a"]]),
                              "C": float(parts[col_index["c"]]),
                              "G": float(parts[col_index["g"]]),
                              "T": float(parts[col_index["t"]])}

    if not dict_matrix:
        raise ValueError("No markov transition matrix data found.")

    # Extract first prefix
    first_prefix = next(iter(dict_matrix))
    # Deducing the order of the model
    order = len(first_prefix)

    return {"matrix": dict_matrix,
            "prefixes_prob": dict_prefixes_prob,
            "order": order}

######################################
#   Function: Sequence probability   #
######################################
def sequence_probability(sequences, model):
    """
    Compute sequence probability of a biological sequence given a Markov background model.

    The probability is computed using an order k Markov chain, where k is defined by the model order.
    For order 0, bases are assumed independent. For higher orders, the probability is computed as:

        P(s) = P(prefix) × Π P(base_i | context_i)

    where context_i is the kmer preceding each base.

    Args:
        sequences (dict): Dictionary of biological sequence {seq_id:sequence(str)}.
        model (dict): Markov model containing:
            - "matrix" (dict): Transition probabilities such that
            matrix[context][base] = P(base_i | context_i)
            - "prefixes_prob" (dict): Probability of each prefix P(prefix)
            -"order" (int): Markov order (length of prefix(kmer))
    Returns:
        dict: A dictionary containing result for each sequence.
        {seq_id: {"probability": scientific notation string (str),
                  "log_proba: log10 probability (float)}
                  }

    Raises:
        ValueError: if a prefix is missing from the background model.

    Notes:
        - If order == 0, bases are assumed independent (Bernoulli).
        - All sequences are normalized to uppercase internally.
        - Probability is computed in log10 space for numerical stability.
    """

    # Extract transition matrix
    matrix = model["matrix"]

    # Extract prior probabilities of the prefixes
    prefixes_prob = model["prefixes_prob"]

    # Extract the order of the model
    order = model["order"]

    dict_proba = {}
    dict_log_proba = {}

    for seq_id, sequence in tqdm(sequences.items(), desc="Computing sequence probabilities"):

    # For computational precision and efficiency, we first compute log_proba = log10(p) and then derive p = 10^log_proba

        # Sequence normalisation
        sequence = sequence.upper()

        # If the sequence is too short or an empty sequence, the probability is null
        if len(sequence) < order or len(sequence) == 0:

            dict_proba[seq_id] = 0.0
            dict_log_proba[seq_id] = float("-inf")


        ## Special case: Bernoulli model (Markov order==0)
        if order == 0:
            log_proba = 0

            for base in sequence:
                # Independently multiplies each base
                residue_proba = matrix[""].get(base.upper(), 0.0)

                if residue_proba == 0:

                    dict_proba[seq_id] = 0.0
                    dict_log_proba[seq_id] = float("-inf")

                log_proba += math.log10(residue_proba)

            dict_proba[seq_id] = 10**log_proba
            dict_log_proba[seq_id] = float("-inf")

        # Extract first prefix
        initial_prefix = sequence[:order]

        # Check that the initial sequence prefix is present in the Markov model
        if initial_prefix not in prefixes_prob:
            raise ValueError(
                f"Initial sequence prefix '{initial_prefix}' is absent from the transition matrix."
            )

        # Initialize log_proba
        prefix_prob = prefixes_prob[initial_prefix]
        if prefix_prob == 0:
            dict_proba[seq_id] = 0.0
            dict_log_proba[seq_id] = float("-inf")
        log_proba = math.log10(prefix_prob)

        # Iterate over the residues following the prefix, and aggregate the transition probabilities
        # JvH: to check: are the start and en indices correct ?
        for i_index in range(order, len(sequence)):
            # Extract prefix
            prefix_start = i_index - order
            prefix_end = i_index
            prefix = sequence[prefix_start:prefix_end]

            # Next base
            base = sequence[i_index]

            # If no prefix
            if prefix not in matrix:
                raise ValueError(
                    f"Prefix '{prefix}' absent from the transition matrix but found at index {i_index} of sequence '{sequence}'."
                )

            # Markov multiplication P(base∣context)
            residue_proba = matrix[prefix].get(base, 0.0)

            if residue_proba == 0:
                dict_proba[seq_id] = 0.0
                dict_log_proba[seq_id] = float("-inf")

            log_proba += math.log10(residue_proba)

#        scientific_prob = utils.engineer_mode(log_proba)

        dict_proba[seq_id] = 10 ** log_proba
        dict_log_proba[seq_id] = log_proba

    return dict_proba, dict_log_proba
