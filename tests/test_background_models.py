#################
#   Libraries   #
#################
import pytest
import math

############################
#   Internal libraries     #
############################
import bioseq_kmers.background_models as bg
from bioseq_kmers.background_models import sequence_probability


#######################
#   Functional test   #
#######################
def test_expected_frequencies():
    freqs = {"A": 0.3,
             "C": 0.2,
             "G": 0.3,
             "T": 0.2}

    result = bg.expected_frequencies("ATG", freqs)

    assert result == 0.018

def test_expected_frequencies_invalid_bases():
    freqs = {"A": 0.3,
             "C": 0.2,
             "G": 0.3,
             "T": 0.2}

    with pytest.raises(KeyError):
        bg.expected_frequencies("ATN", freqs)


def test_markov_from_kmers():
    kmer_counts = {"AA": 2,
                   "AC": 1,
                   "CA": 1}

    matrix, total, context_counts = bg.markov_from_kmers(kmer_counts, order=1)

    assert matrix == {"A": {"A": 2/3, "C": 1/3},
                      "C":{"A": 1.0}}
    assert total == 4
    assert dict(context_counts) == {"A": 3, "C": 1}

def test_markov_model():
    sequence = {"seq1": "ATCG"}

    matrix, total, context_counts = bg.markov_model(sequence, order=1)

    assert matrix == {"A": {"T": 1.0},
                      "T": {"C": 1.0},
                      "C": {"G": 1.0}}
    assert total == 3
    assert dict(context_counts) == {"A": 1,
                                    "T": 1,
                                    "C": 1}

def test_sequence_probability_order_0():
    model = {"matrix": {"": {"A": 0.20,
                             "C": 0.15,
                             "G": 0.35,
                             "T": 0.30}
                             },
             "prefixes_prob": {""},
             "order": 0}

    proba, log_proba = sequence_probability("ATGCCG", model)

    expected = 0.20 * 0.30 * 0.35 * 0.15 * 0.15 * 0.35

    assert proba == pytest.approx(expected)
    assert log_proba == pytest.approx(math.log10(expected))

def test_sequence_probability_order_1():
    model = {"matrix": {"A": {"A": 0.2, "C": 0.1, "G": 0.4, "T": 0.4},
                        "T": {"G": 0.2},
                        "G": {"C": 0.2},
                        "C": {"A": 0.1, "C": 0.6, "G": 0.2, "T": 0.1}},
             "prefixes_prob": {"A": 0.5, "C": 0.2, "G": 0.1, "T": 0.2},
             "order": 1}

    proba, log_proba = sequence_probability("ATGCCG", model)

    expected = 0.5 * 0.4 * 0.2 * 0.2 * 0.6 * 0.2

    assert proba == pytest.approx(expected)
    assert log_proba == pytest.approx(math.log10(expected))

def test_sequence_probability_order_2():
    model = {"matrix": {"AT": {"A": 0.2, "C": 0.2, "G": 0.4, "T": 0.2},
                        "TG": {"A": 0.4, "C": 0.2, "G": 0.2, "T": 0.2},
                        "GC": {"C": 0.5},
                        "CC": {"G": 0.2}},
             "prefixes_prob": {"AA": 0.02, "AT": 0.03, "CA": 0.02},
#             "prefixes_prob": {"AT": 0.3},
             "order": 2}

    proba, log_proba = sequence_probability("ATGCCG", model)

    expected = 0.03 * 0.4 * 0.2 * 0.5 * 0.2

    assert proba == pytest.approx(expected)
    assert log_proba == pytest.approx(math.log10(expected))

def tests_sequence_probability_empty_seq_order_0():
    model = {"matrix": {"": {"A": 0.20,
                             "C": 0.15,
                             "G": 0.35,
                             "T": 0.30}
                             },
             "prefixes_prob": {""},
             "order": 0}

    proba, log_proba = sequence_probability("", model)

    assert proba == 0
    assert log_proba == float("-inf")

def tests_sequence_probability_empty_seq_order_1():
    model = {"matrix": {"A": {"A": 0.2, "C": 0.1, "G": 0.4, "T": 0.4},
                        "T": {"G": 0.2},
                        "G": {"C": 0.2},
                        "C": {"A": 0.1, "C": 0.6, "G": 0.2, "T": 0.1}},
             "prefixes_prob": {"A": 0.5, "C": 0.2, "G": 0.1, "T": 0.2},
             "order": 1}

    proba, log_proba = sequence_probability("", model)

    assert proba == 0
    assert log_proba == float("-inf")

def tests_sequence_probability_empty_seq_order_2():
    model = {"matrix": {"AT": {"A": 0.2, "C": 0.2, "G": 0.4, "T": 0.2},
                        "TG": {"A": 0.4, "C": 0.2, "G": 0.2, "T": 0.2},
                        "GC": {"C": 0.5},
                        "CC": {"G": 0.2}},
             "prefixes_prob": {"AA": 0.02, "AT": 0.03, "CA": 0.02},
             "order": 2}


    proba, log_proba = sequence_probability("", model)

    assert proba == 0
    assert log_proba == float("-inf")

def test_sequence_probability_prefix_only():
    model = {"matrix": {"AT": {"A": 0.2, "C": 0.2, "G": 0.4, "T": 0.2},
                        "TG": {"A": 0.4, "C": 0.2, "G": 0.2, "T": 0.2},
                        "GC": {"C": 0.5},
                        "CC": {"G": 0.2}},
             "prefixes_prob": {"AA": 0.02, "AT": 0.03, "CA": 0.02},
             "order": 2}


    proba, log_proba = sequence_probability("AT", model)

    assert proba == pytest.approx(0.03)
    assert log_proba == math.log10(0.03)

def tests_sequence_probability_seq_to_shorter_than_order():
    model = {"matrix": {},
             "prefixes_prob": {},
             "order": 3}

    proba, log_proba = sequence_probability("AT", model)

    assert proba == 0.0
    assert log_proba == float("-inf")
