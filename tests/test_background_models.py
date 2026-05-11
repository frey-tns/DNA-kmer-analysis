#################
#   Libraries   #
#################
import pytest

############################
#   Internal libraries     #
############################
import bioseq_kmers.background_models as bg

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