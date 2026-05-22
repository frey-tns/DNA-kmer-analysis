#################
#   Libraries   #
#################
from collections import Counter
import pytest
import math

############################
#   Internal libraries     #
############################
import bioseq_kmers.kmer_stats as stat

####################
#   Function test  #
####################
def test_reverse_complement_basic():
    assert stat.reverse_complementary("ACTG") == "CAGT"

def test_reverse_complement_invalid_character():
    with pytest.raises(KeyError):
        stat.reverse_complementary("ACTX")

def test_canonic_kmer_reverse_complement():
    assert stat.canonic_kmer("ATCG") == stat.canonic_kmer("CGAT")

def test_count_kmer_single_strand():
    sequence = {"seq1":"ATAT"}

    observed = stat.counts_kmer(sequence, 2, "single")

    expected = Counter({"AT":2,
                        "TA":1})

    assert observed == expected


def test_count_kmer_both_strand():
    sequence = {"seq1": "ATGC"}

    observed = stat.counts_kmer(sequence, 2, "both")

    expected = Counter({"AT": 1,
                        "CA": 1,
                        "GC": 1})

    assert observed == expected

def test_counts_kmer_ignore_n():
    sequences = {"seq1": "ATNA"}

    observed = stat.counts_kmer(sequences, 2, "single")

    expected = Counter({"AT": 1,})

    assert observed == expected

def test_counts_kmer_sequence_shorter_than_k():
    sequences = {"seq1": "AT"}

    observed = stat.counts_kmer(sequences, 3, "single")

    assert observed == Counter()

def test_nucleotide_frequencies():
    sequences = {"seq1": "AATT",
                 "seq2": "CCGG"}

    frequencies = stat.nucleotide_frequencies(sequences)

    assert frequencies == {"A": 0.25,
                           "T": 0.25,
                           "G": 0.25,
                           "C": 0.25}

def test_nucleotide_frequencies_ignor_n():
    sequences = {"seq1": "AATN"}

    frequencies = stat.nucleotide_frequencies(sequences)

    assert frequencies == {"A": 2/3,
                           "T": 1/3,
                           "G": 0.0,
                           "C": 0.0}

def test_nucleotide_frequencies_sum_to_one():
    sequences = {"seq1": "AATCG"}

    frequencies = stat.nucleotide_frequencies(sequences)

    assert sum(frequencies.values()) == pytest.approx(1.0)

@pytest.mark.parametrize("lam", [10, 20, 50, 3000])
def test_poisson_statistics_lambda(lam):
    occ = 10
    nb_tests = 1000

    stats = stat.poisson_statistics(occ, lam, nb_tests)

    assert math.isclose(stats["occ_E"], stats["occ_P"] * nb_tests)
    assert round(stats["occ_sig"], 3) == round(-math.log10(stats["occ_E"]), 3)


@pytest.mark.parametrize(("occ", "occ_sig"),
                         [(10, 3.953),
                          (50, 61.909),
                          (100, 155.400),
                          (150, 260.188),
                          (160, 282.105),
                          (170, 304.293),
                          (180, 320)])
def test_poisson_statistics_occ_series(occ, occ_sig):
    lam = 1
    nb_tests = 1000

    stats = stat.poisson_statistics(occ, lam, nb_tests)

    assert round(stats["occ_sig"], 3) == occ_sig

