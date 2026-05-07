#################
#   Libraries   #
#################
import pytest

############################
#   Internal libraries     #
############################
from bioseq_kmers.sequences import read_fasta

#######################
#   Functional test   #
#######################
def test_read_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nACGGNA\n"
                     ">seq2\nACGTACGN\n"
                     ">seq3\nAAAAAA\n")

    sequences, total_length, seq_number = read_fasta(fasta)

    assert sequences == {"seq1": "ACGGNA",
                         "seq2": "ACGTACGN",
                         "seq3": "AAAAAA"}
    assert total_length == 20
    assert seq_number == 3

def test_read_fasta_invalid_character(tmp_path):
    fasta = tmp_path / "bad.fasta"
    fasta.write_text(">seq1\nACGRXA\nAGTCNT")

    with pytest.raises(ValueError):
        read_fasta(fasta)