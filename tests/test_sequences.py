"""
Unit tests for the functions defined in bioseq_kmers/sequences.py

USAGE

pytest -v tests/test_sequences.py

"""

# External libraries
import pytest

# Internal libraries
from bioseq_kmers.sequences import read_fasta

#######################
#   Functional test   #
#######################
def test_read_fasta(tmp_path):
    """Check outputs for valid input sequences"""
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

def test_read_fasta_with_spaces(tmp_path):
    """Check that spaces are suppressed when reading fasta file"""
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nACGGNA \n"
                     ">seq2\n ACGTACGN\n"
                     ">seq3\nA A A A A A\n")

    sequences, total_length, seq_number = read_fasta(fasta)

    assert sequences == {"seq1": "ACGGNA",
                         "seq2": "ACGTACGN",
                         "seq3": "AAAAAA"}
    assert total_length == 20
    assert seq_number == 3


def test_read_fasta_empty_file(tmp_path):
    """Return empty outputs for an empty FASTA file"""
    fasta = tmp_path / "test.fasta"
    fasta.write_text("")

    sequences, total_length, seq_number = read_fasta(fasta)

    assert sequences == {}
    assert total_length == 0
    assert seq_number == 0

def test_read_fasta_no_header(tmp_path):
    """Return empty outputs if no header is present."""
    fasta = tmp_path / "test.fasta"
    fasta.write_text("ATGCTAGCATG\nATGGTGATG\n")

    sequences, total_length, seq_number = read_fasta(fasta)

    assert sequences == {}
    assert total_length == 0
    assert seq_number == 0

def test_read_fasta_invalid_character(tmp_path):
    """Reject sequences containing non-DNA characters."""
    fasta_file = tmp_path / "bad.fasta"
    fasta_file.write_text(">seq1\nACGRXA\nAGTCNT")

    with pytest.raises(ValueError):
        read_fasta(fasta_file)
