import pytest
from bioseq_kmers.sequences import read_fasta

def test_read_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nACGTACGN\n"
                     ">seq2\nACGTACGN\n"
                     ">seq3\nACGTACGN\n")

    sequences, total_length, seq_number = read_fasta(fasta)

    assert sequences == {"seq1": "ACGTACGN",
                         "seq2": "ACGTACGN",
                         "seq3": "ACGTACGN"}
    assert total_length == 24
    assert seq_number == 3