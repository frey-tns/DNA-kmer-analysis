import subprocess
from pathlib import Path

import pytest


@pytest.mark.parametrize("markov_order", range(0, 8))
def test_markov_from_seq(markov_order):
    # ------------------------------------------------------------------
    # Input / output files
    # ------------------------------------------------------------------

    output_dir = Path("results/bg-models")
    output_dir.mkdir(parents=True, exist_ok=True)

    input_fasta = "data/seq/yeast_all-upstream-noorf.fasta.gz"

#    background_model = (
#        f"data/bg-models/"
#        f"yeast_all-upstream-noorf_Markov_mkv{markov_order}.tsv"
#    )

    output_tsv = (
            output_dir /
            f"yeast_all-upstream-noorf_mkv{markov_order}.tsv"
    )

    # ------------------------------------------------------------------
    # Command : markov-from-seq
    # ------------------------------------------------------------------

    cmd = [
        "markov-from-seq",
        "-i", input_fasta,
        "-m", str(markov_order),
        "-o", str(output_tsv),
    ]

    result = subprocess.run(cmd,
                            capture_output=True,
                            text=True)

    print(result.stdout)
    print(result.stderr)

    assert result.returncode == 0
    assert output_tsv.exists()
