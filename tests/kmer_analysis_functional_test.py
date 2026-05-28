import subprocess
from pathlib import Path

import pytest


@pytest.mark.parametrize("markov_order", range(0, 6))
def test_kmer_analysis_markov_from_bgfiles(markov_order):
    # ------------------------------------------------------------------
    # Input / output files
    # ------------------------------------------------------------------

    output_dir = Path("results/kmer-analysis-markov-series/bg-from-file")
    output_dir.mkdir(parents=True, exist_ok=True)

    input_fasta = "data/seq/yeast_MET_upstream.fasta"

    background_model = (
        f"data/bg-models/"
        f"yeast_all-upstream-noorf_mkv{markov_order}.tsv"
    )

    output_tsv = (
            output_dir /
            f"yeast_MET_upstream_k6_mkv{markov_order}_enriched.tsv"
    )

    filtered_tsv = (
            output_dir /
            f"yeast_MET_upstream_k6_mkv{markov_order}_enriched_sig0.tsv"
    )

    # ------------------------------------------------------------------
    # Command 1 : kmer-analysis
    # ------------------------------------------------------------------

    cmd1 = [
        "kmer-analysis",
        "-i", input_fasta,
        "-k", "6",
        "-s", "both",
        "--background", background_model,
        "--return",
        "occ,exp_occ,obs_freq,exp_freq,occ_P,occ_E,occ_sig",
        "--sort", "alpha",
        "-o", str(output_tsv),
    ]

    result1 = subprocess.run(
        cmd1,
        capture_output=True,
        text=True,
    )

    print(result1.stdout)
    print(result1.stderr)

    assert result1.returncode == 0
    assert output_tsv.exists()

    # ------------------------------------------------------------------
    # Command 2 : awk + sort
    # ------------------------------------------------------------------

    shell_cmd = (
        f"awk '$9 > 0' {output_tsv} "
        f"| sort -n -k 9 "
        f"> {filtered_tsv}"
    )

    result2 = subprocess.run(
        shell_cmd,
        shell=True,
        capture_output=True,
        text=True,
    )

    print(result2.stdout)
    print(result2.stderr)

    assert result2.returncode == 0
    assert filtered_tsv.exists()

    # ------------------------------------------------------------------
    # Optional sanity checks
    # ------------------------------------------------------------------

    assert filtered_tsv.stat().st_size > 0

@pytest.mark.parametrize("markov_order", range(1, 5))
def test_kmer_analysis_markov_from_input(markov_order):
    # ------------------------------------------------------------------
    # Input / output files
    # ------------------------------------------------------------------

    output_dir = Path("results/kmer-analysis-markov-series/bg-from-input")
    output_dir.mkdir(parents=True, exist_ok=True)

    input_fasta = "data/seq/yeast_MET_upstream.fasta"


    output_tsv = (
        output_dir /
        f"yeast_MET_upstream_k6_mkv{markov_order}_bginput_enriched.tsv"
    )

    filtered_tsv = (
        output_dir /
        f"yeast_MET_upstream_k6_mkv{markov_order}_bginput_enriched_sig0.tsv"
    )

    # ------------------------------------------------------------------
    # Command 1 : kmer-analysis
    # ------------------------------------------------------------------

    cmd1 = [
        "kmer-analysis",
        "-i", input_fasta,
        "-k", "6",
        "-s", "both",
        "--markov-order", str(markov_order),
        "--return",
        "occ,exp_occ,obs_freq,exp_freq,occ_P,occ_E,occ_sig",
        "--sort", "alpha",
        "-o", str(output_tsv),
    ]

    print(cmd1)

    result1 = subprocess.run(
        cmd1,
        capture_output=True,
        text=True,
    )

    print(result1.stdout)
    print(result1.stderr)

    assert result1.returncode == 0
    assert output_tsv.exists()

    # ------------------------------------------------------------------
    # Command 2 : awk + sort
    # ------------------------------------------------------------------

    shell_cmd = (
        f"awk '$9 > 0' {output_tsv} "
        f"| sort -n -k 9 "
        f"> {filtered_tsv}"
    )

    result2 = subprocess.run(
        shell_cmd,
        shell=True,
        capture_output=True,
        text=True,
    )

    print(result2.stdout)
    print(result2.stderr)

    assert result2.returncode == 0
    assert filtered_tsv.exists()

    # ------------------------------------------------------------------
    # Optional sanity checks
    # ------------------------------------------------------------------

    assert filtered_tsv.stat().st_size > 0


@pytest.mark.parametrize("markov_order", range(1, 5))
def test_oligo_analysis_markov_from_input(markov_order):
    # ------------------------------------------------------------------
    # Input / output files
    # ------------------------------------------------------------------


    output_dir = Path("results/oligo-analysis-markov-series/bg-from-input")
    output_dir.mkdir(parents=True, exist_ok=True)

    input_fasta = "data/seq/yeast_MET_upstream.fasta"


    output_tsv = (
        output_dir /
        f"oligos_yeast_MET_upstream_k6_mkv{markov_order}_bginput_enriched.tsv"
    )

    # ------------------------------------------------------------------
    # Command 1 : oligo-analysis
    # ------------------------------------------------------------------

    cmd1 = [
        "source ~/packages/rsat/RSAT_config.bashrc && "
        "rsat oligo-analysis -v 1"
        f" -i {input_fasta}"
        " -format fasta"
        " -l 6"
        " -2str"
        f" -markov {markov_order}"
        " -return freq,occ,proba,rank"
        " -seqtype dna"
        " -pseudo 0"
        f" -o {output_tsv}"
    ]

    print("\nCommand:\n" + " ".join(cmd1))

    result1 = subprocess.run(
        cmd1,
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )

    print(result1.stdout)
    print(result1.stderr)

#    assert result1.returncode == 0
    assert output_tsv.exists()

@pytest.mark.parametrize("markov_order", range(1, 5))
def test_oligo_analysis_markov_from_file(markov_order):
    # ------------------------------------------------------------------
    # Input / output files
    # ------------------------------------------------------------------


    output_dir = Path("results/oligo-analysis-markov-series/bg-from-file")
    output_dir.mkdir(parents=True, exist_ok=True)

    input_fasta = "data/seq/yeast_MET_upstream.fasta"

    background_model = (
        f"data/bg-models/"
        f"yeast_all-upstream-noorf_mkv{markov_order}.tsv"
    )

    output_tsv = (
        output_dir /
        f"oligos_yeast_MET_upstream_k6_mkv{markov_order}_bginput_enriched.tsv"
    )

    # ------------------------------------------------------------------
    # Command 1 : oligo-analysis
    # ------------------------------------------------------------------

    cmd1 = [
        "source ~/packages/rsat/RSAT_config.bashrc && "
        "rsat oligo-analysis -v 1"
        f" -i {input_fasta}"
        " -format fasta"
        " -l 6"
        " -2str"
        f" -markov {markov_order}"
        " -return freq,occ,proba,rank"
        " -seqtype dna"
        " -pseudo 0"
        f" -o {output_tsv}"
    ]

    print("\nCommand:\n" + " ".join(cmd1))

    result1 = subprocess.run(
        cmd1,
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )

    print(result1.stdout)
    print(result1.stderr)

#    assert result1.returncode == 0
    assert output_tsv.exists()
