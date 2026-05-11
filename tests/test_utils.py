#################
#   Libraries   #
#################
import shlex
import os
import pytest

############################
#   Internal libraries     #
############################
import bioseq_kmers.utils as utils

####################
#   Function test  #
####################
def test_format_command_line_relative_paths(tmp_path, monkeypatch):
    project = tmp_path / "project"
    project.mkdir()

    monkeypatch.chdir(project)

    input_rel = os.path.join("data", "input.fasta")
    output_rel = os.path.join("results", "out.tsv")

    argv = ["python",
            str(project/"script.py"),
            "-i",
            str(project/"data"/"input.fasta"),
            "-o",
            str(project/"results"/"out.tsv")]

    command = utils.format_command_line(argv)

    assert f"-i {shlex.quote(input_rel)}" in command
    assert f"-o {shlex.quote(output_rel)}" in command