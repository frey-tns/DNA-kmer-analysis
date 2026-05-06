# kmer-analysis

RISCH Anouk  
2026-05-06

## Table of contents

- [Description](#-description-)
- [Pipeline workflow](#-pipeline-workflow)
- [Input](#-input)
- [Output](#-output)
- [Requirements](#-requirements)
- [Installation](#-installation)
- [Quick start](#-quick-start)
- [Usage](#-usage-)
- [Examples](#examples)
- [Related tools](#-related-tools)

### 🧬 Description 

`kmer-analysis` is a command-line tool for analyzing k-mers in biological sequences.
It allows downstream analyses such as k-mers statistics, background model estimation (Bernoulli/Markov models) and sequence probability computation.  
This bioinformatics software identifies enriched or depleted k-mers in DNA sequence sets by comparing observed counts with expected frequencies under background probabilistic models.

This software is a mimic of the `oligo-analysis` tool from [RSAT](https://rsat.eead.csic.es/plants/oligo-analysis_form.cgi).

### 🔄 Pipeline workflow
The analysis consists of the following steps:

1. Read input DNA sequences from a FASTA file
2. Compute sequence-level statistics:
   - number of sequences
   - total sequence length
3. Count observed k-mer occurrences
4. Compute observed k-mer frequencies
5. Estimate expected frequencies using a background model:
   - Bernoulli model
   - Markov model
6. Compute enrichment statistics:
   - observed/expected ratio
   - P-value
   - E-value
   
### 📂 Input

- FASTA files 
- Sequences are expected to contain standard nucleotide symbols (`A`, `C`, `G`, `T`, `N`)

### 📈 Output

The program generates a tab-separated value file (extension  `.tsv`) with one row per observed k-mer in the input sequences,
and one column per k-mer statistics.

Available output fields include:
- `seq` — k-mer sequence 
- `id` — Strand mode identifier (forward + reverse) 
- `occ` — Observed occurrences 
- `exp_occ` — Expected occurrences 
- `obs_freq` — Observed frequency (frequency/position)
- `exp_freq` — Expected frequency under the selected background models (Bernoulli or Markov) 
- `ratio` — Ratio of observed to expected occurrences 
- `occ_P` — P-value 
- `occ_E` — E-value 
- `sig` — significance score: -log10(E-value)

Example output:
```
#seq	id	    exp_freq        	    obs_freq    	    occ	        exp_occ
TAA	TAA|TTA	    0.03224340262535573	    0.049433682650576904    467 	304.60
GTA	GTA|TAC	    0.01752288814005347	    0.03186196676193501	    301 	165.54
ACA	ACA|TGT	    0.019012088273685668    0.03916587276384037	    370         179.61
CAA	CAA|TTG	    0.019012088273685668    0.03916587276384037	    370 	179.61
```
⚠️ The output columns depend on the fields selected with `--return`.

### ⚙️ Requirements

The requirements are described in the files:

- [environment.txt](environment.txt): user requirements
- [environment-dev.txt](environment-dev.txt): additional requirements for developers

### 📎 Installation

```bash
pip install .
```
### 🚀 Quick start

```bash
kmer-analysis -i data/example.fasta -k 6 -s both -o results/output.tsv
```

### ▶️ Usage 

SYNOPSIS USAGE
- Display command synopsis:
```bash
kmer-analysis
```
- Print help message:
```bash
kmer-analysis --help
```
- usage:

```bash
kmer-analysis [-h] --input FASTA_FILE --kmer KMER_LENGTH --output OUTPUT_FILE [options]
```

💡 Options

    -h, --help
        Display this help message and exit.

    -i, --input FASTA_FILE
        Input sequence file in fasta format.

    -k, --kmer KMER_LENGTH
        Length of k-mer sequences to analyze (integer ≥ 1).

    -s, --strand {single,both}
        Strand mode to compute counts for.
        Handling of the strands for computing the occurrences and derived statistics.

        single: count occurrences on the forward strand only
        both: regroup k-mers by pairs of reverse complements
            This actually makes the counts and derived statistics strand-insensitive

    -o, --output OUTPUT_FILE
        Output file in tab-separated values, with recommended extension .tsv

    -r, --return FIELD[,FIELD...]
        Comma-separated list of output fields.



Available values for `--return`:

- `occ`       : observed occurrences
- `obs_freq`  : observed relative frequency
- `exp_occ`   : expected occurrences
- `exp_freq`  : expected frequency


- proba     : statistical significance measures, including:  
  - `p_value` : binomial P-value according to the specified background model
  - `e_value` : E-value, i.e. P-value × T, where T is the number of k-mers tested
  - `sig`     : significance score, defined as -log10(E-value)


### Examples

```bash
kmer-analysis -i data/yeast_MET_upstream.fasta  -k 6 -s both -o results/yeast_MET_upstream_6nt_2str.tsv

kmer-analysis -i data/yeast_MET_upstream.fasta  -k 6 -s both \
    --return occ,obs_freq,exp_occ,exp_freq \
    -o results/yeast_MET_upstream_6nt_2str.tsv
```

### 🧩 Related tools


This program is part of a small toolkit for sequence analysis.

- `markov-from-kmers` — estimate background Markov models from k-mer frequencies (from RSAT: [convert-background-model](https://rsat.eead.csic.es/plants/convert-background-model_form.cgi) )
- `seq-proba` — compute sequence probabilities under Bernoulli or Markov models (from RSAT: [seq-proba](https://rsat.eead.csic.es/plants/seq-proba_form.cgi))