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
- [Usage and options](#-usage-and-options)

### 🧬 Description 

`kmer-analysis` is a command-line tool for analyzing k-mers in biological sequences.
It allows downstream analyses such as k-mers statistics, background model estimation (Bernoulli/Markov models) and sequence probability computation.  
This bioinformatics software identifies enriched or depleted k-mers in DNA sequence sets by comparing observed counts with expected frequencies under background probabilistic models.

This software is a re-implementation of the `oligo-analysis` tool from [RSAT](https://rsat.eead.csic.es/plants/oligo-analysis_form.cgi).

### Package contents

| Tool | Purpose                                                                                                                                                                                                               |
|:----------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `markov-from-seq` | estimate the parameters of a Markov model of order $m$ (in the form of a transition matrix) from a set of background sequences. By extension, a Bernoulli model can be computed by setting $m=0$.                     |
| `markov-from-kmers` | estimate a background Markov model of order $m$ from a table of k-mer occurrences, where $k=m+1$ (corresponds to RSAT [convert-background-model](https://rsat.eead.csic.es/plants/convert-background-model_form.cgi)) |
| `kmer-analysis` | compute k-mer occurrences in a set of input sequences, and derive over- or under-representation statistics given a user-specified background model                                                                    |
| `seq-proba` | compute sequence probabilities under a Bernoulli or a Markov model (corresponds to RSAT [seq-proba](https://rsat.eead.csic.es/plants/seq-proba_form.cgi))                                                             |


### 🔄 Pipeline workflow

1. `markov-from-seq`: compute a background model based on a set of background sequences. 
2. `kmer-anlysis`: compute  occurrences and derived statistics of the sequences of interest, given the background model. 

### `kmer-analysis` steps

The analysis consists of the following steps:

1. Read input DNA sequences from a FASTA file
2. Compute sequence-level statistics:
   - number of sequences
   - total sequence length
3. Count observed k-mer occurrences
4. Compute observed k-mer relative frequencies
5. Estimate expected frequencies using a background model (Markov or Bernoulli)
6. Compute enrichment/depletion statistics:
   - observed/expected ratio
   - P-value
   - E-value
   - significance

### 📂 Input

- A FASTA-formated file containing one or more sequences. The statistics are computed on the sequence set as a whole and not on each sequence separately.
- Sequences are expected to contain standard nucleotide symbols (`A`, `C`, `G`, `T`, `N`)

### 📈 Output

The program `kmer-analysis` generates a tab-separated value file (extension  `.tsv`) with one row per observed k-mer in the input sequences,
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

Output example:

```
#seq	id	    exp_freq        	    obs_freq    	    occ	        exp_occ
TAA	TAA|TTA	    0.03224340262535573	    0.049433682650576904    467 	304.60
GTA	GTA|TAC	    0.01752288814005347	    0.03186196676193501	    301 	165.54
ACA	ACA|TGT	    0.019012088273685668    0.03916587276384037	    370         179.61
CAA	CAA|TTG	    0.019012088273685668    0.03916587276384037	    370 	179.61
```
⚠️ The output columns depend on the fields selected with  the option `--return`.

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
## Compute the background model from the set of all non-coding upstream sequences from the yeast Saccharomyces cerevisiae
mkdir -p results/bg-models
markov-from-seq -i data/seq/yeast_all-upstream-noorf.fasta -m 4 -o results/bg-models/yeast_allup-noorf_m4.tsv

## Analyse kmers in the upstream sequences of yeast genes involved in methionine metabolism
mkdir -p results/kmer-statistics
kmer-analysis -i data/seq/yeast_MET_upstream.fasta -k 6 -s both -o results/kmer-statistics/yeast_MET_upstream_6nt_2str.tsv
```

### ▶️ Usage and options

When called without any argument, each tool of the package displays a usage line with the supported arguments

```bash
kmer-analysis
```

A more detailed description of the tool and the options can be obtained with the option `--help`. 

```bash
kmer-analysis --help
```

