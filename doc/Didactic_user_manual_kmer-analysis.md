# Tutorial  — Using `kmer-analysis`

RISCH Anouk  
2026-04-23

## Table of contents 

- [Introduction](#introduction)
- [Input format](#input-format)
- [Output format](#output-format)
- [1. Basic usage: Counts Occurrences](#1-basic-usage-counts-occurrences) 
- [2. How k-mers are Counted](#2-how-k-mers-are-counted)
- [3. Computing frequencies](#3-computing-frequencies)
- [4. Exploring Results as TSV](#4-exploring-result-as-tsv)
- [5. Statistical Significance](#5-statistical-significance)
- [6. Background Models](#6-background-models)
- [7. Oligonucleotide size](#7-oligonucleotide-size)
- [8.Reverse Complements - Strand option](#8-reverse-complements---strand-options)
- [9. Understanding Output Columns](#9-understanding-output-columns)
- [10. Getting Help](#10-getting-help)
- [11. Pydoc](#11-pydoc)
- [12. Probabilities](#12-probabilities)

### Introduction

`kmer-analysis` computes statistics about k-mer occurrences in DNA sequences.

It counts the occurrences of all k-mers in a set of input sequences (provided as a fasta-formatted file),
and derives different statistics:

- observed k-mer occurrence counts and relative frequencies
- expected occurrences and frequencies based on either a Bernoulli or a Markov model
- over- or under-representation statistics (P-value, E-value, significance)

The computed statistics are adapted according to user-selected return fields. 

Over-representation statistics rely on the comparison between the k-mer occurrences observed in the input sequence set 
and those expected under a specified background model (Bernoulli or Markov). This background model is typically 
estimated from a set of background sequences. 

### Input format

The program takes as input a fasta-formatted sequence files, which can contains one or several sequences. 

### Output format

The output is a tabular text file (tab-separated values, recommended extension `.tsv`) with one row per k-mer
and one column per statistics. 

The column content is indicated in the header of the output (row prefixed with a `#` character). 

| Column   | Description                  |
|:---------|:-----------------------------|
| seq      | oligomer sequence            |
| id       | oligomer identifier          |
| exp_freq | expected relative frequency  |
| obs_freq | 	observed relative frequency |
| occ      | observed occurrences         |
| exp_occ  | expected occurrences         |

Depending on the specified level of verbosity, the output can contain additional information in comment lines, 
prefixed by a semicolumn character (`;`)

### 1. Basic usage: Counts Occurrences

The simplest use of the program is to count the number of occurrences of each k-mer in the input sequences.

```commandline
python3 scripts/kmer_analysis.py -i sequence.fa -k 2
```

#### Parameters

| Option | Meaning              |
|--------|----------------------|
| `-i`   | Input Fasta file     |
| `-k`   | Length of the k-mers |


#### Output example

```
#seq    occ
AAA     12
AAC     4
AAG     7
```

Interpretation: 

- `AAA` appears 12 times in the sequence
- `AAC` appears 4 times

### 2. How K-mers are counted

The tool scans each sequence using a sliding window of size $k$ and counts all valid k-mers. 
The compute time thus increases linearly with the total sequence size (i.e. the sum of sequence sizes in the 
input sequence file). 

#### Example

For `ATTCG` sequence  with `K=2`.

|  Position  | k-mer  |
|:----------:|:------:|
|     0      |  ATT   |
|     1      |  TTC   |
|     2      |  TCG   |

The total numbers of k-mer positions depends on the k-mer size and the sequence lengths    

$N = \sum_{i=1}^{S} (L_i - k + 1)$

Where $N$ is the number of positions where a k-mer can be found in the input sequences. 

where 
- $S$ is the number of sequence in the input file
- $L_i$ is the length of the $i^{th}$ sequence
- $k$ is the k-mer size


### 3. Computing relative frequencies

The option `-return freq` prints the observed relative frequencies, in addition to the observed occurrences. 

**Output Example**
```
#seq    occ         obs_freq
AAA     occ=12      freq=0.031
AAC     occ=2       freq=0.010
...    
```

#### Interpretation

$F(K) = \frac{N(K)}{\sum_{i=1}^{W}N(i)}$

Where
- $K$ is a given k-mer of size $k$
- $N(K)$ is the number of occurrences of $K$
- $F(K)$ is the relative frequency
- $W$ is the number of k-mers (words) of size $k$ found in the sequence

### 4. Exploring Result as TSV

To save results in tabular format :

```commandline
python3 scripts/kmer_analysis.py -i sequences.fa -k 2 -o results.tsv
```
TSV files can be opened in :  
- Excel
- LibreOffice Calc
- R 
- Notepad
- Others

### 5. Statistical Significance
Some k-mers appear often simply by chance. 

The program can compare : 
- Observed occurrences
- Excepted occurrences

**Example Output**  
```
kmer    occ     exp_occ
TAAA    18      4.2
```
**Interpretation**  

`TAAA` appears much more often than expected.
This may indicate a biologically meaningful motif.

### 6. Background Models
Background models are probabilistic models used to estimate the expected frequency of each oligonucleotide (short DNA motif → k-mer). The choice of the expected frequency can strongly influence the results of the analysis.

**Predefined background frequencies :**  This approach consists of comparing the oligonucleotide frequencies observed in the query sequence with those measured in a reference sequence used as the background model.

**Bernoulli Model**  
Assuming nucleotides are independent (`P(A)`, `P(C)`, `P(T)`, `P(G)`) :
```
P(A) = P(T) = P(C) = P(G) 
```
Simple and fast but not realistic 


**Markov Model**  
Accounts for local dependencies.  
Expected oligonucleotide (k-mer) frequencies are estimated from the observed frequencies of subwords within the input sequence set. This approach accounts for higher-order dependencies between adjacent residues, thereby incorporating local sequence context into the expected frequency model.

**Example** 
```
P(GATAAC) = P(GATAA) * P(C|GATAA)
          = P(GATAA) * P(ATAAC) / P(ATAA)
```
Thus 
```
Expected(GATAAC) = observed(GATAA) * observed(ATAAC) / observed(ATAA) 
```

More realistic for biological sequences.



A 0 order Markov model uses the observed single-residue frequencies to estimate expected oligonucleotide frequencies, assuming no dependencies between adjacent residues. Under this assumption, sequence positions are treated as independent and identically distributed, corresponding to a Bernoulli model.
 
### 7. Oligonucleotide size

The analysis can be performed using oligonucleotides of length ranging from 1 to 10. Setting the size to 1 corresponds to quantifying the usage of the underlying alphabet within the input sequences. For the identification of putative regulatory sites, we recommend initiating the analysis with hexanucleotides (size = 6) and exploring k-mer lengths between 4 and 8. Biologically relevant patterns that are significantly overrepresented are typically consistently detected across analyses performed with different k-mer sizes.

**Command**
```commandline
python3 scripts/kmer_analysis.py -i sequence.fa -k lengh_kmers
```
| Option | Meaning              |
|--------|----------------------|
| `-k`   | Length of the k-mers |

**Small k (1–3)**

Advantages :  
- fast
- robust statistics

Limitations :  
- low specificity

**Medium k (4–8)**

- Good compromise.

**Large k (9+)**

Advantages :
- specific motifs

Limitations :
- many possible k-mers
- sparse counts
- slower computation

### 8. Reverse Complements - Strand Options
**Command** 
```commandline
python3 scripts/kmer_analysis.py -i sequences.fa -k 2 -o results.tsv -s both
```
| Option    | Meaning                       |
|:----------|:------------------------------|
| `single`  | Forward                       |
| `both`    | Forward + reverse complement  |

By selecting the `both` option, the occurrences of each oligonucleotide are aggregated across both strands. This enables the detection of sequence elements that function in an orientation-independent manner. 
### 9. Understanding Output Columns

| Columns                          | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                               |
|:---------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| seq                              | oligomer sequence                                                                                                                                                                                                                                                                                                                                                                                                                     |
| id                               | oligomer identifier                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Expected frequencies (exp_freq)  | Proportion of a given oligonucleotide that is expected within the sequence set under the selected probabilistic model. It is computed from the expected occurrences by normalizing them with respect to the total number of possible oligonucleotide positions in the dataset. This measure reflects the theoretical probability of observing each oligonucleotide given the assumed model of sequence composition and dependencies.  |
| Frequencies (obs_freq)           | 	Relative frequencies are defined as the number of occurrences of each oligonucleotide divided by the total number of occurrences across all oligonucleotides in the dataset.                                                                                                                                                                                                                                                         |
| Occurrences (occ)                | Number of occurrences of each oligonucleotide observed in the dataset. All matches are included, and overlapping occurrences are detected and summed during the count.                                                                                                                                                                                                                                                                |
| Expected occurrences (exp_occ)   | Predicted number of times a given oligonucleotide is expected to appear within the sequence set. This value is computed according to the probabilistic model chosen by the user (see above), which defines how sequence composition and dependencies are taken into account.                                                                                                                                                          |

### 10. Getting Help
List all options :
````commandline
python3 /scripts/kmer_analysis.py --help
````
### 11. Pydoc
Automatic Python code documentation : 
- function 
- parameters 

**Command**
````commandline
pydoc /scripts/kmer_analysis.py
````
### 12. Probabilities

EXPECTED OCCURRENCES  

${Exp\_occ} = p \times T = p \times \sum_{j=1}^{S} (L_j - k + 1)$
	
Where	
- $p$ is the probability of the pattern
- $S$ is the number of sequences in the sequence set. 
- $L_j$ is the length of the jth regulatory region
- $k$ is the length of oligomer
- $T$ is the number of possible matching positions.

NUMBER OF POSSIBLE POSITIONS  

$T = j = \sum_{j=1}^{s} (L_j -k + 1)$


PROBABILITY OF A SEQUENCE SEGMENT 

#### Bernoulli model

$p = \prod_{i=1}^{k} P(r_i)$ 


Where   
- $i = 1..k$ is the index of nucleotide positions
- $ri$ is the residue found at position I   
- $P(ri)$ is the probability of this residue

#### Markov  model

$P(r_i \mid S_{1-m,i-1}) $

where   
- $S$ is the sequence
- $i$ is the current position 
- $m$ is the Markov model order
- $S_{i−m,1−i}$ = m nucleotides previously

Sequence probability given the background mode

$p = P(r_i) \times i = \prod_{i=2}^{k} P(r_i \mid S_{1-m,1-1})$

          