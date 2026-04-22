# Tutorial  — Using `kmer-analysis`

RISCH Anouk  
2026-04-2026

## Table of contents 

- [Introduction](#introduction-)
- [Input format](#input-format-)
- [Output format](#output-format-)
- [1. Basic usage: Counts Occurrences](#1-basic-usage-counts-occurrences-) 
- [2. How k-mers are Counted](#2-how-k-mers-are-counted)
- [3. Computing frequencies](#3-computing-frequencies-)
- [4. Exploring Results as TSV](#4-exploring-result-as-tsv-)
- [5. Statistical Significance](#5-statistical-significance)
- [6. Background Models](#6-background-models-)
- [7. Oligonucleotide size](#7-oligonucleotide-size)
- [8.Reverse Complement Strand option](#8-reverse-complement---strand-options-)
- [9. Understanding Output Columns](#9-understanding-output-columns)
- [10. Getting Help](#10-getting-help)
- [11. Pydoc](#11-pydoc-)

### Introduction 

`kmer-analysis` is a command-line program designed to analyse biological sequences by counting and characterizing **k-mers**.

The software compares the frequencies and occurrences of oligonucleotides between two input sequence files 
and returns the oligonucleotides that are significantly enriched in one file compared to the other.

### Input format 
The program takes as input a pair of sequence files in fasta format.

### Output format 
The output is a tab-delimted file with one row per oligonucleotide,
and one column per statistics.  
The column content is detailed in the header of the output :

| Columns  | Meaning                      |
|----------|------------------------------|
| seq      | oligomer sequence            |
| id       | oligomer identifier          |
| exp_freq | expected relative frequency  |
| obs_freq | 	observed relative frequency |
| occ      | observed occurrences|
| exp_occ  | expected occurrences|

### 1. Basic usage: Counts Occurrences 

The simplest use of the program is to count how many times each k-mer appears.
```commandline
python3 kmer-analysis.py -i sequence.fa -k 2
```

**Parameters**

| Option | Meaning              |
|--------|----------------------|
| `-i`   | Input Fasta file     |
| `-k`   | Length of the k-mers |

**Example output**
```
AAA     12
AAC     4
AAG     7
```
This means:
- `AAA` appears 12 times 
- `AAC` appears 4 times

### 2. How K-mers are Counted
The program uses a sliding window.  
He scans each sequence using a sliding window of size k and counts all valid k-mers.

**Example**  

For `ATTCG` sequence  with `K=2`.

| Position |k-mer|
|----------|-----|
| 0        |ATT|
| 1        |TTC|
| 2        |TCG|

Total numbers of windows :  
`sequence_lenght - k +1 `

### 3. Computing Frequencies 
Raw counts are useful, however frequencies make it easier to compare sequences of different lengths.

**Output Example**
```
AAA     occ=12      freq=0.031
AAC     occ=2       freq=0.010    
```
**Interpretation :**  

Frequency :  
`occurence / total posisble windows`  

Useful when comparing datasets of different sizes.

### 4. Exploring Result as TSV 
To save results in tabular format :
```commandline
python3 kmer-analysis.py -i sequences.fa -k 2 -o results.tsv
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
P(ATG) =  P(A)*P(T)*P(G)
```
Simple and fast 

Use when : 
- quick analyses 
- random-like sequences
- first exploration

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

Use when : 
- genomic DNA
- promoter analysis
- motif enrichment studies  

A 0 order Markov model uses the observed single-residue frequencies to estimate expected oligonucleotide frequencies, assuming no dependencies between adjacent residues. Under this assumption, sequence positions are treated as independent and identically distributed, corresponding to a Bernoulli model.
 
### 7. Oligonucleotide size

The analysis can be performed using oligonucleotides of length ranging from 1 to 10. Setting the size to 1 corresponds to quantifying the usage of the underlying alphabet within the input sequences. For the identification of putative regulatory sites, we recommend initiating the analysis with hexanucleotides (size = 6) and exploring k-mer lengths between 4 and 8. Biologically relevant patterns that are significantly overrepresented are typically consistently detected across analyses performed with different k-mer sizes.

**Command**
```commandline
python3 kmer-analysis.py -i sequence.fa -k lengh_kmers
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

### 8. Reverse Complement - Strand Options 
**Command** 
```commandline
python3 kmer-analysis.py -i sequences.fa -k 2 -o results.tsv -s both
```
| Option   | Meaning                      |
|----------|------------------------------|
| `single` | Forward                      |
|`both`| Forward + reverse complement |

By selecting the `both` option, the occurrences of each oligonucleotide are aggregated across both strands. This enables the detection of sequence elements that function in an orientation-independent manner. 
### 9. Understanding Output Columns

|                          Columns | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                              |
|---------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                              seq | oligomer sequence                                                                                                                                                                                                                                                                                                                                                                                                                    |
|                               id | oligomer identifier                                                                                                                                                                                                                                                                                                                                                                                                                  |
|  Expected frequencies (exp_freq) | Proportion of a given oligonucleotide that is expected within the sequence set under the selected probabilistic model. It is computed from the expected occurrences by normalizing them with respect to the total number of possible oligonucleotide positions in the dataset. This measure reflects the theoretical probability of observing each oligonucleotide given the assumed model of sequence composition and dependencies. |
|           Frequencies (obs_freq) | 	Relative frequencies are defined as the number of occurrences of each oligonucleotide divided by the total number of occurrences across all oligonucleotides in the dataset.                                                                                                                                                                                                                                                        |
|                Occurrences (occ) | Number of occurrences of each oligonucleotide observed in the dataset. All matches are included, and overlapping occurrences are detected and summed during the count.                                                                                                                                                                                                                                                               |
|   Expected occurrences (exp_occ) | Predicted number of times a given oligonucleotide is expected to appear within the sequence set. This value is computed according to the probabilistic model chosen by the user (see above), which defines how sequence composition and dependencies are taken into account.                                                                                                                                                         |

### 10. Getting Help
List all options :
````commandline
python3 kmer-analisys.py --help
````
### 11. Pydoc 
Automatic Python code documentation : 
- function 
- parameters 

**Command**
````commandline
pydoc kmer-analysis
````

