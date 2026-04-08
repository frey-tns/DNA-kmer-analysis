# DNA-kmer-analysis

## Description 
Logiciel d'analyse des k-mers dans les séquences d'ADN

Le script: 
- input: fichier de séquences fasta :
  - Lecture du fichier fasta
  - Stat sur les séquences
    - Nombre de séquences
    - Somme des longueurs
- output : on ajoute progressivement (en testant à chaque étape) différentes statistiques :
  - occurrences observées
  - occurrences attendues
  - fréquence observée (fréquence / position)
  - Lecture d’un modèle probabiliste (modèle de Bernoulli ou Markov, sous forme d’une matrice de transition)
  - fréquence attendue (probabilité a priori, prior) en fonction de différents modèles probabilistes, Bernoulli ou Markov
  - ratio obs / exp
  - P-valeur
  - E-valeur

## Installation 
1. Prérequis
  - Python (version >= 3.x) 
2. Installation des dépendances
... A COMPLETER ...

## Input 
Fichier FASTA contenant les séquences ADN.

→ Modifier chemin d'accès dans le code :
path_outputs = "/chemin/vers/fichier.fasta"

## Utilisation
Lancer le script :
python scipt.py (NOM A MODIFIER)

## Output
Les résultats sont exportés dans un fichier HTML:
results_test.html (NOM A MODIFIER)

Contient :
  - occurrences observées (occ)
  - occurrences attendues (exp_occ)
  - fréquence observée (fréquence / position)
  - Lecture d’un modèle probabiliste (modèle de Bernoulli ou Markov, sous forme d’une matrice de transition)
  - fréquence attendue (probabilité a priori, prior) en fonction de différents modèles probabilistes, Bernoulli ou Markov (exp_freq)
  - ratio obs / exp
  - P-valeur (occ_P)
  - E-valeur (ooc_E)

## Remarque
- Le chemin des fichiers est acutellement codé en dur (Windows) 
