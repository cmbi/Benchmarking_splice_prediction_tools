# Benchmarking splice prediction tools

This repository contains the scripts used for the analysis of the manuscript **"Benchmarking deep learning splice prediction tools using functional splice assays"**.

All variants and splice prediction scores used for the analysis are included in `variant_scores.xlsx`. The datasets are:
- ABCA4 NCSS (non-canonical splice site) variants
- ABCA4 DI (deep intronic) variants
- MYBPC3 NCSS variants


## Analysis scripts

#### `analysis_variants.py`
This script provides information about the variants in the dataset:
- number of variants:  
    - number of non splice altering variants
    - number of splice altering variants
- number of variants located close to the splice donor site (SDS)
    - number of splice altering variants at the SDS
    - number of non splice altering variant at the SDS
- number of variants located close to the acceptor donor site (SAS)
    - number of splice altering variants at the SAS
    - number of non splice altering variant at the SAS
- Number of splice altering an non splice altering variants at each position in the non-canonical splice site

#### `confusion_matrix.py`
This scripts calculates, for each splice prediction tool used in the analysis, the optimal threshold and corresponding confusion matrix. The format of the confusion matrix in python is:

| TN | FP |
|-|-|
| **FN** | **TP** |

#### `roc.py`
The roc.py script creates the receiver operatur curve (ROC) curve for the dataset. Additionally, it also prints the area under the curve (AUC) for each tool. 


## Splice prediction tools

#### Alamut 3/4, GeneSplicer, NNSPLICE, MaxEntScan and SpliceSiteFinder-like

#### CADD

#### DSSP
DSSP required the following input sequences:
- An input for acceptor site should be a 140-mer string with the AG at positions 69 and 70
- An input for donor site should be a 140-mer string with the GT at positions 71 and 72

To reproduce the analysis, the following steps are required:
1. The input sequence in the right format for DSSP can be obtained with `get_input_sequence.ipynb` or `get_input_sequence_di.ipynb` for deep intronic variants.
2. `AS_DSSP.py` and `DS_DSSP.py` are used to calculated the scores for acceptor and donor sites, respectively.

#### MMSplice and MTSplice

#### Spidex

#### SpliceAI

#### SpliceRover



