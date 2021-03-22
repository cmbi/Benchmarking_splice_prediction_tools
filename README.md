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
GeneSplicer, NNSPLICE, MaxEntScan and SpliceSiteFinder-like were accessed through Alamut Visual Software version 2.13. The Alamut 3/4 consensus approach is a consensus of the four mentioned Alamut tools.

#### [CADD](https://cadd.gs.washington.edu/score)

#### [DSSP](https://github.com/DSSP-github/DSSP)

#### [MMSplice and MTSplice](https://github.com/gagneurlab/MMSplice_MTSplice)

#### [Spidex](http://tools.genes.toronto.edu/)

#### [SpliceAI](https://github.com/Illumina/SpliceAI)

#### [SpliceRover](http://bioit2.irc.ugent.be/rover/splicerover)



