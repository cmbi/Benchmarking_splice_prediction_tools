# Benchmarking splice prediction tools

This repository contains the scripts used for the analysis of the manuscript **"Benchmarking deep learning splice prediction tools using functional splice assays"**.

All variants and splice prediction scores used for the analysis are included in `data/variant_scores.xlsx`. The datasets are:
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
- Number of splice altering and non splice altering variants at each position in the non-canonical splice site

#### `confusion_matrix.py`
This scripts calculates, for each splice prediction tool used in the analysis, the optimal threshold and corresponding confusion matrix. The format of the confusion matrix in python is:

| TN | FP |
|-|-|
| **FN** | **TP** |

#### `create_vcffile.py`
This script is used to convert the variants into vcf format. This is required for CADD, MMSplice and SpliceAI. The script makes use of the pyhgvs package nad it required a [reference genome file](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/) and [RefSeq transcripts](https://github.com/counsyl/hgvs/blob/master/pyhgvs/data/genes.refGene).

#### `roc.py`
The roc.py script creates the receiver operatur curve (ROC) curve for the dataset. Additionally, it also prints the area under the curve (AUC) for each tool. 

#### `roc_best5tools.py`
This scripts plots the ROC curve for the 5 best tools for the dataset including the AUC. 

## Splice prediction tools

* Alamut 3/4 consensus (consensus of GeneSplicer, MaxEntScan, NNSPLICE and SpliceSiteFinder-like)
* [CADD](https://cadd.gs.washington.edu/score)
* [DSSP](https://github.com/DSSP-github/DSSP)
* [GeneSplicer](https://ccb.jhu.edu/software/genesplicer/)
* [MaxEntScan](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html)
* [MMSplice](https://github.com/gagneurlab/MMSplice_MTSplice)
* [NNSPLICE](https://www.fruitfly.org/seq_tools/splice.html)
* [Spidex](http://tools.genes.toronto.edu/)
* [SpliceAI](https://github.com/Illumina/SpliceAI)
* [SpliceRover](http://bioit2.irc.ugent.be/rover/splicerover)
* SpliceSiteFinder-like
