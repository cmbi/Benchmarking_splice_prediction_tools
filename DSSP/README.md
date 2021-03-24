## [DSSP](https://github.com/DSSP-github/DSSP)

#### Download reference genome sequences
Reference fasta sequences for chromosome 1 (ABCA4) and chromosome 11(MYBPC3) can be downloaded from the following links into the `references` folder.
- [ABCA4](http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.1.fa.gz)
- [MYBPC3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000011.9?report=fasta)

DSSP requires the following input sequences:
- Input for the SAS should be a 140-mer string with the canonical splice site AG at position 69 and 70
- Input for the SDS should be a 140-mer string with the canonical splice site GT at position 71 and 72

To reproduce the analysis, the following steps are required:
1. Generation of the input sequences for DSSP: Separate scripts for NCSS variants (`DSSP_NCSS_input.py`) and DI variants (`DSSP_DI_input.py`) are available. Both scripts produce two separate output files, one for variants located at the SDS and one for variants located at the SAS. For ABCA4 NCSS variants for instance the files are called `ABCA4_NCSS_acceptor.fa.out` and `ABCA4_NCSS_donor.fa.out`. The output files contain two sequences for each variant, one is the wild type sequence and the second one contains the variant.
2. Calulation of the DSSP scores: `AS_DSSP.py` and `DS_DSSP.py` are used to calculated the scores for acceptor and donor sites, respectively. DSSP can be called from the command line with the following command:

`$ python DS_DSSP.py -I input.fasta -O output.txt`

The DSSP scripts are modified so that they can handle the fasta files with the wild type and variant sequence. The result is written to `ABCA4_NCSS_acceptor.csv` for ABCA4 NCSS SAS variants. The left column contains the score for the wild type sequence and the right column the score for the variant sequence. The final score is calculated as the difference between the score of the two sequences.
