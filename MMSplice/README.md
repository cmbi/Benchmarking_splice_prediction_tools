## [MMSplice and MTSplice](https://github.com/gagneurlab/MMSplice_MTSplice)

1. Install MMSplice and all necessary packages (information on MMSplice github page)
2. Download the necessary reference genome sequences and genome annotation files"
   *  [Human gene annotation](http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)
   *  [Reference fasta sequence chromosome 1](http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.1.fa.gz)
   *  [Reference fasta sequence chromosome 11](https://www.ncbi.nlm.nih.gov/nuccore/NC_000011.9?report=fasta)
3. Create the variant vcf file with `create_vcffile.py'
4. Run `mmsplice_scores.py` on the dataset. The parameter tissue_specificity determines if MMSplice or MTSplice is used.
5. The scores are stored in a csv file called `mmsplice_ABCA4_NCSS.csv` for ABCA4 NCSS variants and tissue_specificity = False.
