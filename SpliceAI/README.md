## [SpliceAI](https://github.com/Illumina/SpliceAI)

1. Install SpliceAI (more details can be found on the SpliceAI github page)
2. Download the reference genome file [`hg19.fa`](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/)
3. Create the variant vcf file with `create_vcffile.py`
4. Run SpliceAI in the command line (example for ABCA4 NCSS variants): 

   `spliceai -I data/ABCA4_NCSS_variants.vcf -O SpliceAI/ABCA4_NCSS_output.vcf -R references/hg19.fa -A grch37 -D 500`
   
5. Run `vcf_to_excel.py` to write the SpliceAI scores to an excel file

