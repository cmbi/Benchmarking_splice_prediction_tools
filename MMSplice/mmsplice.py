#import
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table, predict_save
from mmsplice.utils import max_varEff
import cyvcf2
from gtfparse import read_gtf
import pandas as pd

# Define variables
variants = 'NCSS'
gene = 'MYBPC3'

# Get necessary files
# 1) Standard human gene annotation file in GTF format -> filtering for protein coding genes is recommended
# can be downloaded from http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gtf = 'Homo_sapiens.GRCh37.75.gtf'

# 2) VCF file
vcf = ('../data/' + gene + '_' + variants + '_variants.vcf')

 3) Human reference fasta file from ensemble/gencode, chromosome names have to match with GTF file
if gene == 'ABCA4':
    # reference fasta can be downloaded from http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.1.fa.gz
    fasta = 'Homo_sapiens.GRCh37.dna_rm.chromosome.1.fa'
if gene == 'MYBPC3':
    # reference fasta can be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_000011.9?report=fasta
    fasta = 'chr11.fa'

csv = ('mtsplice_predictions_' + gene + '_' + variants + '.csv')

# Filter the gtf file for variants in the gene of interest
# read all transcripts from the gtf file
df = pd.read_csv(gtf, sep='\t', header=None)

# filter for transcripts of the gene of interest
g = df[8].str.contains(gene)
result = df[g]
print(result.head())

# write the result to a csv file for reuse
new_gtf = ('Homo_sapiens_' + gene + '_all.GRCh37.75.gtf')
result.to_csv(new_gtf , sep='\t', index=False, header=None)