#imports 
import pandas as pd
from datetime import date
import pyhgvs as hgvs
from pyfaidx import Fasta
from pyhgvs.utils import read_transcripts

# Define variables
genename = 'ABCA4'
variants = 'NCSS'

path = 'data/variant_scores.xlsx'
if genename == 'ABCA4':
    gene = 'NM_000350.2'
    chromosome = 1
elif genename == 'MYBPC3':
    gene = 'NM_000256.3'
    chromosome = 11

#Load the data
df = pd.read_excel(io=path, sheet_name=genename + '_' + variants, engine='openpyxl')
data = []
for name in df['cDNA variant']:
    data.append(name)

#read the genome 
# the reference genome can be downloaded from: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/
genome = Fasta('references/hg19.fa')

# Read RefSeq transcripts into a python dict.
# The RefSeq transcripts can be downloaded from: https://github.com/counsyl/hgvs/blob/master/pyhgvs/data/genes.refGene
with open('references/genes.refGene') as infile:
    transcripts = read_transcripts(infile)

# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)

# Store the variant information in a list
vcf = []
for v in data:
    chrom, offset, ref, alt = hgvs.parse_hgvs_name(gene + ':' + v, genome, get_transcript=get_transcript)
    vcf.append([chrom, offset, ref, alt])

# Define reference genome
ref_genome = ['reference=GRCh37/hg19', 'contig=<ID=1,length=249250621>', 'contig=<ID=2,length=243199373>',
            'contig=<ID=3,length=198022430>', 'contig=<ID=4,length=191154276>', 'contig=<ID=5,length=180915260>',
            'contig=<ID=6,length=171115067>', 'contig=<ID=7,length=159138663>', 'contig=<ID=8,length=146364022>',
            'contig=<ID=9,length=141213431>', 'contig=<ID=10,length=135534747>','contig=<ID=11,length=135006516>',
            'contig=<ID=12,length=133851895>', 'contig=<ID=13,length=115169878.', 'contig=<ID=14,length=107349540>',
            'contig=<ID=15,length=102531392>', 'contig=<ID=16,length=90354753>', 'contig=<ID=17,length=81195210>',
            'contig=<ID=18,length=78077248>', 'contig=<ID=19,length=59128983>', 'contig=<ID=20,length=63025520>',
            'contig=<ID=21,length=48129895>', 'contig=<ID=22,length=51304566>', 'contig=<ID=X,length=155270560>',
            'contig=<ID=Y,length=59373566>']

# Write to file
with open(('data/' + genename + '_' + variants + '_variants.vcf'),'w') as file:
    file.write('##fileformat=VCFv4.2\n')

    today = date.today().strftime("%m%d%Y")
    
    file.write('##fileDate=' + str(today) + '\n')
    
    for element in ref_genome:
        file.write('##' + element + '\n')
    
    file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    
    for i in range(len(vcf)):
        file.write(str(chromosome) + '\t')
        file.write(str(vcf[i][1]) + '\t')
        file.write('.' + '\t')
        file.write(str(vcf[i][2]) + '\t')
        file.write(str(vcf[i][3]) + '\t')
        file.write('.' + '\t')
        file.write('.' + '\t')
        file.write('.' + '\n')
