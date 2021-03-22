#imports 
import pandas as pd
from datetime import date
import hgvs
from hgvs.easy import *

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

genome = 'GRCh37'


#Load the data
df = pd.read_excel(io=path, sheet_name=genename + '_' + variants, engine='openpyxl')
variant = []
for name in df['cDNA variant']:
    variant.append(name)

# Convert the cDNA location to genomic location
# create a dictionary to store the values for the vcf file
values = dict()

# get the values for all the variants
for name in variant:
    
    #create an empty dictionary to store the values
    values[name] = {}
    
    # create a hgvs parser to convert cDNA location to genomic location
    #if the variant is a deletion, we also want to write the previous base to the vcf file
    if 'del' in name:
        hp = hgvs.parser.Parser()
        hgvs_var = hp.parse_hgvs_variant(gene + ':' + name)
        hgvs_var.posedit.pos.start.base = hgvs_var.posedit.pos.start.base -1
    
    else:  
        hp = hgvs.parser.Parser()
        hgvs_var = hp.parse_hgvs_variant(gene + ':' + name)

    hdp = hgvs.dataproviders.uta.connect()

    am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=genome)
    var_t = am.c_to_g(hgvs_var)

    # get the chromosome
    chromosome = int((var_t.ac).split('.')[0][-2:])
    
    # replace 23 with X and 24 with Y 
    if chromosome == '23':
        chromosome = 'X'
    elif chromosome == '24':
        chromosome = 'Y'
        
    values[name]['CHROM'] = chromosome
    
    # get the position
    values[name]['POS'] = var_t.posedit.pos.start.base
    
    # get the reference base
    values[name]['REF'] = var_t.posedit.edit.ref
        
    # get the alternative base
    values[name]['ALT'] = var_t.posedit.edit.alt
    
# Convert to Dataframe
vcf_data = pd.DataFrame.from_dict(values, orient='index')

# Define reference genome
genome = ['reference=GRCh37/hg19', 'contig=<ID=1,length=249250621>', 'contig=<ID=2,length=243199373>',
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
    
    for element in genome:
        file.write('##' + element + '\n')
    
    file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    
    for i in range(len(vcf_data)):
        file.write(str(vcf_data.iloc[i]['CHROM']) + '\t')
        file.write(str(vcf_data.iloc[i]['POS']) + '\t')
        file.write('.' + '\t')
        file.write(str(vcf_data.iloc[i]['REF']) + '\t')
        if vcf_data.iloc[i]['ALT'] == None:
            file.write(str(vcf_data.iloc[i]['REF'][0]) + '\t')
        else:
            file.write(str(vcf_data.iloc[i]['ALT']) + '\t')
        file.write('.' + '\t')
        file.write('.' + '\t')
        file.write('.' + '\n')