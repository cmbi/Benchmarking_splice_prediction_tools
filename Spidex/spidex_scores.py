#import
import pandas as pd
from cyvcf2 import VCF

#define variables
gene = 'ABCA4'
variants = 'NCSS'

if gene == 'ABCA4':
    chromosome = 1
elif gene == 'MYBPC3':
    chromosome = 11

#Load the variants
# create a dataframe to store the variants
data = pd.DataFrame(columns=['location', 'ref', 'alt', 'score'])

locations = []

# load the variants and store them in de dataframe
for variant in VCF(('../data/' + gene + '_' + variants + '_variants.vcf')): 
    data = data.append({'location': variant.POS, 'ref': variant.REF, 'alt': variant.ALT[0]}, ignore_index=True)
    locations.append(variant.POS)
    
#Get the scores for chromsome of interest from the spidex.txt file
scores = []

# open the file
# file can be downloaded from http://tools.genes.toronto.edu/
with open('hg19_spidex.txt') as f:
    header = f.readline()
    for line in f:
        content = line.split('\t')
        # check if the variant is located on chromosome 11
        if content[0] == str(chromosome):
            # check if the variant position is included into the ABCA4 variant list
            if int(content[1]) in locations:
                scores.append(content)
        elif content[0] == str(chromosome+1):
            break

# create a dictionary to store for each position the score for each alternative base
dict_scores = {}

for row in scores:
    position = int(row[1])
    alt_base = row[4]
    score = row[5]
    if position in dict_scores.keys():
        dict_scores[position].update({alt_base : score})
    else:
        dict_scores[position] = {alt_base : score}

# assign a score to each variant 
for index, row in data.iterrows():
    try:
        score = dict_scores[row['location']][row['alt']]
    except:
        score = 'n.a.'
    data.at[index,'score'] = score

data.to_csv(('spidex_' + gene + '_' + variants + '.csv'))