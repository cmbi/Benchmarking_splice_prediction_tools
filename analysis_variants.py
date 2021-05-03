# import
import numpy as np
import pandas as pd
from numpy import nan
from collections import Counter
from functions import read_scores_from_excel

# define the variants that should be analyzed (ABCA4_NCSS, ABCA4_DI or MYBPC3_NCSS)
variants = 'ABCA4_NCSS'

# calculate the missing scores for each tool
df = pd.read_excel('data/variant_scores.xlsx', variants, engine='openpyxl')
df.isnull().sum()

# replace missing values with 0
df = df.replace(nan, 0)

# print statistics about the variants
print('# variants: ', df.shape[0])
print('# non splice altering variants: ',
      df[df['% Mutant RNA'] <= 20].count()['cDNA variant'])
print('# splice altering variants: ',
      df[df['% Mutant RNA'] > 20].count()['cDNA variant'])
print('donor: ', df[df['affects'] == 'donor'].count()['cDNA variant'])
print('acceptor: ', df[df['affects'] == 'acceptor'].count()['cDNA variant'])

# print the number of variants that affect the SAS and SDS
donor = df['affects'] == "donor"
acceptor = df['affects'] == "acceptor"
sa = df['% Mutant RNA'] > 20
nsa = df['% Mutant RNA'] <= 20
print('SDS, splice altering: ', df[donor & sa].shape[0])
print('SDS, normal splicing: ', df[donor & nsa].shape[0])
print('SAS, splice altering: ', df[acceptor & sa].shape[0])
print('SAS, normal splicing: ', df[acceptor & nsa].shape[0])

# print the position of the variants in the NCSS
locations = []
# only calculate the distribution for NCSS variants
if 'NCSS' in variants:
    for index in df.index:
        # get the ss positions
        pos = df.at[index, 'position ss']
        # check if the variant alters splicing
        if df.at[index, '% Mutant RNA'] > 20:
            sa = 'splice altering'
        else:
            sa = 'normal splicing'
        # check if the variant affects the donor or acceptor
        affects = df.at[index, 'affects']
        locations.append((str(pos) + ' ' + sa + ' ' + affects))
# count how often a variant is located at a certain position in the NCSS motif and print the result
print(Counter(locations))
