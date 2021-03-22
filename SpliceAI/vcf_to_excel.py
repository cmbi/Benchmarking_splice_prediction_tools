#import
import pandas as pd

#define variables
gene = 'ABCA4'
variants = "NCSS"

vcf_file = gene + '_' + variants + '_output.vcf'
excel_file = gene + '_' + variants + '_spliceai_scores.xlsx'

# define a list to store the scores
scores = []

# open the file and read the scores
with open(vcf_file, 'r') as f:
    for line in f.readlines():
        if 'SpliceAI=' in line:
            variant = []
            # first store the genomic location of the variant
            variant.append(int(line.split('\t')[1]))
            result = line.split('\t')[-1]
            result = result.split('|')
            # store the scores
            for i in range(2,10):
                variant.append(float(result[i]))
            scores.append(variant)

# convert the scores list to a dataframe and write to an excel sheet.
df = pd.DataFrame.from_records(scores)
# add column names
df.columns = ['genomic_location', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
df.to_excel(excel_file)