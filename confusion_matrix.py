#import
import pandas as pd 
from functions import Find_Optimal_Cutoff, read_scores_from_excel
from sklearn.metrics import confusion_matrix

#define the variants that should be analyzed (ABCA4_NCSS, ABCA4_DI or MYBPC3_NCSS)
variants = 'MYBPC3_NCSS'

#read in the data
# Define the column headers that are used in the dataframe. For DI variants MMSplice, MTSPlice and SPIDEX are excluded.
if 'NCSS' in variants:
    column_names = ['RNA','CADD','DSSP','GeneSplicer', 'MaxEntScan', 'MMSplice', 'NNSPLICE', 'SPIDEX', 'SpliceAI', 'SpliceRover', 'SpliceSiteFinder-like']
else:
     column_names = ['RNA','CADD','DSSP','GeneSplicer', 'MaxEntScan', 'NNSPLICE', 'SpliceAI', 'SpliceRover', 'SpliceSiteFinder-like']
df = read_scores_from_excel('data/variant_scores.xlsx', variants)
df.columns = column_names

#set the threshold for the values that are considered to affect splicing. Everything above the threshold is defined to affect splicing.
threshold = []
for name in column_names[1:]:
    threshold.append(Find_Optimal_Cutoff(df['RNA'], df[name])[0])
    print(name, Find_Optimal_Cutoff(df['RNA'], df[name])[0])

# create a new dataframe to store the classification and add the classification of the different tools to the dataframe
classification = pd.DataFrame(df['RNA']) 
i = 0
for name in column_names[1:]:
    classification[name] = (df[name] > threshold[i]).astype('int')
    i += 1
#add the classification of the Alamut 3/4 consensus
classification['consensus'] = ((classification['SpliceSiteFinder-like'] + classification['MaxEntScan'] + classification['NNSPLICE'] + classification['GeneSplicer']) > 2)

#write the confusion matrix (format [[TN FP][FN TP]]) to a csv file
cm = []
i = 0
for name in column_names[1:]:
    c = confusion_matrix(classification.RNA.values, classification[name].values)
    print(name)
    print(c)
    cm.append([name, c[1,1], c[0,1],c[0,0],c[1,0], threshold[i]])
    i += 1
print('consensus')
print(confusion_matrix(classification.RNA.values, classification['consensus'].values))
cm_consensus = confusion_matrix(classification.RNA.values, classification['consensus'].values)
cm.append(['Alamut 3/4 consensus', cm_consensus[1,1], cm_consensus[0,1], cm_consensus[0,0],cm_consensus[1,0]])

statistics = pd.DataFrame.from_records(cm)
statistics.to_csv(variants + '_cm.csv', index=False, header=['Splice prediction tool', 'TP', 'FP', 'TN', 'FN', 'threshold'])

print(df['DSSP'])