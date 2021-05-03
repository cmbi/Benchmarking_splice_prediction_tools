# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from sklearn.metrics import roc_curve
from sklearn import metrics
import matplotlib.font_manager as font_manager
from functions import read_scores_from_excel

# define the variants that should be analyzed (ABCA4_NCSS, ABCA4_DI or MYBPC3_NCSS)
variants = 'MYBPC3_NCSS'

# Define the column headers that are used in the dataframe. For DI variants MMSplice, MTSPlice and SPIDEX are excluded.
if 'NCSS' in variants:
    column_names = ['RNA','CADD','DSSP','GeneSplicer', 'MaxEntScan', 'MMSplice', 'NNSPLICE', 'SPIDEX', 'SpliceAI', 'SpliceRover', 'SpliceSiteFinder-like']
else:
     column_names = ['RNA','CADD','DSSP','GeneSplicer', 'MaxEntScan', 'NNSPLICE', 'SpliceAI', 'SpliceRover', 'SpliceSiteFinder-like']

# Import the scores, calculate delta scores and store them in a dataframe
delta_df = read_scores_from_excel('data/variant_scores.xlsx', variants)
delta_df.columns = column_names

names = column_names[1:]

# prepare the data 
# 1) List with classification (0,1)
label = []
for index in delta_df.index:
    value = delta_df.at[index,'RNA']
    if value > 0.2:
        label.append(1)
    else:
        label.append(0)
label = np.array(label)

# 2) list with probabilities predicted by the splicing prediction program 
probabilities = []
for name in names:
    probabilities.append(np.array(delta_df[name].tolist()))

# 3) Add the alamut consensus
if 'NCSS' in variants:
    loc = [2,3,5,9]
else:
    loc = [2,3,4,7]
alamut3 = []
for i in range(len(probabilities[0])):
    p = [probabilities[j][i] for j in loc]
    largest_integer = max(p) 
    p.remove(largest_integer)
    second_largest_integer = max(p)
    p.remove(second_largest_integer)
    third_largest_integer = max(p)
    alamut3.append((largest_integer + second_largest_integer + third_largest_integer)/3)

probabilities.insert(0,alamut3)
names.insert(0,'Alamut 3/4')

# 4) Define the colors for the lines    
colors = {'Alamut 3/4': '#20E2E7', 'CADD' : '#0B3954', 'DSSP' : '#63B0CD',
          'GeneSplicer' : '#0353A4', 'MaxEntScan' : '#12664F', 'MMSplice' : '#95BF74',
          'NNSPLICE' : '#D87CAC', 'SPIDEX' : '#DB4C40', 'SpliceAI' : '#EEC643',
          'SpliceRover' : '#F4743B', 'SpliceSiteFinder-like' : '#8B1E3F'}

# create a dictionary to store the AUC values
aucs = {}
ftrates = {}

# Get the 5 highest auc values
for i in range(len(probabilities)):
    prob = probabilities[i]
    fper, tper, thresholds = roc_curve(label, prob, pos_label=1) 
    auc = metrics.roc_auc_score(label, prob)
    aucs[names[i]] = auc
    ftrates[names[i]] = [fper, tper]
    
sorted_aucs = sorted((value, key) for (key,value) in aucs.items())
l = len(sorted_aucs)
tools = [i[1] for i in sorted_aucs[l-5:l]]
print(tools)

# Plot the ROC curve
plt.figure(figsize=(10,10))
for t in tools:
    plt.plot(ftrates[t][0], ftrates[t][1], color=colors[t], label=t + ': ' + "{0:0.2f}".format(aucs[t]), linewidth=1)

font_prop = font_manager.FontProperties(size=18)
    
plt.plot([0, 1], [0, 1], color='#CED4DA', linestyle='--')
plt.xlabel('False Positive Rate', size=18)
plt.ylabel('True Positive Rate', size = 18)
plt.title(('ROC Curve ' + variants + ' variants'), size = 20)
plt.tick_params(labelsize=18)
plt.legend(prop=font_prop)
plt.savefig(('figures/ROC_5_' + variants + '.svg'),format='svg', dpi=1200)
plt.show()