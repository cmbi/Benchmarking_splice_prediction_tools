import pandas as pd
from numpy import nan
import numpy as np

def delta_score(df, name, index):
    '''calulate the delta score 
    @param df: dataframe that contains the scores
    @param name: name of the splice prediction tool
    @param index: index of the scores
    '''
    # Define the name of the programs for which a delta score is calculated with corresponding max value
    delta = dict()
    delta['SSFL'] = 100
    delta['MES'] = 12
    delta['NNS'] = 1
    delta['GS'] = 15
    delta['SpliceRover'] = 1
    delta['DSSP'] = 1
    
    wt = name + '_wt'
    var = name + '_var'
    if df.at[index,wt] == 0:
        score = float(df.at[index,var])/delta[name]
    else: 
        score = (float(df.at[index,var])-float(df.at[index,wt]))/float(df.at[index,wt])
    return np.absolute(score)

def read_scores_from_excel(file, sheetname, fillna = True, diall = False):
    ''' This function takes an excel sheet with splice prediction scores, fills missing values with 0, 
    and calculates delta scores if necessary. It stores the resulting primary scores in a dataframe.
    @param file: name of the excel file
    @param sheetname: name of the sheet with the splice prediction scores
    @param fillna: If set to True, missing values are replaced with 0.
    @param diall: If set to True, it includes all tools for DI variants, even the ones that cannot predict scores
    Returns primary (delta) scores for all tools for each variant in the dataset
    '''
    
    # store the scores in a dataframe
    di = pd.read_excel(file, sheetname)

    # replace missing values with 0
    if fillna == True:
        di = di.replace(nan, 0)

    # create a dictionary to store the values
    delta_scores = dict()

    for index in di.index:

        element = []

        # get the % mutant RNA
        value = di.at[index,'% Mutant RNA']
        if value > 20:
            element.append(1)
        else:
            element.append(0)

        # add the absolute value of the CADD score 
        element.append(np.absolute(di.at[index,'CADD']))

        # add the DSSP score 
        element.append(delta_score(di, 'DSSP', index))

        # add the GeneSplicer score 
        element.append(delta_score(di, 'GS', index))

        # add the MaxEntScan score 
        element.append(delta_score(di, 'MES', index))

        if 'NCSS' in sheetname or diall == True:

            # add the absolute value of the MMsplice score 
            element.append(np.absolute(di.at[index,'MMSplice']))

            # add the absolute value of the MMsplice score
            element.append(np.absolute(di.at[index,'MTSplice']))

        # add the NNSPLICE score 
        element.append(delta_score(di, 'NNS', index))

        if 'NCSS' in sheetname or diall == True:

            # add the S-SCAP score
            element.append(np.absolute(di.at[index,'SCAP']))

            # add the absolute value of the SPIDEX score
            element.append(np.absolute(di.at[index,'Spidex']))

        # add the SpliceAI score 
        element.append(di.at[index,'SpliceAI'])

        # add  the SpliceRover score 
        element.append(delta_score(di, 'SpliceRover', index))

        # add the SSFL score 
        element.append(delta_score(di, 'SSFL', index))


        delta_scores[index] = element

    delta_df = pd.DataFrame(delta_scores)
    delta_df = delta_df.transpose()
    return delta_df

def reverse_sequence(s):
    ''' Converts a sequence into the sequence of the complementary strand'''
    new_sequence = ''
    for base in s:
        if base == 'A':
            new_sequence = new_sequence + 'T'
        elif base == 'T':
            new_sequence = new_sequence + 'A'
        elif base == 'G':
            new_sequence = new_sequence + 'C'
        elif base == 'C':
            new_sequence = new_sequence + 'G'
        else:
            new_sequence = new_sequence + base
    return new_sequence[::-1]