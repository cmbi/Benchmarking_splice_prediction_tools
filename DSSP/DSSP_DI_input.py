#import
import pandas as pd
import numpy as np
import pybedtools
from Bio import SeqIO
import hgvs
from hgvs.easy import parser
import sys
sys.path.insert(1, '../')
from functions import reverse_sequence

gene = 'ABCA4'
length = 140

# define variables
excel_file = '../variant_scores.xlsx'
genome = 'GRCh37'
chromosome = 1
dataset = gene + '_DI'
# reference fasta can be downloaded from http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.1.fa.gz
reference_fasta = '../references/Homo_sapiens.GRCh37.dna.chromosome.1.fa'

# Read in the data 
# Read the second column of the excel sheet and store the variants in a list
df = pd.read_excel(excel_file, dataset ,index_col=None, engine='openpyxl')

# store the variant information
info = []
hp = hgvs.parser.Parser()
for i in range(len(df['genomic variant'])):
    variant = df['genomic variant'][i]
    v = hp.parse_hgvs_variant('Chr' + str(chromosome) + 'GRCh37:' + variant)
    
    # Store the variant information so that it can be accessed separately
    var_info = []
    var_info.append(v.posedit.pos.start.base)
    var_info.append(v.posedit.pos.end.base)
    var_info.append(v.posedit.edit.ref)
    var_info.append(v.posedit.edit.alt)
    
    # add if it affects a donor or acceptor
    var_info.append(df['affects'][i])

    # add the distance to the splice site
    var_info.append(df['position ss'][i])

    info.append(var_info)

# 2) Create the BED file
# The BED file defines the sequence range that is written to the fasta file later on
with open ((dataset + '.bed'), 'w') as file:
    for i in range(len(info)):
        # An input for acceptor site should be a 140-mer string with the AG at positions 69 and 70
        if df['affects'][i] == 'acceptor':
            var_loc = info[i][0]
            loc = var_loc - info[i][5] + 1
        
        # An input for donor site should be a 140-mer string with the GT at positions 71 and 72
        else:
            var_loc = info[i][0]
            loc = var_loc + info[i][5]
            
        # If the variant contains a deletion, a longer sequence is needed to end up with a 140nt long sequence
        if 'del' in df['cDNA variant'][i]:
            # get the length of the deletion
            l = info[i][1] - info[i][0]
            # For deletions affecting the donor site additional bases are added to the end of the sequence
            if df['affects'][i] == 'donor':
                file.write('chr' + str(chromosome) + '\t' + str(loc-(length//2+2)) + '\t' + str(loc+(length//2-1)+l) + '\t\t\t' + '-' + '\n')
            # For deletions affecting the acceptor site additonal bases are added to the beginning of the sequence
            else:
                file.write('chr' + str(chromosome) + '\t' + str(loc-(length//2+1)-l) + '\t' + str(loc+(length//2)) + '\t\t\t' + '-' + '\n')
                
        else:     
            file.write('chr' + str(chromosome) + '\t' + str(loc-(length//2+1)) + '\t' + str(loc+(length//2-1)) + '\t\t\t' + '-' + '\n')


#Get the sequence for each variant and store it in a fasta file
a = pybedtools.BedTool((dataset + '.bed'))
a = a.sequence(fi = reference_fasta, fo = (dataset + '.fa.out'))

fasta_sequences = SeqIO.parse(open((dataset + '.fa.out')),'fasta')
# open the new fasta file to save the mutated sequences
with open ((dataset + '_donor.fa.out'), 'w') as file:
    with open ((dataset + '_acceptor.fa.out'), 'w') as file2:
        i = 0
        for fasta in fasta_sequences:
            # get the name and the sequence
            name, sequence = fasta.id, str(fasta.seq)
            wt_sequence = reverse_sequence(sequence)[140:]

            # get the location of the variant
            if df['affects'][i] == 'acceptor':
                loc = length//2 + info[i][5] - 1
            else:
                loc = length//2 - info[i][5]

            # variants where one base is changed
            if info[i][0] == info[i][1] and info[i][2] != '':
                assert sequence[loc] == info[i][2]
                # change the base at the variant position
                l = list(sequence)
                l[loc] = info[i][3]
                s = ''.join(l)
                # test if the base at the variant position in the sequence is now the same as the mutated base
                assert s[loc] == info[i][3]  

            # filter for variants where one single base is deleted
            elif info[i][0] == info[i][1]:
                s = sequence[:loc+1] + sequence[(loc+2):]

            # handle deletions with more bases
            else:
                size = info[i][1] - info[i][0] + 1
                s = sequence[:loc + 1] + sequence[(loc + size + 1):]

            s = reverse_sequence(s)

            # Check if the bases at the donor/acceptor position are correct
            if df['affects'][i] == 'acceptor':
                assert s[68:70] == 'AG'
            else:
                assert s[70:72] in ['GT','GC']

            # write the result to a file
            if reverse == True:
                s = reverse_sequence(s)
                
            if df['affects'][i] == 'acceptor':
                file2.write('>' + df['cDNA variant'][i] + '\n' + wt_sequence + '\n') 
                file2.write('>' + df['cDNA variant'][i] + '_var\n' + s + '\n') 
            else:
                file.write('>' + df['cDNA variant'][i] + '\n' + wt_sequence + '\n') 
                file.write('>' + df['cDNA variant'][i] + '_var\n' + s + '\n') 

            i += 1