#imports
import pandas as pd
import numpy as np
import pybedtools
from Bio import SeqIO
import hgvs
from hgvs.easy import *
import sys
sys.path.insert(1, '../')
from functions import reverse_sequence

#define variables
gene = 'ABCA4'
variant = 'NCSS'

length = 410
excel_file = '../data/variant_scores.xlsx'
genome = 'GRCh37'
dataset = gene + '_' + variant

if gene == 'ABCA4':
    chromosome = 1
    # reference fasta can be downloaded from http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.1.fa.gz
    reference_fasta = '../references/Homo_sapiens.GRCh37.dna.chromosome.1.fa'
        
elif gene == 'MYBPC3':
    chromosome = 11
    # reference fasta can be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_000011.9?report=fasta
    reference_fasta = '../references/chr11.fa'

# Read the second column of the excel sheet and store the variants in a list
df = pd.read_excel(excel_file, dataset ,index_col=None, usecols = 'B', engine='openpyxl')
variant_list = df['genomic variant'].tolist()

# Store the variant information so that it can be accessed separately

info = []

hp = hgvs.parser.Parser()
for i in range(len(variant_list)):
    variant = variant_list[i]
    v = hp.parse_hgvs_variant('Chr' + str(chromosome) + genome + ':' + variant)
    
    var_info = []
    var_info.append(v.posedit.pos.start.base)
    var_info.append(v.posedit.pos.end.base)
    var_info.append(v.posedit.edit.ref)
    var_info.append(v.posedit.edit.alt)
    info.append(var_info)

# The BED file defines the sequence range that is written to the fasta file later on
with open ((dataset + '.bed'), 'w') as file:
    for i in info:
        loc = i[0]
        file.write('chr' + str(chromosome) + '\t' + str(loc-(length//2+1)) + '\t' + str(loc+(length//2)) + '\t\t\t' + '-' + '\n')

#Get the sequence for each variant and store it in a fasta file
a = pybedtools.BedTool((dataset + '.bed'))
a = a.sequence(fi = reference_fasta, fo = (dataset + '.fa.out'))

#Change the mutated base in the fasta sequence
fasta_sequences = SeqIO.parse(open((dataset + '.fa.out')),'fasta')

# open the new fasta file to save the mutated sequences
with open ((dataset + '_var.fa.out'), 'w') as file:
    i = 0
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        
        # write the reference sequence to the file
        file.write('>' + variant_list[i] + '\n' + reverse_sequence(sequence) + '\n')

        # variants where one base is changed
        if info[i][0] == info[i][1] and info[i][2] != '':
            # test if the base at the variant position in the sequence is the same as the reference base
            assert sequence[(length//2)] == info[i][2]  
            # change the base at the variant position
            s = list(sequence)
            s[(length//2)] = info[i][3]
            sequence = ''.join(s)
            # test if the base at the variant position in the sequence is now the same as the mutated base
            assert sequence[(length//2)] == info[i][3]  
            file.write('>' + variant_list[i] + '_var\n' + reverse_sequence(sequence) + '\n')
        # filter for variants where one single base is deleted
        elif info[i][0] == info[i][1]:
            s = sequence[:(length//2)] + sequence[(length//2+1):]
            file.write('>' + variant_list[i] + '_var\n' + reverse_sequence(s) + '\n')
        # handle deletions with more bases
        else:
            size = info[i][1] - info[i][0]
            s = sequence[:(length//2)] + sequence[((length//2+1) + size):]
            file.write('>' + variant_list[i] + '_var\n' + reverse_sequence(s) + '\n')
        i += 1

