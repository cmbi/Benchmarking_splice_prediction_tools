# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
from keras.models import model_from_json
from Bio import SeqIO

BASE_KEY = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

dataset = 'ABCA4_DI'


def check_input(input_seq):

    if not type(input_seq) == str:
        print('The input must be a string.')
        sys.exit(1)

    if not len(input_seq) == 140:
        print('The sequence should be 140-mer.')

    if not input_seq[70:72] == 'GT':
        print('The position of 71 and 72 should be "AG".')


def DS_DSSP(input_seq):
    
    check_input(input_seq)
    
    input_vec = np.zeros((1, 140, 5))
    for i in range(len(input_seq)):
        try:
            input_vec[0][i][BASE_KEY[input_seq[i]]] = 1
        except KeyError:
            print('Each base must be "A", "C", "G", "T", or "N"')
            sys.exit(1)
            
    model = model_from_json(open(os.path.join(os.path.dirname(__file__), 'DS_model.json')).read())
    model.load_weights(os.path.join(os.path.dirname(__file__), 'DS_model.hdf5'))

    predict = model.predict(input_vec, batch_size=1, verbose=0)[0, 0]
    
    print ('Donor site probability: {}'.format(predict))

    return predict


def main():
    fasta_sequences = SeqIO.parse(open((dataset + '_donor.fa.out')),'fasta')
    
    i = 0

    with open ((dataset + '_donor.csv'), 'w') as file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            sequence = sequence.upper()
            predict = DS_DSSP(sequence)
            
            if (i % 2) == 0:
                file.write(name + ',' + str(predict))
            else:
                file.write(',' + str(predict) + '\n')
            
            i += 1


if __name__ == '__main__':
    main()