# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
from keras.models import model_from_json
from Bio import SeqIO

BASE_KEY = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'M': 4}

dataset = 'ABCA4_NCSS'


def check_input(input_seq):
    if not type(input_seq) == str:
        print('The input must be a string.')
        sys.exit(1)
    if not len(input_seq) == 170:
        print('The sequence should be 170-mer.')
    if not input_seq.count('M') == 1:
        print('"M" should be appear only once in the sequence.')


def predict_ESM(input_seq):

    check_input(input_seq)

    input_vec = np.zeros((1, 170, 5))
    for i in range(len(input_seq)):
        try:
            input_vec[0][i][BASE_KEY[input_seq[i]]] = 1
        except KeyError:
            print('Each base must be "A", "C", "G", "T", or "M"')
            sys.exit(1)

    model = model_from_json(open(os.path.join(os.path.dirname(__file__), 'seq_based_DNN.json')).read())
    model.load_weights(os.path.join(os.path.dirname(__file__), 'seq_based_DNN.hdf5'))

    predict = model.predict(input_vec, batch_size=1, verbose=0)[0, 0]

    print('The probability of ESM: {}'.format(predict))

    return predict


def main(argv=sys.argv):
    
    fasta_sequences = SeqIO.parse(open((dataset + '_var.fa.out')),'fasta')

    with open ((dataset + '.csv'), 'w') as file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            sequence = sequence.upper()
            predict = predict_ESM(sequence)
            file.write(name + ',' + str(predict) + '\n')
            

if __name__ == '__main__':
    main()