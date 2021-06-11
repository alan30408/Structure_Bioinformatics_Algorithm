#!/usr/bin/python
# coding=utf-8
import os
import sys
import numpy as np
from argparse import ArgumentParser
import warnings


__author__ = "Chun-Wei Tung"

def parseParams():

    arg_parser = ArgumentParser()
    arg_parser.add_argument('-i', dest='input', type=str, default='b_y_spectrum.txt')
    arg_parser.add_argument('-o', dest='output', type=str, default='de_novo_peptide.txt')
    params = arg_parser.parse_args()
    return params

def calculate_mass_aa():

    H = 1.007825
    C = 12
    N = 14.00307
    O = 15.99491
    S = 31.97207
    mass = {
    'A' : C*3+H*5+O+N,#+2*H+O, #C3H5ON
    'R' : C*6+H*12+N*4+O,#+2*H+O, #C6H12ON4
    'N' : C*4+H*6+N*2+O*2,#+2*H+O, #C4H6O2N2
    'D' : C*4+H*5+N+O*3,#+2*H+O, #C4H5O3N
    'C' : C*3+H*5+N+O+S,#+2*H+O, #C3H5ONS
    'Q' : C*5+H*8+N*2+O*2,#+2*H+O, #C5H8O2N2
    'E' : C*5+H*7+N+O*3,#+2*H+O, #C5H7O3N
    'G' : C*2+H*3+N+O,#+2*H+O, #C2H3ON
    'H' : C*6+H*7+N*3+O,#+2*H+O, #C6H7ON3
    'I' : C*6+H*11+N+O,#+2*H+O, #C6H110N
    'L' : C*6+H*11+N+O,#+2*H+O, #C6H110N
    'K' : C*6+H*12+N*2+O,#+2*H+O, #C6H12ON2
    'M' : C*5+H*9+N+O+S,#+2*H+O, #C5H9ONS
    'F' : C*9+H*9+N+O,#+2*H+O, #C9H9ON
    'P' : C*5+H*7+N+O,#+2*H+O, #C5H7ON
    'S' : C*3+H*5+N+O*2,#+2*H+O, #C3H5O2N
    'T' : C*4+H*7+N+O*2,#+2*H+O, #C4H7O2N
    'W' : C*11+H*10+N*2+O,#+2*H+O, #C11H10ON2
    'Y' : C*9+H*9+N+O*2,#+2*H+O, #C9H9O2N
    'V' : C*5+H*9+N+O,#+2*H+O #C5H9ON
    }

    return mass


def read_file(spetrum):

    fragment = []
    with open(spetrum, 'r') as file:
        for line in file:
            fragment.append(float(line.strip()))

    return fragment

def de_novo_sequence(fragment):

    size = len(fragment)
    distance = np.array([[0]*size]*size).astype(np.float)
    mono_table = calculate_mass_aa()
    for i in range(size):
        for j in range(i, size):
            distance[i][j] = abs(fragment[i] - fragment[j])

    count = 0
    possibility = {}
    for i in range(size):
        for j in range(i, size):
            if i != j:
                for aa in mono_table:
                    if distance[i][j] <= mono_table[aa] + 0.08 and distance[i][j] >= mono_table[aa] - 0.08:
                        if (i, j) not in possibility:
                            possibility[(i, j)] = [aa]
                        else:
                            possibility[(i, j)].append(aa)
    b_ion = 1.007825
    for i in range(2):
        for aa in mono_table:
            if fragment[i] - b_ion <= mono_table[aa] + 0.08 and fragment[i] - b_ion >= mono_table[aa] - 0.08:
                start_point_b = i
                start_aa_b = aa
    y_ion = 3*1.007825 + 15.99491
    for aa in mono_table:
        if fragment[1-start_point_b] - y_ion <= mono_table[aa] + 0.08 and fragment[1-start_point_b] - y_ion >= mono_table[aa] - 0.08:
            start_point_y = 1-start_point_b
            start_aa_y = aa

    order_y = []
    list_y = [start_aa_y]
    y_ions = []
    for i, aa in enumerate(possibility):
        if aa[0] == start_point_y:
            ind = aa[1]
            if len(possibility[aa]) != 1:
                list_y.append(possibility[aa][0])
            else:
                list_y.append(possibility[aa][0])
            order_y.append((start_point_y, aa[1]))
            for j, next_aa in enumerate(possibility):
                if next_aa[0] == ind and next_aa != (12, 14):
                    ind = next_aa[1]
                    if len(possibility[next_aa]) != 1:
                        list_y.append(possibility[next_aa][0])
                    else:
                        list_y.append(possibility[next_aa][0])
                    y_ions.append(next_aa[1])
                    order_y.append(next_aa)

    list_b = [start_aa_b]
    index = []
    order_b = []
    for i, aa in enumerate(possibility):
        if aa[0] == start_point_b:
            ind = aa[1]
            if len(possibility[aa]) != 1:
                index.append((aa, 1))
                list_b.append(possibility[aa][0])
            else:
                list_b.append(possibility[aa][0])
            order_b.append((start_point_b, aa[1]))
            for j, next_aa in enumerate(possibility):
                if next_aa[0] == ind and next_aa[1] not in y_ions:
                    ind = next_aa[1]
                    if len(possibility[next_aa]) != 1:

                        list_b.append(possibility[next_aa][0])
                        index.append((next_aa, len(list_b)-1))

                    else:
                        list_b.append(possibility[next_aa][0])
                    order_b.append(next_aa)
    # print("b ion\ty ions")
    # for i, j in zip(order_b, order_y):
    #     print(i,"\t", j)
    return list_b, list_y, index, possibility

def list_to_seq(list):

    seq = ''
    for i in list:
        seq = seq+i
    return seq

def main():

    params = parseParams()
    warnings.filterwarnings('ignore')

    spectrum = params.input
    fragment = read_file(spectrum)
    list_b, list_y, index, possibility = de_novo_sequence(fragment)

    seq_b = list_to_seq(list_b)
    seq_y = list_to_seq(list_y)

    change_list = []
    for i in index:
        if possibility[i[0]] not in change_list:
            change_list.append(possibility[i[0]])

    with open(params.output, 'w') as file:
        file.write('The result of de novo peptide for %s\n\n' %spectrum)
        file.write('seq of b_ion : %s\n'%seq_b)
        file.write('seq of y_ion : %s\n\n'%seq_y)
        file.write('de novo peptide : %s\n\n'%seq_b)
        file.write('According to the calculation result, here exist more than one peptide:\n')
        for i in change_list:
            file.write('\tAmino acid %s can change to %s\n'%(i[0], i[1]))
        file.write('Thus, we get %s different peptides\n'%np.power(2, len(index)))

if __name__ == "__main__":
    main()
