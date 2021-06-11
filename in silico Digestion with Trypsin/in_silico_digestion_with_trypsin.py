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
    arg_parser.add_argument('-i', dest='input', type=str, default='P07327.fasta')
    params = arg_parser.parse_args()
    return params

def open_file(input):

    seq = ''
    with open(input, "r") as f:
        header = f.readline().strip()
        for line in f:
            seq += line.strip()
    return seq

def cut_trypsin(seq):

    cut_peptide = []
    start = 0
    for ind, res in enumerate(seq):
        if ind+1 < len(seq):
            if res == 'K' and seq[ind+1] != 'P':
                cut_peptide.append(seq[start:ind+1])
                start = ind+1
            elif res == 'R' and seq[ind+1] != 'P':
                cut_peptide.append(seq[start:ind+1])
                start = ind+1
        elif ind+1 == len(seq):
            cut_peptide.append(seq[start:ind+1])

    return cut_peptide

def monoisotopic_masses():

    H = 1.007825
    C = 12
    N = 14.00307
    O = 15.99491
    S = 31.97207

    return H, C, N, O, S

def average_masses():
    H = 1.008
    C = 12.01
    N = 14.01
    O = 16.00
    S = 32.06

    return H, C, N, O, S

def calculate_mass_aa(function):

    H , C, N, O, S = function
    mass = {
    'A' : C*3+H*5+O+N,      #C3H5ON
    'R' : C*6+H*12+N*4+O,   #C6H12ON4
    'N' : C*4+H*6+N*2+O*2,  #C4H6O2N2
    'D' : C*4+H*5+N+O*3,    #C4H5O3N
    'C' : C*3+H*5+N+O+S,    #C3H5ONS
    'Q' : C*5+H*8+N*2+O*2,  #C5H8O2N2
    'E' : C*5+H*7+N+O*3,    #C5H7O3N
    'G' : C*2+H*3+N+O,      #C2H3ON
    'H' : C*6+H*7+N*3+O,    #C6H7ON3
    'I' : C*6+H*11+N+O,     #C6H110N
    'L' : C*6+H*11+N+O,     #C6H110N
    'K' : C*6+H*12+N*2+O,   #C6H12ON2
    'M' : C*5+H*9+N+O+S,    #C5H9ONS
    'F' : C*9+H*9+N+O,      #C9H9ON
    'P' : C*5+H*7+N+O,      #C5H7ON
    'S' : C*3+H*5+N+O*2,    #C3H5O2N
    'T' : C*4+H*7+N+O*2,    #C4H7O2N
    'W' : C*11+H*10+N*2+O,  #C11H10ON2
    'Y' : C*9+H*9+N+O*2,    #C9H9O2N
    'V' : C*5+H*9+N+O       #C5H9ON
    }

    return mass

def calculate_mass_peptide(peptide):

    mono = []
    mono_table = calculate_mass_aa(monoisotopic_masses())
    for seg in peptide:
        mass = 0
        for res  in seg:
            mass += mono_table[res]
        mass = mass + 2 * 1.007825 + 15.99491
        mono.append(mass)
    average = []
    ave_table = calculate_mass_aa(average_masses())
    for seg in peptide:
        mass = 0
        for res in seg:
            mass += ave_table[res]
        mass = mass + 2 * 1.008 + 16.00
        average.append(mass)

    return mono, average

def main():

    params = parseParams()
    warnings.filterwarnings('ignore')

    sequence = open_file(params.input)
    peptide = cut_trypsin(sequence)
    mono, average = calculate_mass_peptide(peptide)

    record = []
    for i in range(len(peptide)):
        print('peptide :\t\t\t', peptide[i])
        print('monoisotopic masses:\t\t',mono[i])
        print('average masses:\t\t\t',average[i])
        print("-------------------------------------------------------------------------------")

if __name__ == "__main__":
    main()
