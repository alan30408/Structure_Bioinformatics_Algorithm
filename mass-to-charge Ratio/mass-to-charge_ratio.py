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
    arg_parser.add_argument('-i', dest='input', type=str, default='HFEEDMGRK')
    # arg_parser.add_argument('-o', dest='output', type=str, default='dssp_matrix.tsv')
    params = arg_parser.parse_args()
    return params

def open_file(input):

    seq = ''
    with open(input, "r") as f:
        header = f.readline().strip()
        for line in f:
            seq += line.strip()
    return seq

def cut_peptide(seq):

    cut_peptide = []
    for ind, res in enumerate(seq):
        cut_peptide.append([seq[:ind+1], seq[ind:]])

    return cut_peptide

def monoisotopic_masses():

    H = 1.007825
    C = 12
    N = 14.00307
    O = 15.99491
    S = 31.97207

    return H, C, N, O, S

def calculate_mass_aa(function):

    H , C, N, O, S = function
    mass = {
    'A' : C*3+H*5+O+N,      #+2*H+O, #C3H5ON
    'R' : C*6+H*12+N*4+O,   #+2*H+O, #C6H12ON4
    'N' : C*4+H*6+N*2+O*2,  #+2*H+O, #C4H6O2N2
    'D' : C*4+H*5+N+O*3,    #+2*H+O, #C4H5O3N
    'C' : C*3+H*5+N+O+S,    #+2*H+O, #C3H5ONS
    'Q' : C*5+H*8+N*2+O*2,  #+2*H+O, #C5H8O2N2
    'E' : C*5+H*7+N+O*3,    #+2*H+O, #C5H7O3N
    'G' : C*2+H*3+N+O,      #+2*H+O, #C2H3ON
    'H' : C*6+H*7+N*3+O,    #+2*H+O, #C6H7ON3
    'I' : C*6+H*11+N+O,     #+2*H+O, #C6H110N
    'L' : C*6+H*11+N+O,     #+2*H+O, #C6H110N
    'K' : C*6+H*12+N*2+O,   #+2*H+O, #C6H12ON2
    'M' : C*5+H*9+N+O+S,    #+2*H+O, #C5H9ONS
    'F' : C*9+H*9+N+O,      #+2*H+O, #C9H9ON
    'P' : C*5+H*7+N+O,      #+2*H+O, #C5H7ON
    'S' : C*3+H*5+N+O*2,    #+2*H+O, #C3H5O2N
    'T' : C*4+H*7+N+O*2,    #+2*H+O, #C4H7O2N
    'W' : C*11+H*10+N*2+O,  #+2*H+O, #C11H10ON2
    'Y' : C*9+H*9+N+O*2,    #+2*H+O, #C9H9O2N
    'V' : C*5+H*9+N+O,      #+2*H+O #C5H9ON
    }

    return mass

def calculate_ion_peptide(peptide):

    b_ion, y_ion = [], []
    mono_table = calculate_mass_aa(monoisotopic_masses())
    for seg in peptide:
        b_mass, y_mass = 0, 0
        for res in seg[0]:
            b_mass += mono_table[res]
        b_mass = b_mass + 1.007825
        b_ion.append(b_mass)
        for res in seg[1]:
            y_mass += mono_table[res]
        y_mass  = y_mass + 3*1.007825 + 15.99491
        y_ion.append(y_mass)

    return b_ion, y_ion

def main():

    params = parseParams()
    warnings.filterwarnings('ignore')

    sequence = params.input
    peptide = cut_peptide(sequence)
    b_ion, y_ion = calculate_ion_peptide(peptide)
    #
    print("Seq\t#\tb_ion\ty_ion\t#(+1)")
    for i in range(len(b_ion)):
        print("%s\t%s\t%.2f\t%.2f\t%s" %(sequence[i], i+1, round(b_ion[i], 2), round(y_ion[i], 2), 9-i))


if __name__ == "__main__":
    main()
