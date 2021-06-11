import os
import sys
import numpy as np
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from Bio.PDB import *
import glob
from fpdf import FPDF
import warnings

__author__ = "Chun-Wei Tung"

# helix p_alpha relative probabilities
p_a = {
    'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12,  # builders
    'K': 1.07, 'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77,  # indifferent
    'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53  # breakers
}

# helix builder/indifferent/breaker class weights
w_a = {aa: (1 if score > 1.1 else -1 if score < 0.75 else 0.5) for aa, score in p_a.items()}

# strand p_beta relative probabilities
p_b = {
    'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23,  'L': 1.22, 'T': 1.20, 'W': 1.19,  # builders
    'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80,  # indifferent
    'K': 0.74, 'S': 0.72, 'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26  # breakers
}

# strand builder/indifferent/breaker class weights
# changed to 2/0/-1 encoding for convenience, so that a window of 5 is a strand core if sum(window) >= 5.
w_b = {aa: (2 if score > 1 else -1 if score < 0.78 else 0) for aa, score in p_b.items()}

# sequence of 5JXV
# seq    = 'MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'

# reference secondary structure for accuracy calculation
ss_ref = '-SSSSSSS----SSSSSSS---HHHHHHHHHHHHHH-----SSSSS----SSSSS-'

def parseParams():

    arg_parser = ArgumentParser()
    arg_parser.add_argument('-i', dest='input', type=str, default='5jxv.pdb')
    arg_parser.add_argument('-o', dest='output', type=str, default='Secondary_Structure_Prediction_CF.txt')
    params = arg_parser.parse_args()
    return params

def open_file(input):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('X', input)
    ppb = PPBuilder()
    model = structure[0]
    for pp in ppb.build_peptides(model):
        seq = pp.get_sequence()

    return seq


# start to find helix
def seq_weight_helix(seq):
    seq_w_a = []
    for i in range(len(seq)):
        seq_w_a.append(w_a[seq[i]])

    return seq_w_a

def seq_weight_sheet(seq):
    seq_w_b = []
    for i in range(len(seq)):
        seq_w_b.append(w_b[seq[i]])

    return seq_w_b

def helix_finding(seq_w_a):
    helix = []
    for i in range(len(seq_w_a)):
        j = i+5
        sum = 0
        if j < len(seq_w_a):
            for x in range(i, j+1):
                sum += seq_w_a[x]
            if sum >= 4:
                start, end = helix_extending(int(i), int(j))
                if [start, end] not in helix:
                    helix.append([start, end])
    return helix

def helix_extending(start, end):

    while (float(sum([p_a[x] for x in seq[end - 3:end + 1]])) / float(4)) > 1 and end < len(seq) - 1:
        end += 1
    while (float(sum([p_a[x] for x in seq[start-1:start + 3]])) / float(4)) > 1 and start > 0:
        start -= 1

    return start, end

def sheet_finding(seq_w_b):
    sheet = []
    for i in range(len(seq_w_b)):
        j = i+4
        h, b = 0, 0
        if j < len(seq_w_b):
            for x in range(i, j+1):
                if seq_w_b[x] == 2: h+=1
                elif seq_w_b[x] == -1: b+=1

            if h >= 3 and b<= 1:
                start, end = sheet_extending(i, j)
                if [start, end] not in sheet:
                    sheet.append([start, end])
    return sheet

def sheet_extending(start, end):

    while (float(sum([p_b[x] for x in seq[end - 3:end + 1]])) / float(4)) > 1 and end < len(seq) - 1:
        end += 1
    while (float(sum([p_b[x] for x in seq[start-1:start + 3]])) / float(4)) > 1 and start > 0:
        start -= 1

    return start, end

def region_overlap(region_a, region_b):

    return (region_a[0] <= region_b[0] <= region_a[1]) or (region_b[0] <= region_a[0] <= region_b[1])

def conflict(helix, sheet):
    for index_i, i in enumerate(helix):
        for index_j, j in enumerate(sheet):
            if region_overlap(i, j):
                segment = [max(i[0], j[0]), min(i[1], j[1])]
                pa = sum([p_a[x] for x in seq[segment[0]:segment[1]+1]])
                pb = sum([p_b[x] for x in seq[segment[0]:segment[1]+1]])

                if pa >= pb:
                    if j[0] < segment[0] and j[1] <= segment[1]:
                        sheet[index_j] = ([j[0], segment[0]-1])
                    elif j[0] >= segment[0] and j[1] > segment[1]:
                        sheet[index_j] = ([segment[1]+1, j[1]])
                    elif j[0] < segment[0] and j[1] > segment[1]:
                        sheet[index_j] = ([j[0], segment[0]-1])
                        sheet.append([segment[1]+1, j[1]])
                    else:
                        sheet[index_j] = ([0, 0])
                else:
                    if i[0] < segment[0] and i[1] <= segment[1]:
                        helix[index_i] = ([i[0], segment[0]-1])
                    elif i[0] >= segment[0] and i[1] > segment[1]:
                        helix[index_i] = ([segment[1]+1, i[1]])
                    elif i[0] < segment[0] and i[1] > segment[1]:
                        helix[index_i] = ([i[0], segment[0]-1])
                        helix[index_i].append([segment[1], i[1]])
                    else:
                        helix[index_i] = ([0, 0])
    return helix, sheet

params = parseParams()


seq = open_file(params.input)
seq_w_a = seq_weight_helix(seq)
helix = helix_finding(seq_w_a)
seq_w_b = seq_weight_sheet(seq)
sheet = sheet_finding(seq_w_b)
helix, sheet = conflict(helix, sheet)

result = ['-' for i in range(len(seq))]
for i in helix:
    if i != [0, 0]:
        for x in range(i[0], i[1]+1):
            result[x] = 'H'

for i in sheet:
    if i != [0, 0]:
        for x in range(i[0], i[1]+1):
            result[x] = 'S'

with open(params.output, 'w') as file:

    file.write(str(seq) + "\n")
    for i in result:
        file.write(i)







