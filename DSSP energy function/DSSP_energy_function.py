#!/usr/bin/python
# coding=utf-8
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

def parseParams():

    arg_parser = ArgumentParser()
    arg_parser.add_argument('-i', dest='input', type=str, default='5jxv.pdb')
    arg_parser.add_argument('-o', dest='output', type=str, default='dssp_matrix.tsv')
    params = arg_parser.parse_args()
    return params

def open_file(input):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('X', input)

    return structure

def DSSP_vector(structure):

    vector = {}
    model = structure[0]
    for chain in model:
        for residue_i in chain:
            vector[residue_i.get_id()[1]] = [residue_i["C"].get_vector(), residue_i["O"].get_vector(),
                                               residue_i["N"].get_vector(), residue_i["H"].get_vector()]
    return vector

def DSSP_energy(vector):

    energy_record = {}
    #distance
    for i in vector:
        for j in vector:
            # distance: C-N
            cn = np.linalg.norm(vector[i][0] - vector[j][2])
            # distance: C-H
            ch = np.linalg.norm(vector[i][0] - vector[j][3])
            # distance: O-N
            on = np.linalg.norm(vector[i][1] - vector[j][2])
            # distance: O-H
            oh = np.linalg.norm(vector[i][1] - vector[j][3])

            energy = 0.084*(1/on+1/ch-1/oh-1/cn)*332*4.2 #kJ/mol
            if i != j:
                energy_record[(i, j)] = energy


    return energy_record

def DSSP_matrix(energy, length):

    matrix = [[0 for i in range(length)] for i in range(length)]
    for i in range(length):
        for j in range(length):
            if i != j:
                matrix[i][j] = energy[(i+1, j+1)]

    return matrix


def DSSP_heatmap(matrix, length):

    xlabel = [i for i in range(1, length+1, 4)]
    ylabel = [i for i in range(1, length+1, 4)]

    fig = plt.figure()
    plt.yticks(ylabel)
    plt.xticks(xlabel)
    plt.imshow(matrix, cmap=plt.cm.jet)
    plt.colorbar()
    plt.title("Heatmap")
    plt.savefig("Heatmap.jpg")
    plt.close()

def main():

    params = parseParams()
    warnings.filterwarnings('ignore')

    structure = open_file(params.input)
    vector = DSSP_vector(structure)
    energy = DSSP_energy(vector)
    matrix = DSSP_matrix(energy, len(vector))
    DSSP_heatmap(matrix, len(vector))

    with open(params.output, 'w') as file:
        for i in range(len(vector)):
            for j in range(len(vector)):
                file.write(str(matrix[i][j])+"\t")
            file.write("\n")













if __name__ == "__main__":
    main()
