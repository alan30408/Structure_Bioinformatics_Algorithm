#!/usr/bin/python
# coding=utf-8
import os
import sys
import numpy as np
from argparse import ArgumentParser
from Bio.PDB import *
import warnings


__author__ = "Chun-Wei Tung"

def parseParams():

    arg_parser = ArgumentParser()
    arg_parser.add_argument('-i', dest='input', type=str, default='homology_model.pdb, 6nh4.1.pdb')
    params = arg_parser.parse_args()
    return params

def open_file(input):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('X', input)

    return structure

def root_mean_square_deviation(homology, reference):

    parser = PDBParser(PERMISSIVE=1)
    chain_list = []
    record_homo = []
    record_ref = []
    homo_model = homology[0]
    ref_model = reference[0]
    for chain in homo_model:
        if chain.id != "_":
            chain_list.append(chain.id)
            for residue in chain:
                if residue.get_id()[0] == " ":
                    record_homo.append(residue)

    for chain in ref_model:
        if chain.id in chain_list:
            for residue in chain:
                if residue.get_id()[0] == " ":
                    record_ref.append(residue)

    length = len(record_homo)
    for i in range(length):
        index = i
        try:
            while record_homo[index].get_resname() != record_ref[i].get_resname() and index <= len(record_homo):
                record_homo.pop(index)
        except:
            continue

    for i in range(length):
        index = i
        try:
            while record_homo[i].get_resname() != record_ref[index].get_resname() and index <= len(record_homo):
                record_ref.pop(index)
        except:
            continue

    vector_homo = {}
    vector_ref = {}
    for i in range(len(record_homo)):
        vector_homo[i] = record_homo[i]["CA"].get_vector()
    for i in range(len(record_ref)):
        vector_ref[i] = record_ref[i]["CA"].get_vector()

    distance = 0
    for i in range(len(vector_homo)):
        x = np.square(vector_homo[i][0] - vector_ref[i][0])
        y = np.square(vector_homo[i][1] - vector_ref[i][1])
        z = np.square(vector_homo[i][2] - vector_ref[i][2])
        distance += (x + y + z)

    RMSD = np.sqrt(distance/len(vector_homo))
    print("RMSD :", RMSD)

def main():

    params = parseParams()
    warnings.filterwarnings('ignore')

    homology = params.input.split()[0]
    reference = params.input.split()[1]
    homology_struc = open_file(homology)
    reference_struc = open_file(reference)
    root_mean_square_deviation(homology_struc, reference_struc)














if __name__ == "__main__":
    main()
