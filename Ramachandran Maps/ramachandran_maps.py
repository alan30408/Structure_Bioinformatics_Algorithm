#!/usr/bin/python
# coding=utf-8
import os
import sys
import numpy as np
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from Bio.PDB import *
import math

__author__ = "Chun-Wei Tung"

def parseParams():

    arg_parser = ArgumentParser()
    arg_parser.add_argument('-i', dest='input', type=str, default=['5ire.pdb', '1igt.pdb', '4ire.pdb'])
    params = arg_parser.parse_args()
    return params

def open_file(input):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('X', input)

    for model in structure:
        record = []
        for chain in model:
            residue_list = chain.get_list()

            index = 0


            # reduce Het item
            for residue in residue_list:
                residue_id = residue.get_id()
                hetfield = residue_id[0]
                if hetfield [ 0 ] == "H" :
                    index += 1
            for i in range(index):
                residue_list.pop(-1)

            # find four point of item
            for residue in chain:
                index = 0

                for atom in residue:
                    if atom.get_name() == 'N':
                        n = atom.get_vector()
                        index += 1
                    elif atom.get_name() == 'CA':
                        ca = atom.get_vector()
                        index += 1
                    elif atom.get_name() == 'C':
                        c = atom.get_vector()
                        index += 1
                if index != 3:
                    break

                # find previous C
                if residue.id[1] - 2 < 0:
                    continue
                elif residue.id[1] > len(residue_list):
                    continue
                else:
                    for atom in residue_list[residue.id[1] - 2]:
                        if atom.get_name() == 'C':
                            cp = atom.get_vector()
                phi = torsion(cp, n, ca, c)

                # find next N
                if residue.id[1] >= len(residue_list):
                    continue
                else:
                    for atom in residue_list[residue.id[1]]:
                        if atom.get_name() == 'N':
                            nn = atom.get_vector()
                psi = torsion(n, ca, c, nn)

                record.append([phi, psi])


    return record

def torsion(v1, v2, v3, v4):

    v1 = np.array(list(v1))
    v2 = np.array(list(v2))
    v3 = np.array(list(v3))
    v4 = np.array(list(v4))
    q1 = np.subtract(v2, v1)
    q2 = np.subtract(v3, v2)
    q3 = np.subtract(v4, v3)

    q12 = np.cross(q1, q2)
    q23 = np.cross(q2, q3)

    n1 = q12 / np.sqrt(np.dot(q12, q12))
    n2 = q23 / np.sqrt(np.dot(q23, q23))

    u1 = n2
    u3 = q2 / (np.sqrt(np.dot(q2, q2)))
    u2 = np.cross(u3, u1)
    cos = np.dot(n1, u1)
    sin = np.dot(n1, u2)

    theta = -math.atan2(sin, cos)
    theta_deg = np.degrees(theta)

    return theta_deg

def plot_figure(record, name):

    name = name.split(".")
    record = np.array(record)
    fig = plt.figure
    plt.title('Ramachandran map - %s.png'%name[0])
    plt.scatter(record[:, 0], record[:, 1], c='b', s = 1)
    plt.xlabel(r'$\Phi$ Angle [degree]')
    plt.xlim(-180, 180)
    plt.xticks([-180, -120, -60, 0, 60, 120, 180])
    plt.ylabel(r'$\Psi$ Angle [degree]')
    plt.ylim(-180, 180)
    plt.yticks([-180, -120, -60, 0, 60, 120, 180])
    plt.grid()
    plt.savefig('Ramachandran map - %s.png'%name[0])
    plt.close()
    return fig


def main():

    params = parseParams()

    for i in range(len(params.input)):
        input = params.input[i]
        torsion_angel = open_file(input)
        figure = plot_figure(torsion_angel, input)



if __name__ == "__main__":
    main()
