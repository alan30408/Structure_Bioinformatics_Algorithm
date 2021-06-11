#!/usr/bin/python
# coding=utf-8
import os
import sys
import numpy as np
from argparse import ArgumentParser




__author__ = "Chun-Wei Tung"

def parseParams():

    arg_parser = ArgumentParser()
    arg_parser.add_argument('-i', dest='rna_input', type=str, required=True, default='testfasta.sec')
    arg_parser.add_argument('--output', dest='output_file', type=str, default='result.bpseq')
    arg_parser.add_argument('--min-loop-length', dest='length', type=int, default=4)
    arg_parser.add_argument('--score-GC', dest='gc', type=int, default=1)
    arg_parser.add_argument('--score-AU', dest='au', type=int, default=1)
    arg_parser.add_argument('--score-GU', dest='gu', type=int, default=0)
    params = arg_parser.parse_args()
    return params


def delta(l, m):

    params = parseParams()
    if l=='A' and m=='U':
        return params.au
    elif l=='U' and m=='A':
        return params.au
    elif l=='G' and m=='C':
        return params.gc
    elif l=='C' and m=='G':
        return params.gc
    elif l=='G' and m=='U':
        return params.gu
    elif l=='U' and m=='G':
        return params.gu
    else:
        return 0


def matrix_builder(seq):

    params = parseParams()
    length = len(seq)
    mat = [i[:] for i in [[0] * length] * length]
    for n in range(1, length):
        for j in range(n, length):
            i = j - n
            case1 = mat[i+1][j-1] + delta(seq[i], seq[j])
            case2 = mat[i+1][j]
            case3 = mat[i][j-1]
            if i < j - params.length:
                record = []
                for k in range(i+1, j-1):
                    record.append(mat[i][k]+mat[k+1][j])
                case4 = max(record)
                mat[i][j] = max(case1, case2, case3, case4)
            else:
                mat[i][j] = max(case1, case2, case3)
    return mat


def traceback_to_alignment(mat, i, j, seq, match):
    if i <= j:
        if mat[i][j] == mat[i+1][j]:
            traceback_to_alignment(mat, i+1, j, seq, match)
        elif mat[i][j] == mat[i][j-1]:
            traceback_to_alignment(mat, i, j-1, seq, match)
        elif mat[i][j] == mat[i+1][j-1] + delta(seq[i], seq[j]):
            match.append([i+1, j+1, str(seq[i]), str(seq[j])])
            traceback_to_alignment(mat, i+1, j-1, seq, match)
        else:
            for k in range(i+1, j-1):
                if mat[i][j] == mat[i][k] + mat[k+1][j]:
                    traceback_to_alignment(mat, i, k, seq, match)
                    traceback_to_alignment(mat, k+1, j, seq, match)
                    break
    return match


def bracket_dot(seq, match):
    dot = ["." for i in range(len(seq))]
    for x in range(len(match)):
        dot[match[x][0]-1] = "("
        dot[match[x][1]-1] = ")"
    return dot


def seq_file_reader(input_file):

    seq = ''
    with open(input_file, "r") as f:
        header = f.readline().strip()
        for line in f:
            seq += line.strip()
    return header, seq


def main():

    params = parseParams()

    header, seq = seq_file_reader(params.rna_input)
    matrix = matrix_builder(seq)
    match = traceback_to_alignment(matrix, 0, len(seq)-1, seq, [])
    dot = bracket_dot(seq, match)
    with open(params.output_file, "w") as file:
        file.write("Filename:" + params.rna_input + "\n")
        file.write("Min-loop: " + str(params.length) + "\n")
        file.write("GC: " + str(params.gc) + "\n")
        file.write("AU: " + str(params.au) + "\n")
        file.write("GU: " + str(params.gu) + "\n")
        file.write("Score: " + str(matrix[0][len(seq)-1]) + "\n")
        for i in range(len(seq)):
            number = i+1

            position = 0
            for x in range(len(match)):
                if number == match[x][0]:
                    position = match[x][1]
                elif number == match[x][1]:
                    position = match[x][0]

            file.write(str(number) + " " + str(seq[i]) + " " + str(position) + "\n")
        for j in range(len(dot)):
            file.write(dot[j])

if __name__ == "__main__":
    main()
