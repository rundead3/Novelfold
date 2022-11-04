
"""(x,y,z) = max(x,y,z) - min(x,y,z)"""
import random
import qdpy as qd
import numpy as np
import math
import matplotlib.pyplot as plt
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser


import os

import subprocess

def fitness_function():
    """COMPUTE BLOBS"""
    os.chdir("C:\\Users\\42077\\Omegaforl\\blobulator\\blobulator-main\\")
    subprocess.run("python compute_blobs.py --fasta C:\\Users\\42077\\Omegafold\\randseq.txt --cutoff 0.4 --minBlob 4 --oname .\\Batch\\Outputs\\", shell=True, stdout=subprocess.DEVNULL)

    "CALCULATE NICHES"
    nicheCoords = []
    for i in range(0, population):
        nicheCoords.append(init_niches(i))

    print(nicheCoords)


    "START OMEGAFOLD"
    os.chdir("C:\\Users\\42077\\OmegaFold")
    subprocess.run("omegafold C:\\Users\\42077\\Omegafold\\randseq.txt C:\\Users\\42077\\Omegaforl\\res1 --num_cycle 1", shell=True)
    "END OMEGAFOLD"

    # Load PDB file
    parser = PDBParser()


    "generate protein sequences"
    bestRes = []
    bestFit = 0
    "get structure data"
    for i in range(0, population):
        structure = parser.get_structure(str(i)+"th chain", "C:\\Users\\42077\\Omegaforl\\res1\\"+str(i)+"th chain.pdb")
        for model in structure:
            for chain in model:
                a = 0
                b = 0
                res = []
                for residue in chain:
                    res.append(residue.get_resname())
                    for atom in residue:
                        b += 1
                        a += atom.bfactor
                fitness = ( a / b)
                if (fitness>bestFit):
                    bestFit = fitness
                    bestRes = res
                print(fitness)

    "convert residues into fasta format"
    res1 = triToFasta(bestRes)

    fasta = open("C:\\Users\\42077\\Omegafold\\randseq.txt", "w")

    for i in range(0, population):
        fasta.write(">"+str(i)+"th chain")
        fasta.write("\n")

        """mutation of residues"""

        for i in res1:
            """mutation of residues 95% of the time we put a normal residue but 5% of the time we put a random residue"""
            if random.random() < 0.95:
                fasta.write(i)
            else:
                fasta.write(random.choice(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))
        fasta.write("\n")

    fasta.close()
    return fitness

def triToFasta(res):
    "convert residues to one-letter code"
    res1 = []
    for i in res:
        if i == 'ALA':
            res1.append('A')
        if i == 'ARG':
            res1.append('R')
        if i == 'ASN':
            res1.append('N')
        if i == 'ASP':
            res1.append('D')
        if i == 'CYS':
            res1.append('C')
        if i == 'GLN':
            res1.append('Q')
        if i == 'GLU':
            res1.append('E')
        if i == 'GLY':
            res1.append('G')
        if i == 'HIS':
            res1.append('H')
        if i == 'ILE':
            res1.append('I')
        if i == 'LEU':
            res1.append('L')
        if i == 'LYS':
            res1.append('K')
        if i == 'MET':
            res1.append('M')
        if i == 'PHE':
            res1.append('F')
        if i == 'PRO':
            res1.append('P')
        if i == 'SER':
            res1.append('S')
        if i == 'THR':
            res1.append('T')
        if i == 'TRP':
            res1.append('W')
        if i == 'TYR':
            res1.append('Y')
        if i == 'VAL':
            res1.append('V')
    return res1

def init_niches(i):
    csv_file = open("C:\\Users\\42077\\Omegaforl\\blobulator\\blobulator-main\\Batch\\Outputs\\"+str(i)+"th.csv", "r")
    csv_list = csv_file.readlines()
    csv_file.close()

    csv_list.pop(0)
    nichx = float(0)
    nichy = float(0)
    for i in csv_list:
        nichx+=float(i.split(',')[5])
        nichy+=float(i.split(',')[6])
    return nichx/len(csv_list), nichy/len(csv_list)



population = 20
fitness=0
genN = 0
while fitness < 100:
    fitness = fitness_function()
    genN = genN+1
    print(genN)

