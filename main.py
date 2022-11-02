
"""(x,y,z) = max(x,y,z) - min(x,y,z)"""
import random

import numpy as np
import math
import matplotlib.pyplot as plt
from Bio.PDB import *
from qdpy import *
from Bio.PDB.PDBParser import PDBParser
"import map_elites.cvt as cvt_map_elites"
"archive = cvt_map_elites.compute(2, 5, box, n_niches=10000, max_evals=1e6, log_file=open('cvt_arm.dat', 'w'), params=px)"

import os

import subprocess

def fitness_function():
    fitness = 0
    "START OMEGAFOLD"

    os.chdir("C:\\Omegafold\\")
    subprocess.run("omegafold fasta.txt C:\\Omegafold --num_cycle 1", shell=True)


    "END OMEGAFOLD"

    # Load PDB file
    parser = PDBParser()
    structure = parser.get_structure('1st chain' ,'C:\\Omegafold\\1st chain.pdb')

    for model in structure:
        for chain in model:
            a = 0
            b = 0
            for residue in chain:

                for atom in residue:

                    b += 1
                    a += atom.bfactor
            fitness = ( a / b)
            print(fitness)


    "generate protein sequences"

    # Load PDB file
    parser = PDBParser()
    structure = parser.get_structure('1st chain', 'C:\\Omegafold\\1st chain.pdb')
    res=[]
    for model in structure:
        for chain in model:
            for residue in chain:
                res.append(residue.get_resname())

    "convert residues to one letter code"
    res1 =[]
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

    "convert residues into fasta format"
    fasta = open("C:\\Omegafold\\fasta.txt", "w")
    fasta.write(">1st chain")
    fasta.write("\n")


    """mutation of residues"""

    for i in res1:
        """mutation of residues 90% of the time we put a normal residue but 10% of the time we put a random residue"""
        if random.random() < 0.9:
            fasta.write(i)
        else:
            fasta.write(random.choice(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))

    fasta.close()
    return fitness

fitness=0
while fitness < 35:

    fitness = fitness_function()