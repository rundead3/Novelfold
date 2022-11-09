"""(x,y,z) = max(x,y,z) - min(x,y,z)"""
import random
import numpy as np
import math
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser
import os
import subprocess
from chainLogic import chainLogic


def init_population(popSize, length):
    fasta = open("C:\\Users\\42077\\Omegafold\\randseq.txt", "w")

    for i in range(0, popSize):
        fasta.write(">" + str(i) + "th chain")
        fasta.write("\n")

        """mutation of residues"""
        fasta.write(str(random_chain(length)))

        fasta.write("\n")


def random_chain(length):
    strSeq = ""
    for i in range(0, length):
        strSeq += str(random.choice(
            ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))
    return strSeq


def blobulate():
    """COMPUTE BLOBS"""
    os.chdir("C:\\Users\\42077\\Omegaforl\\blobulator\\blobulator-main\\")
    subprocess.run(
        "python compute_blobs.py --fasta C:\\Users\\42077\\Omegafold\\randseq.txt --cutoff 0.4 --minBlob 4 --oname .\\Batch\\Outputs\\",
        shell=True, stdout=subprocess.DEVNULL)


def fold():
    "START OMEGAFOLD"
    os.chdir("C:\\Users\\42077\\OmegaFold")
    subprocess.run("omegafold C:\\Users\\42077\\Omegafold\\randseq.txt C:\\Users\\42077\\Omegaforl\\res1 --num_cycle 1",
                   shell=True)
    "END OMEGAFOLD"


def fitness_function():

    # Create individuals
    chainList = []
    bestFit = 0
    bestIndex = 0
    for i in range(0, population):
        chainList.append(chainLogic(i))
        print(chainList[i])
        if chainList[i].fitScore > bestFit:
            bestIndex = i
            bestFit = chainList[i].fitScore




    "convert residues into fasta format"
    res1 = tri_to_fasta(chainList[bestIndex])

    fasta = open("C:\\Users\\42077\\Omegafold\\randseq.txt", "w")

    for i in range(0, population):
        fasta.write(">" + str(i) + "th chain")
        fasta.write("\n")

        """mutation of residues"""
        fasta.write(str(mutate(res1)))

        fasta.write("\n")

    fasta.close()
    return bestFit



def mutate(res):
    resMut = ""
    for i in res:
        """mutation of residues 95% of the time we put a normal residue but 5% of the time we put a random residue"""
        if random.random() < 0.95:
            resMut += str(i)
        else:
            resMut += str(random.choice(
                ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))
    return resMut


def tri_to_fasta(res):
    "convert residues to one-letter code"
    res1 = []
    for i in res.triplets:
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





# MAIN # MAIN # MAIN
population = 20
chain_length = 160
fitness = 0
genN = 0
sphereResolution = 0.2 # 0-1 ratio of atoms measured

init_population(population, chain_length)

while fitness < 100:
    #blobulate()
    #fold()

    genN += 1
    fitness = fitness_function()
    print("---------------------------------------------------------", genN, fitness)
