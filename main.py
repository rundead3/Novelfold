"""(x,y,z) = max(x,y,z) - min(x,y,z)"""
import random
import numpy as np
import math
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser
import os
import subprocess
import triplets
from chainLogic import chainLogic
from nicheSpace import nicheSpace


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
    for i in range(0, population):
        chainList.append(chainLogic(i))
        map.add_entry(chainList[i])


    "convert residues into fasta format"
    fasta = open("C:\\Users\\42077\\Omegafold\\randseq.txt", "w")

    for i in range(0, population):
        fasta.write(">" + str(i) + "th chain")
        fasta.write("\n")

        """mutation of residues"""
        fasta.write(str(mutate(triplets.tri_to_fasta(map.get_random()))))

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


# MAIN # MAIN # MAIN
population = 20
chain_length = 160
fitness = 0
genN = 0
sphereResolution = 0.2 # 0-1 ratio of atoms measured

init_population(population, chain_length)
map = nicheSpace()

while fitness < 100:
    #blobulate()
    #fold()

    genN += 1
    fitness = fitness_function()
    print("---------------------------------------------------------", genN, fitness)
    print(map)
