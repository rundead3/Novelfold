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
        shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


def fold():
    "START OMEGAFOLD"
    os.chdir("C:\\Users\\42077\\OmegaFold")
    subprocess.run("omegafold C:\\Users\\42077\\Omegafold\\randseq.txt C:\\Users\\42077\\Omegaforl\\res1 --num_cycle 1",
                   shell=True)
    "END OMEGAFOLD"


def next_gen():
    # Create individuals
    newChainList = []
    for i in range(0, population):
        newChainList.append(chainLogic(i))

    map.adjust_range(newChainList)
    map.clear_matrix()

    # new
    for i in range(0, population):
        map.add_entry(newChainList[i])
    # old
    for chainRow in map.get_old().tolist():
        for chain in chainRow:
            if chain != 0:
                map.add_entry(chain)


def new_population():
    "convert residues into fasta format"
    fasta = open("C:\\Users\\42077\\Omegafold\\randseq.txt", "w")

    for i in range(0, population):
        fasta.write(">" + str(i) + "th chain")
        fasta.write("\n")

        """mutation of residues"""
        fasta.write(str(mutate(triplets.tri_to_fasta(map.get_random()), triplets.tri_to_fasta(map.get_random()))))

        fasta.write("\n")

    fasta.close()


def mutate(res1, res2):
    crossChance = 0.1
    mutChance = 0.05
    read1 = True
    resMut = ""
    for i in range(0, len(res1)):
        """mutation of residues 95% of the time we put a normal residue but 5% of the time we put a random residue"""
        if random.random() < mutChance:
            resMut += str(random.choice(
                ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))
        else:
            if random.random() < crossChance:
                read1 = not read1

            if read1:
                resMut += str(res1[i])
            else:
                resMut += str(res2[i])
    return resMut


# MAIN # MAIN # MAIN
population = 20
chain_length = 160
generations = 500

fitness = 0
genN = 0
sphereResolution = 0.2  # ratio of atoms measured

init_population(population, chain_length)
map = nicheSpace()

while genN < generations:
    blobulate()
    fold()

    genN += 1
    next_gen()

    new_population()

    print("----------------------------------", genN, "----------------------------------")
    print(map)
    print(map.print_info())

    if genN % 100 == 1:
        map.write_archive_fastas(genN)
