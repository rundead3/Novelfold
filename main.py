"""(x,y,z) = max(x,y,z) - min(x,y,z)"""
import random
from random import randrange
import os
import subprocess

import numpy as np

import triplets
from chainLogic import chainLogic
from nicheSpace import nicheSpace


def init_population(popSize, length):
    fasta = open("C:\\Users\\42077\\Omegafold\\randseq.txt", "w")

    for i in range(0, popSize):
        fasta.write(">" + str(i) + "th_chain")
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
    subprocess.run("omegafold C:\\Users\\42077\\Omegafold\\randseq.txt C:\\Users\\42077\\Omegaforl\\res1 --num_cycle 1 --subbatch_size 256",
                   shell=True)
    "END OMEGAFOLD"


def dssp():
    pdbs = ""
    for i in range(0, population):
        pdbs += " C:\\Users\\42077\\Omegaforl\\res1\\"+str(i) + "th_chain.pdb"
    output = " -o C:\\Users\\42077\\Omegaforl\\res1\\dssps.txt"
    subprocess.run("python C:\\Users\\42077\\Omegaforl\\PyDSSP\\scripts\\pydssp"+pdbs+output)



def next_gen():
    # Create individuals
    newChainList = []
    for i in range(0, population):
        newPeep = chainLogic(i)
        if newPeep.survivable(confidence_cutoff):
            newChainList.append(newPeep)
        else:
            print(newPeep)

    maps[epoch].adjust_range(newChainList)
    maps[epoch].clear_matrix()

    # new
    for i in range(0, len(newChainList)):
        maps[epoch].add_entry(newChainList[i])

    # old only if not new epoch
    if genN % epochGens != 0:
        for chainRow in maps[epoch].get_old().tolist():
            for chain in chainRow:
                if chain != 0:
                    maps[epoch].add_entry(chain)


def new_population():
    "convert residues into fasta format"
    fasta = open("C:\\Users\\42077\\Omegafold\\randseq.txt", "w")

    for i in range(0, population):
        fasta.write(">" + str(i) + "th_chain")
        fasta.write("\n")

        """mutation of residues"""
        newChain = str(mutate(triplets.tri_to_fasta(maps[randrange(epoch+1)].get_random()), triplets.tri_to_fasta(maps[randrange(epoch+1)].get_random())))
        while len(newChain) < chain_length:
            newChain += newChain[::-1]
        fasta.write(newChain)

        fasta.write("\n")

    fasta.close()


def mutate(res1, res2):
    crossChance = 0.1
    mutChance = 0.05
    read1 = True
    resMut = ""
    for i in range(0, min(len(res1), len(res2))):
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
confidence_cutoff = 28
chain_length = 64
generations = 2000
epochGens = 400
epoch = 0

fitness = 0
genN = 0
sphereResolution = 0.2  # ratio of atoms measured

init_population(population, chain_length)
maps = []
maps.append(nicheSpace())

while genN < generations:
    #blobulate()
    fold()
    dssp()

    genN += 1
    next_gen()

    new_population()

    print(maps[epoch])
    maps[epoch].print_info()

    if genN % 50 == 0:
        maps[epoch].write_archive_fastas(genN)
        maps[epoch].fold_archive(genN)

    if genN % epochGens == 0:
        chain_length *= 2
        epoch += 1
        maps.append(nicheSpace())

    print("----------------------------------", genN, "----------------------------------")
