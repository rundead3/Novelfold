import os
import random
import subprocess

import numpy as np

import triplets
from chainLogic import chainLogic

class nicheSpace:

    def __init__(self):

        self.old_archive = None
        self.archive = None
        self.minx = 100
        self.maxx = 0
        self.miny = 100
        self.maxy = 0
        self.boxNo = 8
        self.clear_matrix()


    def __str__(self):
        return str(self.archive)

    def clear_matrix(self):
        self.old_archive = self.archive
        self.archive = np.zeros((self.boxNo, self.boxNo), chainLogic)

    def adjust_range(self, newChains):

        for chain in newChains:
            self.minx = min(self.minx, chain.get_features()[0])
            self.maxx = max(self.maxx, chain.get_features()[0])
            self.miny = min(self.miny, chain.get_features()[1])
            self.maxy = max(self.maxy, chain.get_features()[1])


    def add_entry(self, chain):

        feats = chain.get_features()

        xi = 0
        yi = 0
        for x in np.linspace(self.minx, self.maxx, self.boxNo):
            if feats[0] - x < (self.maxx-self.minx) / self.boxNo:
                for y in np.linspace(self.miny, self.maxy, self.boxNo):
                    if feats[1] - y < (self.maxy - self.miny) / self.boxNo:
                        if self.archive[xi, yi] < chain:
                            self.archive[xi, yi] = chain
                        return
                    yi += 1
                yi = 0
            xi += 1


    def get_random(self):
        choice = 0
        while choice == 0:
            x = random.randrange(0, self.boxNo)
            y = random.randrange(0, self.boxNo)
            choice = self.archive[x, y]
        return choice

    def get_old(self):
        return self.old_archive

    def print_info(self):
        print("X range:[" + str(self.miny) + "," + str(self.maxy) + "]")
        print("Y range:[" + str(self.minx) + "," + str(self.maxx) + "]")


    def write_archive_fastas(self, genNo):

        "convert residues into fasta format"
        fasta = open("C:\\Users\\42077\\Novelfold\\Archive\\fastas\\archive"+str(genNo)+".txt", "w")
        x = 0
        y = 0

        for chainRow in self.archive.tolist():
            for chain in chainRow:
                if chain != 0:
                    fasta.write(">" + str(x) + "_" + str(y) + "_FDB_" + str(int(chain.get_fitness())) + "_" + str(int(100*chain.get_features()[0])) + "_" + str(int(100*chain.get_features()[1])))
                    fasta.write("\n")

                    """mutation of residues"""
                    fasta.write("".join(triplets.tri_to_fasta(chain)))

                    fasta.write("\n")
                y += 1
            x += 1
            y = 0

        fasta.close()


    def fold_archive(self, genNo):
        directory = "C:\\Users\\42077\\Novelfold\\Archive\\pdbs\\"+str(genNo)
        if not os.path.exists(directory):
            os.mkdir(directory)
        subprocess.run(
            "omegafold C:\\Users\\42077\\Novelfold\\Archive\\fastas\\archive"+str(genNo)+".txt "+directory+" --num_cycle 2", shell=True)