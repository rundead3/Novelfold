import numpy as np
from Bio.PDB.PDBParser import PDBParser
import math


def get_bounding_sphere(coords):
    maxDist = 0

    resolution = 0.1  # every 10th atom

    for i in range(0, len(coords), int(len(coords) * resolution)):
        for j in range(i, len(coords), int(len(coords) * resolution)):
            maxDist = max(maxDist, math.dist(coords[i], coords[j]))

    return maxDist


def get_gyration_radius(coords):
    # get the center of mass
    center = np.mean(coords, axis=0)

    # get the squared distances from the center of mass squared
    squared_distances = np.sum((coords - center) ** 2, axis=1)

    # get the gyration radius
    gyration_radius = np.sqrt(np.mean(squared_distances))

    return gyration_radius


class chainLogic:

    def __init__(self, i):
        self.gRadius = None
        self.confidence = None
        self.chainLength = None
        self.triplets = []
        self.index = i
        # self.blobulation = self.load_blobulator()
        self.load_omegafold()
        self.secondaryStruct = self.load_dssp()

    def __str__(self):
        msg = "Conf:" + str(self.confidence) + " Density:" + str(self.gRadius) + " secondary:" + str(
            self.secondaryStruct) + " Lengh:" + str(self.chainLength)
        return msg

    def __repr__(self):
        return str(int(self.confidence)) + "|" + str(int(self.gRadius))

    def __gt__(self, other):
        if other == 0:
            return True
        return self.get_fitness() > other.get_fitness()

    def load_blobulator(self):
        csv_file = open(
            "C:\\Users\\42077\\Omegaforl\\blobulator-main\\Batch" + str(self.index) + "th.csv",
            "r")
        csv_list = csv_file.readlines()
        csv_file.close()

        csv_list.pop(0)
        nichx = float(0)
        nichy = float(0)
        for i in csv_list:
            nichx += float(i.split(',')[5])
            nichy += float(i.split(',')[6])

        return nichx / len(csv_list), nichy / len(csv_list)

    def load_omegafold(self):
        parser = PDBParser()
        structure = parser.get_structure(str(self.index) + "th chain",
                                         "C:\\Users\\42077\\Omegaforl\\res1\\" + str(self.index) + "th_chain.pdb")

        for model in structure:
            for chain in model:
                a = 0
                b = 0
                l = 0
                res = []
                coords = []
                for residue in chain:
                    l += 1
                    res.append(residue.get_resname())
                    for atom in residue:
                        b += 1
                        a += atom.bfactor
                        coords.append(atom.get_coord())

                avg = (a / b)
                self.chainLength = l
                self.confidence = avg

                self.triplets.append(res)

                ##self.denseScore = 1 / (get_bounding_sphere(coords) / math.pow(self.chainLength, 1/3))

                self.gRadius = get_gyration_radius(coords)

    def load_dssp(self):
        file = open("C:\\Users\\42077\\Omegaforl\\res1\\dssps.txt", "r")
        linelist = file.readlines()
        file.close()

        ourline = linelist[self.index]
        alphas = 0
        betas = 0
        for i in range(0, self.chainLength):
            if ourline[i] == "H":
                alphas += 1
            elif ourline[i] == "E":
                betas += 1
        return alphas / self.chainLength, betas / self.chainLength

    def get_fitness(self):
        return self.confidence/self.gRadius

    def get_features(self):
        return self.secondaryStruct[0], self.secondaryStruct[1]

    def get_triplets(self):
        return self.triplets

    def survivable(self, cutoff):
        # return self.confidence > cutoff
        return True
