import numpy as np
from Bio.PDB.PDBParser import PDBParser
import math
import pickle

def get_bounding_sphere(coords):
    maxDist = 0

    resolution = 0.1  # every 10th atom

    for i in range(0, len(coords), int(len(coords) * resolution)):
        for j in range(i, len(coords), int(len(coords) * resolution)):
            maxDist = max(maxDist, math.dist(coords[i], coords[j]))

    return maxDist


def get_start_end_distance(coords):
    return math.dist(coords[0], coords[-1])


def get_rmsd(chain1, chain2):
    # get the coordinates of the two chains
    coords1 = np.array([chain1.get_coord() for atom in chain1.get_atoms()])
    coords2 = np.array([atom.get_coord() for atom in chain2.get_atoms()])

    # get the rmsd
    rmsd = np.sqrt(np.mean((coords1 - coords2) ** 2))

    return rmsd


def get_gyration_radius(coords):
    # get the center of mass
    center = np.mean(coords, axis=0)

    # get the squared distances from the center of mass squared
    squared_distances = np.sum((coords - center) ** 2, axis=1)

    # get the gyration radius
    gyration_radius = np.sqrt(np.mean(squared_distances))

    return gyration_radius


class chainLogic:

    def __init__(self, gi, i):
        self.global_index = gi
        self.denseScore = None
        self.confidence = None
        self.chainLength = None
        self.triplets = []
        self.index = i
        # self.blobulation = self.load_blobulator()
        self.sasa = self.load_sasa()
        self.load_omegafold()
        self.secondaryStruct = self.load_dssp()


    def __str__(self):
        msg = "Conf:" + str(self.confidence) + " Fitness:" + str(self.get_fitness()) + " Secondary:" + str(
            self.secondaryStruct) + " Length:" + str(self.chainLength)
        return msg

    def __repr__(self):
        return str(int(self.confidence)) + "|" + str(int(self.denseScore))

    def __gt__(self, other):
        if other == 0:
            return True
        return self.get_fitness() < other.get_fitness()

    def __lt__(self, other):
        if other == 0:
            return True
        return self.get_fitness() > other.get_fitness()
    def load_blobulator(self):
        csv_file = open(
            "C:\\Users\\Rundead\\Omegaforl\\blobulator-main\\Batch" + str(self.index) + "th.csv",
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
                                         "C:\\Users\\Rundead\\Omegaforl\\res1\\" + str(self.index) + "th_chain.pdb")

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

                #self.denseScore = get_gyration_radius(coords) * get_start_end_distance(coords)
                self.denseScore = self.sasa

    def load_dssp(self):
        file = open("C:\\Users\\Rundead\\Omegaforl\\res1\\dssps.txt", "r")
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

    def load_sasa(self):
        with open('sasa_values.pkl', 'rb') as f:
            sasa_values = pickle.load(f)
        return sasa_values[self.index]

    def get_fitness(self):
        return self.denseScore

    def get_features(self):
        return self.secondaryStruct[0], self.secondaryStruct[1]

    def get_triplets(self):
        return self.triplets

    def survivable(self, cutoff):
        # return self.confidence > cutoff
        return True
    def get_index(self):
        return self.global_index
