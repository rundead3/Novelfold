from Bio.PDB.PDBParser import PDBParser
import math


def get_bounding_sphere(coords):
    maxDist = 0

    resolution = 0.1 # every 10th atom

    for i in range(0, len(coords), int(len(coords) * resolution)):
        for j in range(i, len(coords), int(len(coords) * resolution)):
            maxDist = max(maxDist, math.dist(coords[i], coords[j]))

    return maxDist


class chainLogic:

    def __init__(self, i):
        self.denseScore = []
        self.fitScore = None
        self.chainLength = None
        self.triplets = []
        self.index = i
        self.blobulation = self.load_blobulator()
        self.load_omegafold()

    def __str__(self):
        msg = "Fit:" + str(self.fitScore) + " Density:" + str(self.denseScore) + " blob:" + str(
            self.blobulation) + " Lengh:" + str(self.chainLength)
        return msg

    def load_blobulator(self):
        csv_file = open(
            "C:\\Users\\42077\\Omegaforl\\blobulator\\blobulator-main\\Batch\\Outputs\\" + str(self.index) + "th.csv",
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
                                         "C:\\Users\\42077\\Omegaforl\\res1\\" + str(self.index) + "th chain.pdb")
        for model in structure:
            for chain in model:
                a = 0
                b = 0
                res = []
                coords = []
                for residue in chain:
                    res.append(residue.get_resname())
                    for atom in residue:
                        b += 1
                        a += atom.bfactor
                        coords.append(atom.get_coord())

                fitness = (a / b)
                self.chainLength = b
                self.fitScore = fitness

                self.triplets.append(res)

                self.denseScore.append(get_bounding_sphere(coords))
