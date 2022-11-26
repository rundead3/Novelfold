"""(x,y,z) = max(x,y,z) - min(x,y,z)"""
import random
import os
import subprocess
import triplets
from chainLogic import chainLogic
from nicheSpace import nicheSpace
import blosum as bl

def init_population(popSize, length):
    fasta = open("C:\\Users\\Rundead\\Omegafold\\randseq.txt", "w")

    for i in range(0, popSize):
        fasta.write(">" + str(i) + "th_chain")
        fasta.write("\n")

        """mutation of residues"""
        fasta.write(str(triplets.random_chain(length)))

        fasta.write("\n")


def random_chain(length):
    strSeq = ""
    for i in range(0, length):
        strSeq += str(random.choice(
            ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))
    return strSeq


def blobulate():
    """COMPUTE BLOBS"""
    #os.chdir(blobpath)
    subprocess.run(
        "python compute_blobs.py --fasta C:\\Users\\Rundead\\Omegafold\\randseq.txt --cutoff 0.4 --minBlob 4 --oname C:\\Users\\Rundead\\Omegaforl\\blobulator-main",
        shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


############################################
##def netsurf(fasta):
##    nsp3 = biolib.load('DTU/NetSurfP-3')
##    nsp3_results = nsp3.cli(args='-i' + str(fasta))
##    nsp3_results.save_files("C:\\Users\\Rundead\\Omegafold\\netsurf")


def feature_extraction():
    """extract features from pdb files"""
    os.chdir("C:\\Users\\Rundead\\Videos\\PyDSSP-0.9.0\\scripts")
    subprocess.run("python pydssp -i C:\\Users\\Rundead\\Omegaforl\\res1\\ -o C:\\Users\\Rundead\\Omegaforl\\res1\\")

def fold():
    "START OMEGAFOLD"
    os.chdir("C:\\Users\\Rundead\\OmegaFold")
    subprocess.run("omegafold C:\\Users\\Rundead\\Omegafold\\randseq.txt C:\\Users\\Rundead\\Omegaforl\\res1 --num_cycle 1",
                   shell=True)
    "END OMEGAFOLD"


def dssp():
    pdbs = ""
    for i in range(0, population):
        pdbs += " C:\\Users\\Rundead\\Omegaforl\\res1\\"+str(i) + "th_chain.pdb"
    output = " -o C:\\Users\\Rundead\\Omegaforl\\res1\\dssps.txt"
    subprocess.run("python C:\\Users\\Rundead\\Omegaforl\\PyDSSP-0.9.0\\scripts\\pydssp"+pdbs+output)



def next_gen():
    # Create individuals
    newChainList = []
    for i in range(0, population):
        newPeep = chainLogic(i)
        if newPeep.survivable(confidence_cutoff):
            newChainList.append(newPeep)
        else:
            print(newPeep)

    map.adjust_range(newChainList)
    map.clear_matrix()

    # new
    for i in range(0, len(newChainList)):
        map.add_entry(newChainList[i])
    # old
    for chainRow in map.get_old().tolist():
        for chain in chainRow:
            if chain != 0:
                map.add_entry(chain)


def new_population():
    "convert residues into fasta format"
    fasta = open("C:\\Users\\Rundead\\Omegafold\\randseq.txt", "w")

    for i in range(0, population):
        fasta.write(">" + str(i) + "th_chain")
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
    dic = bl.BLOSUM(45)
    for i in range(0, len(res1)):
        """mutation of residues 95% of the time we put a normal residue but 5% of the time we put a random residue"""
        ##if random.random() < mutChance:
            ##val = dic(str(i)+ str(for i in ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
        if random.random() < mutChance:
            resMut += str(triplets.random_triplet())
            #resMut += str(random.choice(
            #    ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))
        else:
            if random.random() < crossChance:
                read1 = not read1

            if read1:
                resMut += str(res1[i])
            else:
                resMut += str(res2[i])
    return resMut


# MAIN # MAIN # MAIN

## settings ##
population = 20
confidence_cutoff = 0
chain_length = 144
generations = 100000


# init
fitness = 0
genN = 0
sphereResolution = 0.4  # ratio of atoms measured

init_population(population, chain_length)
map = nicheSpace()

while genN < generations:
    #blobulate()
    fold()
    dssp()

    genN += 1
    next_gen()

    new_population()

    print(map)
    map.print_info()

    if genN % 25 == 0 and genN != 0:
        map.write_archive_fastas(genN)
        map.fold_archive(genN)

    print("----------------------------------", genN, "----------------------------------")
