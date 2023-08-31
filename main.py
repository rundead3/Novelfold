"""(x,y,z) = max(x,y,z) - min(x,y,z)"""
import random
import os
import subprocess
import triplets
import config
from chainLogic import chainLogic
from nicheSpace import nicheSpace
import freesasa
import pickle
import math


def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')


def init_population(popSize, length):
    fasta = open(config.get_new_fastas_path(), "w")

    for i in range(0, popSize):
        fasta.write(">" + str(i) + "th_chain")
        fasta.write("\n")

        """mutation of residues"""
        fasta.write(str(triplets.random_chain(length)))

        fasta.write("\n")


def blobulate():
    """COMPUTE BLOBS"""
    subprocess.run(
        "python compute_blobs.py --fasta " + config.get_new_fastas_path() + " --cutoff 0.4 --minBlob 4 --oname " + config.get_blobulator_path(),
        shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


############################################
##def netsurf(fasta):
##    nsp3 = biolib.load('DTU/NetSurfP-3')
##    nsp3_results = nsp3.cli(args='-i' + str(fasta))
##    nsp3_results.save_files("C:\\Users\\Rundead\\Omegafold\\netsurf")

def fold():
    "START OMEGAFOLD"
    os.chdir(config.get_omegafold_path())
    subprocess.run("omegafold " + config.get_new_fastas_path() + " " + config.get_pdbs_path() + " --num_cycle 1",
                   shell=True)
    "END OMEGAFOLD"


def dssp():
    pdbs = ""
    for i in range(0, population):
        pdbs += " " + config.get_pdbs_path() + str(i) + "th_chain.pdb"
    output = " -o " + config.get_dssp_output_path()
    subprocess.run("python " + config.get_dssp_path() + pdbs + output)


def normalize_to_probabilities(data):
    """Normalize free energy differences to positive scores between 0 and 1."""

    # Convert free energy differences to positive scores
    score_data = {key: math.exp(-value) for key, value in data.items()}

    # Normalize scores for each index
    index_data = {}
    for key, score in score_data.items():
        index = key[:-1]
        if index not in index_data:
            index_data[index] = 0
        index_data[index] += score

    normalized_data = {}
    for key, score in score_data.items():
        index = key[:-1]
        normalized_data[key] = score / index_data[index]

    return normalized_data


def convert_to_normalized_probabilities(input_file, output_file):
    """Convert the values in the input file to normalized probabilities and save to output file."""

    data = {}
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            key = parts[0]
            value = float(parts[1])
            data[key] = value

    normalized_data = normalize_to_probabilities(data)

    with open(output_file, 'w') as file:
        for key, value in normalized_data.items():
            file.write(f"{key} {value:.4f}\n")


def pythia():
    batch_file_path = r"C:\Users\Rundead\Omegaforl\res1\new.bat"
    subprocess.run(batch_file_path, shell=True)
    # Run the convert to normalized probabilities function on all the files just created
    for i in range(0, population):
        input_file = r"C:\Users\Rundead\Omegaforl\res1\\" + str(i) + "th_chain_pred_mask.txt"
        output_file = r"C:\Users\Rundead\Omegaforl\res1\\" + str(i) + "th_chain.txt"
        convert_to_normalized_probabilities(input_file, output_file)


def foldx():
    pdbsi = []
    for i in range(0, population):
        pdbsi.append(str(i) + "th_chain.pdb")
    os.chdir(config.get_pdbs_path())
    for i in pdbsi:
        result = subprocess.run("Foldx_5.exe -cStability --pdb=" + str(i), shell=True, capture_output=True, text=True)
        # Get the captured output
        output = result.stdout

        # Process the output to extract the desired information
        # Here's an example to extract the stability value
        start_index = output.find("Total          = ")
        end_index = output.find("\n", start_index)
        stability = output[start_index + len("Total          = "):end_index]

        print("Stability:", stability)


def sasa():
    pdbs = []
    for i in range(0, population):
        pdbs.append(config.get_pdbs_path() + str(i) + "th_chain.pdb")
    print(pdbs[0])
    print(pdbs)
    for pdb in pdbs:
        structure = freesasa.Structure(pdb)
        result = freesasa.calc(structure)
        area_classes = freesasa.classifyResults(result, structure)

        print("Total : %.2f A2" % result.totalArea())
        for key in area_classes:
            print(key, ": %.2f A2" % area_classes[key])


def calculate_sasa():
    pdbs = [config.get_pdbs_path() + str(i) + "th_chain.pdb" for i in range(population)]
    sasa_values = []
    for pdb in pdbs:
        structure = freesasa.Structure(pdb)
        result = freesasa.calc(structure)
        total_area = result.totalArea()
        sasa_values.append(total_area)
    return sasa_values


def next_gen():
    # Create individuals
    global index_global
    newChainList = []
    for i in range(0, population):
        newPeep = chainLogic(index_global, i)
        index_global += 1
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
    map.get_new_entries()


def new_population():
    "convert residues into fasta format"
    fasta = open(config.get_new_fastas_path(), "w")

    for i in range(0, population):
        fasta.write(">" + str(i) + "th_chain")
        fasta.write("\n")

        """mutation of residues"""

        fasta.write(str(mutate(triplets.tri_to_fasta_py(map.get_random()), triplets.tri_to_fasta_py(map.get_random()))))

        fasta.write("\n")

    fasta.close()


def mutate(res1, res2):
    crossChance = 0.05
    mutChance = 0.0
    read1 = True
    resMut = ""
    for i in range(0, len(res1)):
        """mutation of residues 95% of the time we put a normal residue but 5% of the time we put a random residue"""
        if random.random() < mutChance:
            resMut += str(triplets.random_triplet())
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
population = 8
confidence_cutoff = 0
chain_length = 100
generations = 100000000000

# init
fitness = 0
genN = 0

index_global = 0
init_population(population, chain_length)
map = nicheSpace()

while genN < generations:
    # blobulate()
    fold()
    pythia()
    dssp()
    # foldx()
    sasa_values = calculate_sasa()
    with open('sasa_values.pkl', 'wb') as f:
        pickle.dump(sasa_values, f)

    genN += 1

    next_gen()

    new_population()

    print(map)
    map.print_info()

    if genN % 25 == 0 and genN != 0:
        map.write_archive_fastas(genN)
        map.fold_archive(genN)
        clear_console()

    print("----------------------------------", genN, "----------------------------------")
