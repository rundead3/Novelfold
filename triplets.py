import random


def random_triplet():
    aminochance = [8.25, 5.53, 4.06, 5.45, 1.37, 3.93, 6.75, 7.07, 2.27, 5.96, 9.66, 5.84, 2.42, 3.86, 4.70, 6.56, 5.34,
                   1.08, 2.92, 6.87]
    aminotr = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    floor = 0
    choice = random.randrange(0, 100)
    index = 0
    for amino in aminochance:
        if choice < floor + amino:
            return aminotr[index]
        floor += amino
        index += 1


def random_random_chain(length):
    strSeq = ""
    for i in range(0, length):
        strSeq += str(random.choice(
            ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))
    return strSeq


def random_chain(length):
    chain = ""
    for i in range(length):
        chain += random_triplet()
    return chain


def tri_to_fasta(res):
    "convert residues to one-letter code"
    res1 = []
    for i in res.get_triplets()[0]:
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


def tri_to_fasta_py(res):
    """
    Convert residues to one-letter code and mutate based on provided mutation rates.

    Parameters:
    - res_obj: Object containing residue and mutation information
               (assumed to have get_triplets and get_mutations methods)

    Returns:
    - List of one-letter amino acid codes after considering mutations
    """

    # Map for three-letter to one-letter amino acid codes
    amino_acid_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'
    }

    # Extract the residue triplets and mutation data using provided methods
    triplets = res.get_triplets()[0]
    mutation_data = res.get_mutations()

    mutated_sequence = []
    for idx, i in enumerate(triplets):
        # Default is the original residue
        chosen_residue = amino_acid_map[i]

        # 95% of the time, decide based on mutation data
        if random.random() <= 0.95:
            mutation_info = mutation_data[idx] if idx < len(mutation_data) else {}

            # Decide the mutated residue probabilistically
            rand_value = random.random()
            accumulated_rate = 0
            for mutation, rate in mutation_info.items():
                accumulated_rate += rate
                if rand_value <= accumulated_rate:
                    chosen_residue = mutation[-1]  # Take the last character which is the mutated residue
                    break

        mutated_sequence.append(chosen_residue)

    return mutated_sequence
