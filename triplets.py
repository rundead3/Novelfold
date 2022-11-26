import random
def random_triplet():
    aminochance=[8.25,5.53,4.06,5.45,1.37,3.93,6.75,7.07,2.27,5.96,9.66,5.84,2.42,3.86,4.70,6.56,5.34,1.08,2.92,6.87]
    aminotr=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    floor = 0
    choice = random.randrange(0,100)
    index = 0
    for amino in aminochance:
        if choice < floor + amino:
            return aminotr[index]
        floor += amino
        index += 1

def random_chain(length):
    chain = ""
    for i in range(length):
        chain += random_triplet()
    return chain

def random_random_chain(length):
    strSeq = ""
    for i in range(0, length):
        strSeq += str(random.choice(
            ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']))
    return strSeq

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
