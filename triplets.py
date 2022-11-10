

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
