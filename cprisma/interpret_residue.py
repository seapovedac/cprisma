### Function to interpret standard three letters to standard one letter aminoacid name

def interpret_residue(am):

    global residue

    if am == 'ALA':
        residue='A'
    if am == 'ARG':
        residue='R'
    if am == 'ASN':
        residue='N'
    if am == 'ASP':
        residue='D'
    if am == 'CYS':
        residue='C'
    if am == 'GLN':
        residue='Q'
    if am == 'GLU':
        residue='E'
    if am == 'GLY':
        residue='G'
    if am == 'HIS':
        residue='H'
    if am == 'ILE':
        residue='I'
    if am == 'LEU':
        residue='L'
    if am == 'LYS':
        residue='K'
    if am == 'MET':
        residue='M'
    if am == 'PHE':
        residue='F'
    if am == 'PRO':
        residue='P'
    if am == 'SER':
        residue='S'
    if am == 'THR':
        residue='T'
    if am == 'TRP':
        residue='W'
    if am == 'TYR':
        residue='Y'
    if am == 'VAL':
        residue='V'
    if am == 'UNK':
        residue='Z'
    if am == 'XXX':
        residue='X'
    if am == '000':
        residue='0'

    return residue

def interpret_residue3(am):

    global residue

    if am == 'A':
        residue='ALA'
    if am == 'R':
        residue='ARG'
    if am == 'N':
        residue='ASN'
    if am == 'D':
        residue='ASP'
    if am == 'C':
        residue='CYS'
    if am == 'Q':
        residue='GLN'
    if am == 'E':
        residue='GLU'
    if am == 'G':
        residue='GLY'
    if am == 'H':
        residue='HIS'
    if am == 'I':
        residue='ILE'
    if am == 'L':
        residue='LEU'
    if am == 'K':
        residue='LYS'
    if am == 'M':
        residue='MET'
    if am == 'F':
        residue='PHE'
    if am == 'P':
        residue='PRO'
    if am == 'S':
        residue='SER'
    if am == 'T':
        residue='THR'
    if am == 'W':
        residue='TRP'
    if am == 'Y':
        residue='TYR'
    if am == 'V':
        residue='VAL'
    if am == 'Z':
        residue='UNK'
    if am == 'X':
        residue='XXX'
    if am == '0':
        residue='000'

    return residue
