"""bio_info.py : Analyse un brin d'adn

Ce programme permet à partir d'un brin d'adn de trouver toutes les protéines
encodées. Dans une premier temps, le brin et son complément sont analysés pour
identifier les gènes qu'ils contiennent et ces gènes sont ensuite traduits en
chaines d'acides aminés qui forment les protéines.

@Date: 28 mars 2024
@Authors: Mathieu Ducharme & Loic Buisson-Fechter
@Contacts: mathieu.ducharme@umontreal.ca & loic.buisson-fechter@umontreal.ca
@Matricules: 20297456 &
"""

from turtle import *
from turtle import pu as penup
from turtle import pd as pendown


DNA = "TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGC\
CGATACCCAGCCAGCCAGCCAGCGACGGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCC\
AGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATA\
TTCAGCCAGCCAGCCAGCGAACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAA\
CTCGACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTG\
CCAGCCAGCATCCCAGCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAA\
CTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTCGTCTGCGTTCGACAGCCA\
GCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGATT\
GCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAA\
GCCAGCCAGCCGTATGCCAGCCAGCATCCCAGCGA"

AMINO_ACIDS_FULLNAME = {
    "UUU": "Phénylalanine",
    "UUC": "Phénylalanine",
    "UUA": "Leucine",
    "UUG": "Leucine",
    "CUU": "Leucine",
    "CUC": "Leucine",
    "CUA": "Leucine",
    "CUG": "Leucine",
    "AUU": "Isoleucine",
    "AUC": "Isoleucine",
    "AUA": "Isoleucine",
    "AUG": "Méthionine (Start)",
    "GUU": "Valine",
    "GUC": "Valine",
    "GUA": "Valine",
    "GUG": "Valine",
    "UCU": "Sérine",
    "UCC": "Sérine",
    "UCA": "Sérine",
    "UCG": "Sérine",
    "CCU": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "ACU": "Thrénine",
    "ACC": "Thrénine",
    "ACA": "Thrénine",
    "ACG": "Thrénine",
    "GCU": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "UAU": "Tyrosine",
    "UAC": "Tyrosine",
    "UAA": "Stop",
    "UAG": "Stop",
    "CAU": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "AAU": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "GAU": "Aspartate",
    "GAC": "Aspartate",
    "GAA": "Glutamate",
    "GAG": "Glutamate",
    "UGU": "Cystéine",
    "UGC": "Cystéine",
    "UGA": "Stop",
    "UGG": "Tryptophane",
    "CGU": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "AGU": "Sérine",
    "AGC": "Sérine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GGU": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine"
}

AMINO_ACIDS_ABBREVIATED = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}


# Équivalent à la spécification 4.1: antisens(brinAdn)
def get_reverse_complement(dna):
    """Génère le brin inversé avec les nucléotides complémentaires

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
    Returns:
        (str): le 2e brin d'adn complémentaire et inversé
    """
    reverse_complement = ''

    for nucleotide in dna.upper():
        if nucleotide == 'A':
            reverse_complement = 'T' + reverse_complement
        elif nucleotide == 'T':
            reverse_complement = 'A' + reverse_complement
        elif nucleotide == 'C':
            reverse_complement = 'G' + reverse_complement
        elif nucleotide == 'G':
            reverse_complement = 'C' + reverse_complement
        else:
            raise ValueError("L'adn contient des caracteres autres que 'A', 'T', 'C' ou 'G'")

    return reverse_complement


def test_get_reverse_complement():
    assert get_reverse_complement('') == ''
    assert get_reverse_complement('TG') == 'CA'
    assert get_reverse_complement('ATCG') == 'CGAT'
    assert get_reverse_complement('CGATTGCCAGCCAGCCAGCCAG') == 'CTGGCTGGCTGGCTGGCAATCG'
    try:
        get_reverse_complement('XYZ$%#')
    except ValueError:
        pass


# Équivalent à la spécification 4.2: trouveDebut(brinAdn)
def get_position_of_starting_codon(dna):
    """Trouve dans un brin d'adn la position de tous les codons de début ('TAC')

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
    Returns:
        (list[int]): liste avec la position de tous les codons de débuts
    """
    starting_codon_position = []
    
    for i in range(len(dna)-2):
        codon = dna[i] + dna[i+1] + dna[i+2]
        if codon == 'TAC':
            starting_codon_position.append(i)

    return starting_codon_position


def test_get_position_of_starting_codon():
    assert get_position_of_starting_codon('') == []
    assert get_position_of_starting_codon('TAC') == [0]
    assert get_position_of_starting_codon('GGGGGGTAC') == [6]
    assert get_position_of_starting_codon('76w898bzfhga') == []
    assert get_position_of_starting_codon('TACTACGGGGTAC') == [0, 3, 10]


# Équivalent à la spécification 4.3: trouveFin(brinAdn)
def get_position_of_ending_codon(dna):
    """Trouve dans un brin d'adn la position de tous les codons de fin ('ATT', 'ATC, 'ACT')

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
    Returns:
        (list[int]): liste avec la position de tous les codons de fin
    """
    starting_codon_position = []

    for i in range(len(dna)-2):
        codon = dna[i] + dna[i+1] + dna[i+2]
        if codon in ['ATT', 'ATC', 'ACT']:
            starting_codon_position.append(i)

    return starting_codon_position


def test_get_position_of_ending_codon():
    assert get_position_of_ending_codon('') == []
    assert get_position_of_ending_codon('ATT') == [0]
    assert get_position_of_ending_codon('ACGCATGCA') == []
    assert get_position_of_ending_codon('GGGGGATTC') == [5]
    assert get_position_of_ending_codon('ATTATCACT') == [0, 3, 6]


# Équivalent à la spécification 4.4: trouveGene(debut, fin)
def get_genes_coordinates(starting_codons, ending_codons):
    """Trouve la position de tous les gènes valides parmis plusieurs couples possibles
    
    Args:
        starting_codons (list): les positions avec le codon de début 'TAC'
        ending_codons (list): les positions avec les codons de fin 'ATT', 'ATC', 'ACT'
    Returns
        (list[tuple]): liste contenant tous les couples de gènes valides
    """
    genes_coordinates = []

    for starting_codon in starting_codons:
        for ending_codon in ending_codons:
            lenght = ending_codon - starting_codon
            if lenght < 0:
                continue
            if lenght % 3 == 0:
                genes_coordinates.append((starting_codon, ending_codon))
                break
        else:
            continue
    
    return genes_coordinates


def test_get_genes_coordinates():
    assert get_genes_coordinates([], []) == []
    assert get_genes_coordinates([0], [3]) == [(0, 3)]
    assert get_genes_coordinates([0], [1, 2, 3, 4, 5]) == [(0, 3)]
    assert get_genes_coordinates([0, 10], [6, 40]) == [(0, 6), (10, 40)]
    assert get_genes_coordinates([2, 11, 30, 41], [5, 33]) == [(2, 5), (30, 33)]


def get_genes(dna, genes_coordinates):
    """Renvoie tous les gènes d'un brin d'adn

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
        genes_coordinates (list[tuple]): positions (début, fin) de tous les gènes valides
    Returns:
        (list[str]): contient tous les gènes sous forme d'ADN
    """
    genes = []

    for gene_coordinate in genes_coordinates:
        start = gene_coordinate[0]
        end = gene_coordinate[1]
        genes.append(dna[start:end+3])

    return genes


def test_get_genes():
    assert get_genes('', []) == []
    assert get_genes('ATCGATACGTCAG', []) == [] 
    assert get_genes('GGTACATTC', [(2, 5)]) == ['TACATT']
    assert get_genes('TACATCGGGGTACACTGGGG', [(0, 3), (10, 13)]) == ['TACATC', 'TACACT']
    assert get_genes('AAABBBCCCDDDEEEFFFGGGHHH', [(2, 6), (0, 21)]) == ['ABBBCCC', 'AAABBBCCCDDDEEEFFFGGGHHH']


# Équivalent à la spécification 4.5: transcrire(brinAdn)
def get_arn_sequence(gene):
    """Transforme une séquence d'ADN (un gène) en ARN

    Args:
        gene (str): gène représenté en adn (suite de A, T, C, G)
    Returns:
        (str): gène représenté en arn (suite de U, A, G, C)
    """
    arn = ''

    for nucleotide in gene:
        if nucleotide == 'A':
            arn += 'U'
        elif nucleotide == 'T':
            arn += 'A'
        elif nucleotide == 'C':
            arn += 'G'
        elif nucleotide == 'G':
            arn += 'C'
        else:
            raise ValueError("L'adn contient des nucléotides non-reconnus")

    return arn


def test_get_arn_sequence():
    assert get_arn_sequence('') == ''
    assert get_arn_sequence('ATCG') == 'UAGC'
    assert get_arn_sequence('AGCCAGCGAA') == 'UCGGUCGCUU'
    assert get_arn_sequence('AGCCGAGTGCCAGC') == 'UCGGCUCACGGUCG'
    try:  # S'assure qu'une erreur est bien détectée
        get_arn_sequence('XYZ$%#')
    except ValueError:
        pass


def get_amino_acids_chain_by_gene(gene):
    """Sépare en codons de 3 nucléotides c'est-à-dire en acides aminés le gène

    Args:
        gene (str): gène représenté en arn (suite de U, A, G, C)
    Returns:
        (list[str]): liste de séquences d'acides aminés sous forme de triplet de nucléotide (ex: ['AUG', 'GGU', 'UAG'])
    """
    amino_acids_chain = []

    for i in range(0, len(gene), 3):
        amino_acids_chain.append(
            gene[i] + gene[i+1] + gene[i+2]
        )

    return amino_acids_chain


def test_get_amino_acids_chain_by_gene():
    assert get_amino_acids_chain_by_gene('') == []
    assert get_amino_acids_chain_by_gene('ATG') == ['ATG']
    assert get_amino_acids_chain_by_gene('X$%') == ['X$%']
    assert get_amino_acids_chain_by_gene('AAGCCA') == ['AAG', 'CCA']
    assert get_amino_acids_chain_by_gene('ATTGCCAGCCAGCCAGCC') == ['ATT', 'GCC', 'AGC', 'CAG', 'CCA', 'GCC']


def get_fullname_amino_acids_chain_by_gene(amino_acids):
    """Trouve la séquence d'acide aminé sous leur nom complet sans le codon d'arret

    Args:
        amino_acids (list[str]): liste de séquences d'acides aminés sous forme de triplet de nucléotide (ex: ['AUG', 'GGU', 'UAG'])
    Returns:
        (list[str]): liste de séquences d'acide aminé sous leur nom complet 
    """

    return list(
        map(
            lambda amino_acid: AMINO_ACIDS_FULLNAME[amino_acid],
            amino_acids[:-1]
        )
    )


def test_get_fullname_amino_acids_chain_by_gene():
    assert get_fullname_amino_acids_chain_by_gene([]) == []
    assert get_fullname_amino_acids_chain_by_gene(['UUU']) == []
    assert get_fullname_amino_acids_chain_by_gene(['CGC', 'UGA']) == ['Arginine']
    assert get_fullname_amino_acids_chain_by_gene(['UUU', 'UAA']) == ['Phénylalanine']
    assert get_fullname_amino_acids_chain_by_gene(['ACA', 'GAA', 'UGC', 'UAG']) == ['Thrénine', 'Glutamate', 'Cystéine']


def get_abbreviated_amino_acids_chain_by_gene(amino_acids):
    """Trouve la séquence d'acide aminé sous leur nom abrégé sans le codon d'arret

    Args:
        amino_acids (list[str]): liste de séquences d'acides aminés sous forme de triplet de nucléotide (ex: ['AUG', 'GGU', 'UAG'])
    Returns:
        (list[str]): liste de séquences d'acide aminé abrégés
    """

    return list(
        map(
            lambda amino_acid: AMINO_ACIDS_ABBREVIATED[amino_acid],
            amino_acids[:-1]
        )
    )


def test_get_abbreviated_amino_acids_chain():
    assert get_abbreviated_amino_acids_chain_by_gene([]) == []
    assert get_abbreviated_amino_acids_chain_by_gene(['UUU']) == []
    assert get_abbreviated_amino_acids_chain_by_gene(['CGC', 'UGA']) == ['R']
    assert get_abbreviated_amino_acids_chain_by_gene(['UUU', 'UAA']) == ['F']
    assert get_abbreviated_amino_acids_chain_by_gene(['ACA', 'GAA', 'UGC', 'UAG']) == ['T', 'E', 'C']


def get_stringify_protein(amino_acids_chain):
    """Transforme une liste d'acides aminés en string

    Args:
        amino_acids_chain (list[str]): liste de séquences d'acide aminé sous leur nom complet
    Returns:
        (str): séquences d'acide aminé sous leur nom complet
    """
    return '-'.join(amino_acids_chain)


def test_get_stringify_protein():
    assert get_stringify_protein([]) == ''
    assert get_stringify_protein(['a', 'b', 'c']) == 'a-b-c'
    assert get_stringify_protein(['a    a', 'b', '?@#']) == 'a    a-b-?@#'
    assert get_stringify_protein(['Valine', 'Tryptophane']) == 'Valine-Tryptophane'
    assert get_stringify_protein(['Méthionine (Start)', 'Glycine', 'Arginine']) == 'Méthionine (Start)-Glycine-Arginine'


def test_get_genes_by_coordinate():
    assert get_genes_by_coordinates([{'dna': '', 'genes_positions': []}]) == []
    assert get_genes_by_coordinates([{'dna': 'TACATT', 'genes_positions': [(0, 3)]}])  == ['TACATT']
    assert get_genes_by_coordinates([{'dna': 'AATACTTTACT', 'genes_positions': [(2,8)]}]) == ['TACTTTACT']
    assert get_genes_by_coordinates([{'dna': 'TACGGGGGGATCGGTACGGATT', 'genes_positions': [(0,9)]}]) == ['TACGGGGGGATC']
    assert get_genes_by_coordinates([{'dna': 'TACGGGGGGATCGGTACGGGATT', 'genes_positions': [(0,9), (14,20)]}]) == ['TACGGGGGGATC', 'TACGGGATT']


def get_genes_coordinates_from_dna(dna):
    """Trouve toutes les positions de gènes sur un brin d'ADN donné

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
    Returns:
        (dict): contient le brin d'adn utilisé et les couples de positions de gènes 
    """
    return {
        "dna": dna,
        "genes_positions": get_genes_coordinates(
            get_position_of_starting_codon(dna),
            get_position_of_ending_codon(dna)
        )
    }


def test_get_genes_coordinate_from_dna():
    assert get_genes_coordinates_from_dna('') == {'dna': '', 'genes_positions': []}
    assert get_genes_coordinates_from_dna('TACATT') == {'dna': 'TACATT', 'genes_positions': [(0, 3)]}
    assert get_genes_coordinates_from_dna('AATACTTTACT') == {'dna': 'AATACTTTACT', 'genes_positions': [(2,8)]}
    assert get_genes_coordinates_from_dna('TACGGGGGGATCGGTACGGATT') == {'dna': 'TACGGGGGGATCGGTACGGATT', 'genes_positions': [(0,9)]}
    assert get_genes_coordinates_from_dna('TACGGGGGGATCGGTACGGGATT') == {'dna': 'TACGGGGGGATCGGTACGGGATT', 'genes_positions': [(0,9), (14,20)]}


# class BranchDto:
#     def __init__(self, dna: str, genes_positions: list[tuple[int, int]]):
#         self.dna = dna
#         self.genes_positions = genes_positions


def get_genes_by_coordinates(branchs): # TODO: make branchs a DTO
    """Trouve la séquence de nucléotide d'un gènes à partir de plusieurs brins
    
    Args:
        branchs (list[{
            'dna': str,
            'genes_positions': list[tuple]
        }]): contient plusieurs brins d'adn et les couples de position de gènes associés
    Returns:
        (list[str]): les gènes représentés sous forme de suite de nucléotide
    """
    genes = []

    for branch in branchs:
        for gene in get_genes(branch['dna'], branch['genes_positions']):
            genes.append(gene)

    return genes


# Équivalent à la spécification 4.6: traduire(brinArn)
def draw_protein(amino_acids_chain):
    """Permet de dessiner une suite d'acides aminés avec Turtle
    Args:
        amino_acids_chain (str): La suite d'acides aminés abrégés qui forme la protéine
    Returns:
        None
    """
    side_length = 15

    for i, amino_acid in enumerate(amino_acids_chain):
        draw_square(i, side_length)
        draw_letter(amino_acid, side_length)


def draw_proteins(amino_acids_chains):
    """Permet de dessiner plusieurs suites d'acides aminés avec Turtle une à la suite de l'autre
    Args:
        amino_acids_chains (list): Les suites d'acides aminés abrégés qui forme les protéines
    Returns:
        None
    """
    x = -112.5; y = 300
    penup(); goto(x, y); pendown()
    for amino_acids_chain in amino_acids_chains:
        draw_protein(amino_acids_chain)
        penup(); rt(90); fd(50); lt(90)
        bk(len(amino_acids_chain) % 15 * 15); pendown()


def draw_letter(amino_acid, side_length):
    """Dessine une lettre pour etre placée à l'intérieur d'un carré

    Args:
        amino_acid (str): Lettre à dessiner
        side_length (int): La taille du carré représentant un acide aminé
    Returns:
        None
    """
    if is_codeboot():
        penup(); bk(side_length/2)
        lt(90); fd(side_length/2); rt(90)
        write(amino_acid)
        lt(90); bk(side_length/2); rt(90)
        fd(side_length/2); pendown()
    else:
        penup(); bk(side_length)
        write('  '+ amino_acid, move=False, align='left', font=('Arial', 8, 'normal'))
        fd(side_length); pendown()


# Équivalent à la spécification 4.7: carre(longueur, nombre)
def draw_square(index, side_length):
    """Dessine un carré avec une lettre à l'intérieur

    Args:
        index (int): Indice du carré à dessiner
        side_length (int): La taille du carré représentant un acide aminé
    Returns:
        None
    """
    if index % 15 == 0 and index != 0:
        penup()
        rt(90); fd(side_length); lt(90)
        bk(side_length * 15)
        pendown()

    for _ in range(4):
        fd(side_length), lt(90)

    penup(); fd(side_length); pendown()


def is_codeboot():
    """Renvoie True si l'interpréteur utilisé est CodeBoot"""
    try:
        import sys
        return False
    except ModuleNotFoundError:
        return True


def run_tests():
    test_get_reverse_complement()
    test_get_position_of_starting_codon()
    test_get_position_of_ending_codon()
    test_get_genes_coordinates()
    test_get_genes()
    test_get_arn_sequence()
    test_get_amino_acids_chain_by_gene()
    test_get_fullname_amino_acids_chain_by_gene()
    test_get_abbreviated_amino_acids_chain()
    test_get_stringify_protein()
    test_get_genes_by_coordinate()
    test_get_genes_coordinate_from_dna()


def main():
    clear(800, 600) if is_codeboot() else speed(0)

    dna_strands = [
        DNA,
        get_reverse_complement(DNA)
    ]

    genes_coordinates = list(map(get_genes_coordinates_from_dna, dna_strands))
    genes = get_genes_by_coordinates(genes_coordinates)

    amino_acids_chains_fullname = []
    amino_acids_chains_shorten = []

    for gene in genes:
        amino_acids = get_amino_acids_chain_by_gene(get_arn_sequence(gene))
        amino_acids_chains_fullname.append(
            get_fullname_amino_acids_chain_by_gene(amino_acids)
        )
        amino_acids_chains_shorten.append(
            get_abbreviated_amino_acids_chain_by_gene(amino_acids)
        )

    for amino_acids_chain in amino_acids_chains_fullname:
        stringify_protein = get_stringify_protein(amino_acids_chain)
        print(stringify_protein + '\n')

    draw_proteins(amino_acids_chains_shorten)

    mainloop() if not is_codeboot() else None


if __name__ == '__main__':
    run_tests()
    main()
