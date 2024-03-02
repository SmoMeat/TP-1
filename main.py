# TP1 - IFT 1015

# le premier brin a 3 gènes
# et le deuxieme en a un seul

# 1) générer le brin complémentaire  OK!
# 2) trouver les gènes (commence par TAC et fini par ATT, ATC ou ACT)
# * le premier brin se lit de gauche à droite, le second de droit à gauche
# 3) transformer les gènes en ARN (A->U, T->A, C->G, G->C)
# 4) transformer la suite d'ARN en acide-aminé
# 5) traduire la suite d'acide-aminé en abréviation
# 6) afficher avec turtle

from bioInfo import adn as dna
from bioInfo import codons_aa as amino_acids_template
from bioInfo import lettreAa as amino_acids_letters
from turtle import pu as penup
from turtle import pd as pendown
from turtle import *
import math


def get_reverse_complement(dna):
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
            raise ValueError ("L'adn contient des caracteres autres que 'A', 'T', 'C' ou 'G'")

    return reverse_complement
    
def get_position_of_starting_codon(dna):
    """Trouve dans un brin d'adn la position de tous les codons de début ('TAC')

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
    Returns:
        starting_codon_position (list[int]): liste avec la position de tous les codons de débuts
    """
    starting_codon_position = []
    for i in range(len(dna)-2):
        codon = dna[i] + dna[i+1] + dna[i+2]
        if codon == 'TAC':
            starting_codon_position.append(i)
    return starting_codon_position

def get_position_of_ending_codon(dna):
    """Trouve dans un brin d'adn la position de tous les codons de fin ('ATT', 'ATC, 'ACT')

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
    Returns:
        ending_codon_position (list[int]): liste avec la position de tous les codons de fin
    """
    starting_codon_position = []
    for i in range(len(dna)-2):
        codon = dna[i] + dna[i+1] + dna[i+2]
        if codon in ['ATT', 'ATC', 'ACT']:
            starting_codon_position.append(i)
    return starting_codon_position

def get_genes_coordinates(starting_codons, ending_codons):
    """Trouve la position de tous les gènes valides parmis plusieurs couples possibles
    
    Args:
        starting_codons (list): les positions avec le codon de début 'TAC'
        ending_codons (list): les positions avec les codons de fin 'ATT', 'ATC', 'ACT'
    Returns
        genes_coordinates (list[tuple]): liste contenant tous les couples de gènes valides
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

def get_genes(dna, genes_coordinates):
    """Renvoie tous les gènes d'un brin d'adn

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
        genes_coordinates (list[tuple]): positions (début, fin) de tous les gènes valides
    Returns:
        genes (list[str]): contient tous les gènes sous forme d'ADN
    """
    genes = []
    for gene_coordinate in genes_coordinates:
        start = gene_coordinate[0]
        end = gene_coordinate[1]
        genes.append(dna[start:end+3])
    return genes

def get_arn_sequence(gene):
    arn = ''
    for nucleotide in gene:
        if nucleotide == 'A':  arn += 'U'
        elif nucleotide == 'T': arn += 'A'
        elif nucleotide == 'C': arn += 'G'
        elif nucleotide == 'G': arn += 'C'
        else: raise ValueError
    return arn

def get_amino_acids(gene):
    amino_acids_chain = []
    for i in range(0, len(gene), 3):
        amino_acids_chain.append(
            gene[i] + gene[i+1] + gene[i+2]
        )
    return amino_acids_chain

def get_amino_acids_chain(gene):
    amino_acids = get_amino_acids(gene)
    amino_acids_chain = []
    for amino_acid in amino_acids[:-1]:
        amino_acids_chain.append(amino_acids_template[amino_acid])
    return amino_acids_chain

def get_abbreviated_amino_acids_chain(gene):
    amino_acids = get_amino_acids(gene)
    amino_acids_names = []
    for amino_acid in amino_acids[:-1]:
        amino_acids_names.append(amino_acids_letters[amino_acid])
    return amino_acids_names

def get_stringify_protein(amino_acids_chain):
    return '-'.join(amino_acids_chain)

def get_stringify_protein_by_gene(gene):
    amino_acids_chain = get_amino_acids_chain(gene)
    return '-'.join(amino_acids_chain)

def draw_protein(x, y, amino_acids_chain, side_length=15):
    """Permet de dessiner une suite d'acides aminés avec Turtle
    Args:
        x (float): Position horizontale où la protéine doit etre dessiné
        y (float): Position verticale où la protéine doit etre dessiné
        amino_acids_chain (str): La suite d'acides aminés abrégés qui forme la protéine
        side_length (int): La taille des carrés représentant un acide aminé (défautl=15)
    Returns:
        None
    """
    _x, _y = x, y
    for i, amino_acid in enumerate(amino_acids_chain):
        if i % 15 == 0 and i != 0:
            _y -= side_length
            _x = x
        draw_square(_x, _y, amino_acid, side_length)
        _x += side_length

def draw_proteins(amino_acids_chains):
    """Permet de dessiner plusieurs suites d'acides aminés avec Turtle une à la suite de l'autre
    Args:
        amino_acids_chains (list): Les suites d'acides aminés abrégés qui forme les protéines
    Returns:
        None
    """
    x = -112.5 ; y = 300
    for amino_acids_chain in amino_acids_chains:
        draw_protein(x, y, amino_acids_chain)
        if is_codeboot(): y -= 15 * (len(amino_acids_chain) // 15) + 50
        else: y = ycor() - 50
        
def get_genes_by_coordinate(branchs):
    """Trouve la séquence de nucléotide d'un gènes à partir de plusieurs brins
    
    Args:
        branchs (list[{
            'dna': str,
            'genes_positions': list[tuple]
        }]): contient plusieurs brins d'adn et les couples de position de gènes associés
    Returns:
        genes (list[str]): les gènes représentés sous forme de suite de nucléotide
    """
    genes = []
    for branch in branchs:
        for gene in get_genes(branch['dna'], branch['genes_positions']):
            genes.append(gene)

    return genes


def get_genes_coordinate_from_dna(dna):
    """Trouve toutes les positions de gènes sur un brin d'ADN donné

    Args:
        dna (str): l'adn sous forme de string (suite de A, T, C, G)
    Returns:
        ? (dict): contient le brin d'adn utilisé et les couples de positions de gènes 
    """
    return {
        "dna": dna,
        "genes_positions": get_genes_coordinates(
            get_position_of_starting_codon(dna),
            get_position_of_ending_codon(dna)
        )
    }

def draw_square(x, y, amino_acid, side_length=15):
    """Dessine un carré avec une lettre à l'intérieur

    Args:
        x (float): Position horizontale où le carré doit etre dessiné
        y (float): Position verticale où le carré doit etre dessiné
        amino_acid (str): La lettre représentant l'acide aminé
        side_length (int): La taille du carré représentant un acide aminé (défautl=15)
    Returns:
        None
    """
    if is_codeboot():
        penup(), goto(x+side_length/2, y+side_length/2), write(amino_acid)
        goto(x,y), pendown()
        for side in range(4):
            fd(side_length), lt(90)
    else:
        penup(), goto(x, y), pendown()
        write('  '+ amino_acid, move=False, align='left', font=('Arial', 8, 'normal'))
        for side in range(4):
            fd(side_length), lt(90)

def is_codeboot():
    """Renvoie True si l'interpréteur utilisé est CodeBoot"""
    try:
        import sys
        return False
    except:
        return True


if __name__ == '__main__':
    clear(800, 600) if is_codeboot() else speed(0)

    dna_strands = [
        dna,
        get_reverse_complement(dna)
    ]

    genes_coordinates = list(map(get_genes_coordinate_from_dna, dna_strands))
    genes = get_genes_by_coordinate(genes_coordinates)

    amino_acids_chains_fullname = []
    amino_acids_chains_shorten = []

    for gene in genes:
        arn = get_arn_sequence(gene)
        amino_acids_chains_fullname.append(get_amino_acids_chain(arn))
        amino_acids_chains_shorten.append(get_abbreviated_amino_acids_chain(arn))

    # print(amino_acids_chains_fullname)
    # print(amino_acids_chains_shorten)


    for amino_acids_chain in amino_acids_chains_fullname:
        xyz = amino_acids_chain
        x = get_stringify_protein(amino_acids_chain)
        print(x)

    draw_proteins(amino_acids_chains_shorten)

    mainloop() if not is_codeboot() else None