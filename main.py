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
from turtle import *

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
            raise ValueError

    return reverse_complement
    
def get_position_of_starting_codon(dna):
    starting_codon_position = []
    for i in range(len(dna)-2):
        codon = dna[i] + dna[i+1] + dna[i+2]
        if codon == 'TAC':
            starting_codon_position.append(i)
    return starting_codon_position

def get_position_of_ending_codon(dna):
    starting_codon_position = []
    for i in range(len(dna)-2):
        codon = dna[i] + dna[i+1] + dna[i+2]
        if codon in ['ATT', 'ATC', 'ACT']:
            starting_codon_position.append(i)
    return starting_codon_position

def get_genes_coordinates(starting_codons, ending_codons):
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

def get_genes(abc):
    dna = abc['dna']
    genes_coordinates = abc['genes_position']
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
    # x = y = 0
    # side_length = 15
    _x, _y = x, y
    for i, amino_acid in enumerate(amino_acids_chain):
        if i % 15 == 0 and i != 0:
            _y -= side_length
            _x = x
        draw_square(_x, _y, amino_acid, side_length)
        _x += side_length

def draw_proteins(amino_acids_chains):
    x = y = 0
    for amino_acids_chain in amino_acids_chains:
        draw_protein(x, y, amino_acids_chain)
        y = pos()[1] - 20

def draw_square(x, y, amino_acid, side_length=15):
    """For visual studio"""
    penup(), goto(x, y), pendown()
    write('  '+ amino_acid, move=False, align='left', font=('Arial', 8, 'normal'))
    for side in range(4):
        fd(side_length), lt(90)

def get_genes_coordinate_from_dna(dna):
    return {
        "dna": dna,
        "genes_position": get_genes_coordinates(
            get_position_of_starting_codon(dna),
            get_position_of_ending_codon(dna)
        )
    }

# def draw_square(x, y, amino_acid, side_length=15):
#     """For code boot"""
#     penup(), goto(x+side_length/2, y+side_length/2), write(amino_acid)
#     goto(x,y), pendown()
#     for side in range(4):
#         fd(side_length), lt(90)


if __name__ == '__main__':
    speed(0)
    # first_dna = dna
    # second_dna = get_reverse_complement(first_dna)

    dna_strands = [
        dna,
        get_reverse_complement(dna)
    ]
    genes_coordinates = list(map(get_genes_coordinate_from_dna, dna_strands))
    genes = list(map(get_genes, genes_coordinates))

    # TODO: mettre l'adn et les genes_coordinates dans list[tuple] et faire un map(get_genes, that_list)

    # first_dna_genes_coordinates = get_genes_coordinates(get_position_of_starting_codon(first_dna),
    #                                                     get_position_of_ending_codon(first_dna))
    # second_dna_genes_coordinates = get_genes_coordinates(get_position_of_starting_codon(second_dna),
    #                                                      get_position_of_ending_codon(second_dna))

    # first_dna_genes = get_genes(first_dna, first_dna_genes_coordinates)
    # second_dna_genes = get_genes(second_dna, second_dna_genes_coordinates)
    
    # genes = first_dna_genes + second_dna_genes

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

    # draw_proteins(amino_acids_chains_shorten)

    # mainloop()