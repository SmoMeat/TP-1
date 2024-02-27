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

from bioInfo import adn

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

def get_genes(dna, genes_coordinates):
    genes = []
    for gene_coordinate in genes_coordinates:
        start = gene_coordinate[0]
        end = gene_coordinate[1]
        genes.append(dna[start:end+3])
    return genes

if __name__ == '__main__':
    first_dna = adn
    second_dna = get_reverse_complement(first_dna)

    first_dna_genes_coordinates = get_genes_coordinates(get_position_of_starting_codon(first_dna),
                                                      get_position_of_ending_codon(first_dna))
    second_dna_genes_coordinates = get_genes_coordinates(get_position_of_starting_codon(second_dna),
                                                       get_position_of_ending_codon(second_dna))

    first_dna_genes = get_genes(first_dna, first_dna_genes_coordinates)
    second_dna_genes = get_genes(second_dna, second_dna_genes_coordinates)
    print(get_genes(first_dna, first_dna_genes_coordinates))

