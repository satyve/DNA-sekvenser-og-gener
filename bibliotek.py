import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U")

def oversett(rna):
    amino_oversikt = { 
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': ' ', 'UAG': ' ', 'UGA': ' ',  # Stopp-kodoner
        'UGC': 'C', 'UGU': 'C', 'UGG': 'W',
    }

    protein = []
    lagre_aminosyre = ""
    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]
        if len(codon) < 3:
            continue
        if codon in amino_oversikt:
            amino_acid = amino_oversikt[codon]
            if amino_acid == ' ':  # Stopp kodon funnet
                protein.append(lagre_aminosyre)
                lagre_aminosyre = ""  # Start en ny sekvens
            else:
                lagre_aminosyre += amino_acid
        else:
            protein.append('?')  # Marker uventede kodoner
    return protein


proteiner = {
    "Hemoglobin Beta": "MVLSLCAKVEITERNL", 
    "Insulin": "MPLALGRAERPSVEVVKVVELV",
    "Myosin Light Chain": "MEAERGGGWKVGACWLLNL", 
    "Collagen Type I Alpha 1": "MPEGQRKERGGSLGCGGSVLN", 
    "Lactase": "MTHSSTLLPRRLMSVSHPCK", 
    "Actin": "MDDIYETEQFVDDGVTPES", 
    "Cytochrome C Oxidase": "MVDCACFCGGFSA", 
    "Amylase": "MAGGAEFLQGLS", 
    "Keratin": "MSPEALQVEAGARAGSDP", 
    "Elastin": "MAPALPACGAPQAGPP", 
    "Chymotrypsin": "MAGQGARTLGCGA", 
    "Immunoglobulin Heavy Chain": "MSGPGKLWVVGEGGLE", 
    "Tubulin": "MTRSVVLPARFSWFD"   
}

def find_protein(protein_sequence):
    de_funnet = []
    for amino_acid in protein_sequence:  # Går gjennom aminosyresekvenser fra RNA-oversettelsen
        for name, sequence in proteiner.items():  # Sjekker mot hver kjent proteinsekvens
            if amino_acid in sequence:  # Ser etter match mellom aminosyrer og proteinsekvenser
                print(f"Aminosyrene {amino_acid} tilsvarer: {name}")
                de_funnet.append(name)  # Legger til proteinet i listen over funn
    return de_funnet

# Funksjon for å sette farge basert på proteinindeks
def fargen(i, protein):
    if i == len(protein) - 1:
        return 'darkred'
    elif i % 2 == 0:
        return 'darkgreen'
    else:
        return 'olivedrab'

def visualiser_protein_sekvens(protein, navn, aldri_mett):
    vinkel = np.linspace(0, 2 * np.pi * len(protein), len(protein))
    z = np.linspace(0, 10, len(protein))
    r = 1
    x = r * np.cos(vinkel)
    y = r * np.sin(vinkel)
 
    figur = plt.figure()
    ax = figur.add_subplot(111, projection='3d')
 
    for i in range(len(protein)):
        if aldri_mett == "ja":
            farge = fargen(i, protein)  # Kall fargen-funksjonen
            ax.scatter(x[i], y[i], z[i], color=farge, s=800)
            ax.text(x[i], y[i], z[i], protein[i], size=8, color=farge)
        else:
            ax.scatter(x[i], y[i], z[i], s=25)  # Plotter punkter på spiralen
            ax.text(x[i] + 0.1, y[i], z[i] + 0.1, protein[i], size=10, color='red')

 
    ax.set_xlabel('X-akse')
    ax.set_ylabel('Y-akse')
    ax.set_zlabel('Z-akse')
    ax.set_title(f'3D Spiral av Proteinsekvens: {navn}')
 
    plt.show()

def analyser_protein(seq, dna_to_rna, oversett, find_protein, proteiner):
    # Rydder. Fjerner linjeskift 
    seq = seq.replace("\n", "").replace("\r", "")
    
    # utfører oversettelsen
    rna = dna_to_rna(seq)
    protein_sequence = oversett(rna)
    
    # finner proteinene
    de_funnet = find_protein(protein_sequence)
    print("Funnet proteiner:", de_funnet)

    # sjekker om noen er funnet
    if not de_funnet:
        print("Ingen proteiner ble funnet som samsvarer med sekvensen.")
    else:
        # Lar folk (Rasmus og oss) velge hvilket protein vi ønkser å se
        print("Velg et protein fra listen:")
        for i, navn in enumerate(de_funnet):
            print(f"{i + 1}. {navn}")
        
        valg = int(input("\nSkriv nummeret til proteinet du vil visualisere: ")) - 1
        protein_navn = de_funnet[valg]
        protein_sekvens = proteiner[protein_navn]
        
        return protein_navn, protein_sekvens