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

def visualiser_protein_sekvens(protein, navn):
    # Genererer en spiral ved hjelp av trigonometri og linære avstander
    vinkel = np.linspace(0, 2 * np.pi * len(protein), len(protein))  # Vinklene for spiral
    z = np.linspace(0, 10, len(protein))  # Z-koordinater for høyden
    r = 1  # radius 
    x = r * np.cos(vinkel) 
    y = r * np.sin(vinkel)  
    
    figur = plt.figure()  # Oppretter en ny figur for 3D-plot
    ax = figur.add_subplot(111, projection='3d')  # Legger til en 3D-akse, funksjon hentes fra mpl_toolkits.mplot3d
    # 111 = én rad, én kolonne, og er det første (og eneste) plottet i denne matrisen
    
    for i in range(len(protein)):  # Itererer gjennom aminosyrer i sekvensen
        ax.scatter(x[i], y[i], z[i], s=25)  # Plotter punkter på spiralen
        ax.text(x[i] + 0.1, y[i], z[i] + 0.1, protein[i], size=10, color='red')  # Legger til etiketter
    
    # Setter aksetitler og diagramtittel
    ax.set_xlabel('X-akse')
    ax.set_ylabel('Y-akse')
    ax.set_zlabel('Z-akse')
    ax.set_title(f'3D Spiral av Proteinsekvens: {navn}')
    
    plt.show()  # Viser diagrammet