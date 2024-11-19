import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Les DNA-sekvens fra fil
inputfile = "Larven.txt"
with open(inputfile, "r") as f:
    seq = f.read()

# Fjerner linjeskift
seq = seq.replace("\n", "").replace("\r", "")

# Funksjon for å konvertere DNA-sekvens til RNA
def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U")

# Funksjon for å oversette RNA-sekvens til en aminosyresekvens (protein)
def translate(rna): 
    amino_oversikt = { 
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T', 
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P', 
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R', 
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G', 
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L', 
        'UAC':'Y', 'UAU':'Y', 'UAA':' ', 'UAG':' ', 
        'UGC':'C', 'UGU':'C', 'UGA':' ', 'UGG':'W', 
    } 

    protein = []
    lagre_aminosyre = ""
    
    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]
        
        if len(codon) < 3:
            continue
        
        if codon in amino_oversikt:
            amino_acid = amino_oversikt[codon]
            if amino_acid == ' ':
                protein.append(lagre_aminosyre)
                lagre_aminosyre = ""
            else:
                lagre_aminosyre += amino_acid
        else:
            protein.append('?')
    
    return protein

# Liste over kjente proteiner og deres aminosyresekvenser
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

# # Velg et protein fra listen
# print("Velg et protein fra listen:")
# for i, navn in enumerate(proteiner.keys()):
#     print(f"{i + 1}. {navn}")

# Utfør oversettelsen
rna = dna_to_rna(seq)
protein_sequence = translate(rna)
# print(rna)
# print(protein_sequence)

de_funnet = []
def find_protein():
    for amino_acid in protein_sequence:
        for name, sequence in proteiner.items():
            # Check if the amino acid is in the current protein sequence from the dictionary
            if amino_acid in sequence:
                print(f"Aminosyrene {amino_acid} tilsvarer: {name}")
                #de_funnet.append(name)

find_protein()

# La brukeren velge et protein fra listen
print("Velg et protein fra listen:")
for i, navn in enumerate(de_funnet):
    print(f"{i + 1}. {navn}")

valg = int(input("\nSkriv nummeret til proteinet du vil visualisere: ")) - 1  # Brukeren velger et nummer
protein_navn = de_funnet[valg] #list(proteiner.keys())[valg]  # Finner navnet på proteinet basert på valget
protein_sekvens = proteiner[protein_navn]

# Funksjon for å visualisere proteinsekvensen
def visualize_protein_sequence(protein, navn):
    theta = np.linspace(0, 2 * np.pi * len(protein), len(protein))
    z = np.linspace(0, 10, len(protein))
    r = 1
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    figur = plt.figure()
    ax = figur.add_subplot(111, projection='3d')

    for i in range(len(protein)):
        if i == len(protein) - 1:
            farge = 'darkred'
        elif i % 2 == 0:
            farge = 'darkgreen'
        else:
            farge = 'olivedrab'

        #ax.scatter(1.15, 0, 10, color='yellow', s=60)  # Første prikk
        #ax.scatter(0.85, 0, 10, color='blue', s=60)  # Andre prikk på samme sted

        ax.scatter(x[i], y[i], z[i], color=farge, s=800)
        ax.text(x[i], y[i], z[i], protein[i], size=8, color=farge)

    ax.set_xlabel('X-akse')
    ax.set_ylabel('Y-akse')
    ax.set_zlabel('Z-akse')
    ax.set_title(f'3D Spiral av Proteinsekvens: {navn}')

    plt.show()

# Visualiserer det valgte proteinet
visualize_protein_sequence(protein_sekvens, protein_navn)
