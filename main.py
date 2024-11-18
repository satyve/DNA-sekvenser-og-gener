import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Les DNA-sekvens fra fil
inputfile = "DNA.txt"
with open(inputfile, "r") as f:
    seq = f.read()

# Fjerner linjeskift
seq = seq.replace("\n", "").replace("\r", "")

# Funksjon for å konvertere DNA til RNA
def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U")

# Funksjon for å oversette RNA til aminosyresekvens
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'', 'UAG':'', 'UGA':'',  # Stopp-kodoner
        'UGC':'C', 'UGU':'C', 'UGG':'W',
    }
    
    protein = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i + 3]
        if len(codon) < 3:
            continue
        amino_acid = amino_oversikt.get(codon, '?')
        if amino_acid == '':  # Stopp-kodon funnet
            break
        protein.append(amino_acid)
    return ''.join(protein)

# Global variabel for proteiner
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

# Utfør oversettelsen
rna = dna_to_rna(seq)
protein_sequence = translate(rna)
print(f"Proteinseksvens: {protein_sequence}")

# Funksjon for å finne proteiner i sekvensen
def find_protein(translated_sequence):
    funnet_proteiner = []
    for navn, sekvens in proteiner.items():
        if sekvens in translated_sequence:
            funnet_proteiner.append(navn)

    if funnet_proteiner:
        for protein in funnet_proteiner:
            print(f"Fant proteinet: {protein}")
    else:
        print("Ingen kjente proteiner funnet.")

# Kjør funksjonen for å finne proteiner
find_protein(protein_sequence)

# La brukeren velge et protein fra listen
print("Velg et protein fra listen:")
for i, navn in enumerate(proteiner.keys()):
    print(f"{i + 1}. {navn}")

valg = int(input("\nSkriv nummeret til proteinet du vil visualisere: ")) - 1  # Brukeren velger et nummer
protein_navn = list(proteiner.keys())[valg]  # Finner navnet på proteinet basert på valget
protein_sekvens = proteiner[protein_navn]  # Henter proteinsekvensen

# Funksjon for å visualisere proteinsekvensen som en 3D spiral med etiketter
def visualize_protein_sequence(protein, navn):
    theta = np.linspace(0, 2 * np.pi * len(protein), len(protein))
    z = np.linspace(0, 10, len(protein))
    r = 1
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    
    figur = plt.figure()
    ax = figur.add_subplot(111, projection='3d')
    
    for i in range(len(protein)):
        ax.scatter(x[i], y[i], z[i], s=25)
        ax.text(x[i] + 0.1, y[i], z[i] + 0.1, protein[i], size=10, color='red')
    
    ax.set_xlabel('X-akse')
    ax.set_ylabel('Y-akse')
    ax.set_zlabel('Z-akse')
    ax.set_title(f'3D Spiral av Proteinsekvens: {navn}')
    
    plt.show()

# Visualiserer det valgte proteinet
visualize_protein_sequence(protein_sekvens, protein_navn)
