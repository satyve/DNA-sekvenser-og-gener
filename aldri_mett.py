import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from bibliotek import dna_to_rna, oversett, find_protein, proteiner, visualiser_protein_sekvens

 
# Les DNA-sekvens fra fil
inputfile = "Larven.txt"
with open(inputfile, "r") as f:
    seq = f.read()
 
# Fjerner linjeskift
seq = seq.replace("\n", "").replace("\r", "")
 
# Funksjon for å konvertere DNA-sekvens til RNA
def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U")
 

# Utfør oversettelsen
rna = dna_to_rna(seq)
protein_sequence = oversett(rna)
# print(rna)
# print(protein_sequence)

# Utfør oversettelsen
rna = dna_to_rna(seq)
protein_sequence = oversett(rna)
 
 
# Finn proteiner
de_funnet = find_protein(protein_sequence)
print("Funnet proteiner:", de_funnet)
find_protein(protein_sequence)
 
# Sjekk om noe ble funnet
if not de_funnet:
    print("Ingen proteiner ble funnet som samsvarer med sekvensen.")
else:
    # La brukeren velge et protein fra listen
    print("Velg et protein fra listen:")
    for i, navn in enumerate(de_funnet):
        print(f"{i + 1}. {navn}")
 
    valg = int(input("\nSkriv nummeret til proteinet du vil visualisere: ")) - 1
    protein_navn = de_funnet[valg]
    protein_sekvens = proteiner[protein_navn]
    
 
# Funksjon for å visualisere proteinsekvensen
def visualiser_protein_sekvens(protein, navn):
    vinkel = np.linspace(0, 2 * np.pi * len(protein), len(protein))
    z = np.linspace(0, 10, len(protein))
    r = 1
    x = r * np.cos(vinkel)
    y = r * np.sin(vinkel)
 
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
visualiser_protein_sekvens(protein_sekvens, protein_navn)