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

# Utfør oversettelsen
rna = dna_to_rna(seq)
protein_sequence = oversett(rna)
# print(rna)
# print(protein_sequence)

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
  
# Visualiserer det valgte proteinet
visualiser_protein_sekvens(protein_sekvens, protein_navn, "ja")