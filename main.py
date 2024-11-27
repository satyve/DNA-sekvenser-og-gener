
from bibliotek import dna_to_rna, oversett, find_protein, proteiner, visualiser_protein_sekvens, analyser_protein

# Les DNA-sekvens fra fil
inputfile = "DNA.txt"
with open(inputfile, "r") as f:
    seq = f.read()

protein_navn, protein_sekvens = analyser_protein(seq, dna_to_rna, oversett, find_protein, proteiner)

# Visualiserer den valgte proteinsekvensen
visualiser_protein_sekvens(protein_sekvens, protein_navn, "nei")