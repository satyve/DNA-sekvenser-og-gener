
from bibliotek import Utfør_funksjoner

# Les DNA-sekvens fra fil
inputfile = "Larven.txt"
with open(inputfile, "r") as f:
    seq = f.read()

protein_navn, protein_sekvens = analyser_protein(seq, dna_to_rna, oversett, find_protein, proteiner)

# Visualiserer det valgte proteinet
visualiser_protein_sekvens(protein_sekvens, protein_navn, "ja")