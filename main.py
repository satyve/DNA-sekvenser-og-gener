
from bibliotek import Utfør_funksjoner, Proteinvisualisering

# Les DNA-sekvens fra fil
inputfile = "DNA.txt"
with open(inputfile, "r") as f:
    seq = f.read()

start = Utfør_funksjoner(seq)

protein_navn, protein_sekvens = start.analyser_protein()

#visualiser = Proteinvisualisering(protein_navn, protein_sekvens)

#visualiser.visualiser_protein_sekvens(protein_navn, protein_sekvens, "nei")

#protein_navn, protein_sekvens = analyser_protein(seq, dna_to_rna, oversett, find_protein, proteiner)

# Visualiserer den valgte proteinsekvensen
#visualiser_protein_sekvens(protein_sekvens, protein_navn, "nei")