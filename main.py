
from bibliotek import Utfør_funksjoner #, Proteinvisualisering

# Les DNA-sekvens fra fil
inputfile = "DNA.txt"
with open(inputfile, "r") as f:
    seq = f.read()

start = Utfør_funksjoner(seq)

# Visualiserer den valgte proteinsekvensen
protein_navn, protein_sekvens = start.analyser_protein()