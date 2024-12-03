
from bibliotek import Utfør_funksjoner #, Proteinvisualisering

# Les DNA-sekvens fra fil
inputfile = "Larven.txt"
with open(inputfile, "r") as f:
    seq = f.read()

start = Utfør_funksjoner(seq, aldri_mett = "Ja")

# Visualiserer det valgte proteinet
protein_navn, protein_sekvens = start.analyser_protein()
