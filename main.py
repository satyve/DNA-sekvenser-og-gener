from bibliotek import DNAhåndtering

# Les DNA-sekvens fra fil
inputfile = "DNA.txt"
with open(inputfile, "r") as f:
    seq = f.read()

analysering = DNAhåndtering(seq)
protein_navn, protein_sekvens = analysering.analyser_og_hent_protein()

# Visualiserer den valgte proteinsekvensen
analysering.visualiser_protein_sekvens(protein_sekvens, protein_navn, "nei")
