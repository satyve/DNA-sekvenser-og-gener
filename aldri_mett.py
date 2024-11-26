from bibliotek import DNAhåndtering
 
# Les DNA-sekvens fra fil
inputfile = "Larven.txt"
with open(inputfile, "r") as f:
    seq = f.read()

analysering = DNAhåndtering(seq)
protein_navn, protein_sekvens = analysering.analyser_og_hent_protein()

# Visualiserer det valgte proteinet
analysering.visualiser_protein_sekvens(protein_sekvens, protein_navn, "ja")