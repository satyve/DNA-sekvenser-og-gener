from bibliotek import dna_to_rna, oversett, find_protein, proteiner, visualiser_protein_sekvens

# Les DNA-sekvens fra fil
inputfile = "DNA.txt"
with open(inputfile, "r") as f:
    seq = f.read()

# Fjerner linjeskift
seq = seq.replace("\n", "").replace("\r", "")

# Konverter DNA til RNA og oversett
rna_sequence = dna_to_rna(seq)
protein_sequence = oversett(rna_sequence)

# Finn proteiner
de_funnet = find_protein(protein_sequence)
print("Funnet proteiner:", de_funnet)
find_protein(protein_sequence)

# Lar brukeren velge et protein fra listen over funn
print("Velg et protein fra listen:")
for i, navn in enumerate(de_funnet):  # Viser tilgjengelige valg med nummerering
    print(f"{i + 1}. {navn}")

valg = int(input("\nSkriv nummeret til proteinet du vil visualisere: ")) - 1  # Brukerens valg
protein_navn = de_funnet[valg]  # Henter navnet på proteinet basert på brukerens valg
protein_sekvens = proteiner[protein_navn]  # Finner proteinsekvensen for det valgte proteinet

# Funksjon for å visualisere proteinsekvensen som en 3D-spiral med etiketter

# Visualiserer den valgte proteinsekvensen
visualiser_protein_sekvens(protein_sekvens, protein_navn, "nei")
