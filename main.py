# Oppskrift 
# - finne dna datasett
# - lage program som leser datasett
# - lage prpogroam som krypterer og lager komplimentert RNA-molekyl som tekst eller liste
# - oversikt over aminosyrer i riktig rekkefølge 
# - stoppe programmet om man finner stopp kodon. 
inputfile = "Chimpanzee.txt"
with open(inputfile, "r") as f:
    seq = f.read()

# Fjerner linjeskift
seq = seq.replace("\n", "").replace("\r", "")

def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U").replace("A", "U").replace("C", "G").replace("G", "C")
print(dna_seq)

def translate(seq): 
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    for i in range(0, len(seq) - len(seq) % 3, 3): 
        codon = seq[i:i + 3] 
        protein += table.get(codon, "?")  # '?' hvis codon ikke finnes i tabellen
    return protein

# Utfør oversettelsen og skriv ut resultatet
protein_sequence = translate(seq)
print(protein_sequence)
