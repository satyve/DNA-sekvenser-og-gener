# Oppskrift 
# - finne dna datasett - X
# - lage program som leser datasett - X
# - lage prpogroam som krypterer og lager komplimentert RNA-molekyl som tekst eller liste - X
# - oversikt over aminosyrer i riktig rekkefølge 
# - stoppe programmet om man finner stopp kodon. 
inputfile = "Chimpanzee.txt"
with open(inputfile, "r") as f:
    seq = f.read()

# Fjerner linjeskift
seq = seq.replace("\n", "").replace("\r", "")


def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U") #.replace("A", "U").replace("C", "G").replace("G", "C")


def translate(rna): 
    aminosyrer = { 
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T', 
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',                 
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P', 
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R', 
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G', 
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L', 
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W', 
    } 
    protein = ""
    for i in range(0, len(seq), 3):
        codon = rna[i:i + 3]
        if len(codon) < 3:
            print(f"Incomplete codon skipped: {codon}")
            continue
        if codon in aminosyrer:
            amino_acid = aminosyrer[codon]
            if amino_acid == '_':  # Stopp kodon funnet
                protein += ' '
            protein += amino_acid
        else:
            print(f"Unknown codon encountered: {codon}")
            protein += '?'  # Håndter uventede kodoner
    return protein


# Utfør oversettelsen og skriv ut resultatet
rna = dna_to_rna(seq)
protein_sequence = translate(rna)
print(protein_sequence)
