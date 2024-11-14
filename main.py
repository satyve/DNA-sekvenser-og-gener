# Oppskrift 
# - finne dna datasett - X
# - lage program som leser datasett - X
# - lage prpogroam som krypterer og lager komplimentert RNA-molekyl som tekst eller liste - X
# - oversikt over amino_oversikt i riktig rekkefølge 
# - stoppe programmet om man finner stopp kodon. X

inputfile = "DNA.txt"
with open(inputfile, "r") as f:
    seq = f.read()

# Fjerner linjeskift
seq = seq.replace("\n", "").replace("\r", "")


def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U") #.replace("A", "U").replace("C", "G").replace("G", "C")



def translate(rna): 
    amino_oversikt = { 
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
        'UAC':'Y', 'UAU':'Y', 'UAA':' ', 'UAG':' ', 
        'UGC':'C', 'UGU':'C', 'UGA':' ', 'UGG':'W', 
    } 
    protein = []
    lagre_aminosyre = ""
    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]
        if len(codon) < 3:
            continue
        if codon in amino_oversikt:
            amino_acid = amino_oversikt[codon]
            if amino_acid == ' ':  # Stopp kodon funnet
                protein.append(lagre_aminosyre)
                lagre_aminosyre = ""  # Start en ny sekvens
            else:
                lagre_aminosyre += amino_acid 
        else:
            protein.append('?')  # Marker uventede kodoner
    return protein


"""
Tidligere verson
if lagre_aminosyre:  # Legg bare til hvis proteinet ikke er tomt
                    protein.append(lagre_aminosyre)
                lagre_aminosyre = ""  # Start en ny sekvens
"""



# Utfør oversettelsen og skriv ut resultatet
rna = dna_to_rna(seq)
#print(rna)
protein_sequence = translate(rna)
print(protein_sequence)


def find_protein():
    proteiner = {
        "Hemoglobin Beta" : "MVLSLCAKVEITERNL", 
        "Insulin" : "MPLALGRAERPSVEVVKVVELV",
        "Myosin Light Chain" : "MEAERGGGWKVGACWLLNL", 
        "Collagen Type I Alpha 1" : "MPEGQRKERGGSLGCGGSVLN", 
        "Lactase" : "MTHSSTLLPRRLMSVSHPCK", 
        "Actin" : "MDDIYETEQFVDDGVTPES", 
        "Cytochrome C Oxidase" : "MVDCACFCGGFSA", 
        "Amylase" : "MAGGAEFLQGLS", 
        "Keratin" : "MSPEALQVEAGARAGSDP", 
        "Elastin" : "MAPALPACGAPQAGPP", 
        "Chymotrypsin" : "MAGQGARTLGCGA", 
        "Immunoglobulin Heavy Chain" : "MSGPGKLWVVGEGGLE", 
        "Tubulin" : "MTRSVVLPARFSWFD"   
    }

    #translated_protein = ''.join(translate(rna))

    translated_protein = translate(rna)

    # Iterate through each amino acid in the translated protein
    for amino_acid in translated_protein:
        for name, sequence in proteiner.items():
            # Check if the amino acid is in the current protein sequence from the dictionary
            if amino_acid in sequence:
                print(f"Fant aminosyrene {amino_acid} som tilsvarer: {name}")

find_protein()