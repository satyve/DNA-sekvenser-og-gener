import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Dna:
    def __init__(self, seq):
        self.seq = seq

    def dna_to_rna(self):
        return self.seq.replace("T", "U")

    def oversett_rna_til_amino(self):
        amino_oversikt = { 
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UAA': ' ', 'UAG': ' ', 'UGA': ' ',  # Stopp-kodoner
            'UGC': 'C', 'UGU': 'C', 'UGG': 'W',
        }

        protein = []
        lagre_aminosyre = ""
        for i in range(0, len(self.rna), 3):
            codon = self.rna[i:i + 3]
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

class Proteinhåndtering:

    def __init__(self, protein_sequence):
        self.protein_sequence = protein_sequence
        
    proteiner = {
        "Hemoglobin Beta": "MVLSLCAKVEITERNL", 
        "Insulin": "MPLALGRAERPSVEVVKVVELV",
        "Myosin Light Chain": "MEAERGGGWKVGACWLLNL", 
        "Collagen Type I Alpha 1": "MPEGQRKERGGSLGCGGSVLN", 
        "Lactase": "MTHSSTLLPRRLMSVSHPCK", 
        "Actin": "MDDIYETEQFVDDGVTPES", 
        "Cytochrome C Oxidase": "MVDCACFCGGFSA", 
        "Amylase": "MAGGAEFLQGLS", 
        "Keratin": "MSPEALQVEAGARAGSDP", 
        "Elastin": "MAPALPACGAPQAGPP", 
        "Chymotrypsin": "MAGQGARTLGCGA", 
        "Immunoglobulin Heavy Chain": "MSGPGKLWVVGEGGLE", 
        "Tubulin": "MTRSVVLPARFSWFD"   
    }

    def find_protein(self):
        de_funnet = []
        for amino_acid in self.protein_sequence:  # Går gjennom aminosyresekvenser fra RNA-oversett_rna_til_aminoelsen
            for name, sequence in self.proteiner.items():  # Sjekker mot hver kjent proteinsekvens
                if amino_acid in sequence:  # Ser etter match mellom aminosyrer og proteinsekvenser
                    print(f"Aminosyrene {amino_acid} tilsvarer: {name}")
                    de_funnet.append(name)  # Legger til proteinet i listen over funn
        return de_funnet


class Proteinvisualisering:

    def __init__(self, protein, name):
        self.protein = protein
        self.name = name

    # Funksjon for å sette farge basert på proteinindeks
    def fargen(self, i, protein):
        if i == len(protein) - 1:
            return 'darkred'
        elif i % 2 == 0:
            return 'darkgreen'
        else:
            return 'olivedrab'

    def visualiser_protein_sekvens(self, protein, navn, aldri_mett):
        vinkel = np.linspace(0, 2 * np.pi * len(protein), len(protein))
        z = np.linspace(0, 10, len(protein))
        r = 1
        x = r * np.cos(vinkel)
        y = r * np.sin(vinkel)

        figur = plt.figure()
        ax = figur.add_subplot(111, projection='3d')

        for i in range(len(protein)):
            if aldri_mett == "ja":
                farge = self.fargen(i, protein)  # Kall fargen-funksjonen
                ax.scatter(x[i], y[i], z[i], color=farge, s=800)
                ax.text(x[i], y[i], z[i], protein[i], size=8, color=farge)
            else:
                ax.scatter(x[i], y[i], z[i], s=25)  # Plotter punkter på spiralen
                ax.text(x[i] + 0.1, y[i], z[i] + 0.1, protein[i], size=10, color='red')


        ax.set_xlabel('X-akse')
        ax.set_ylabel('Y-akse')
        ax.set_zlabel('Z-akse')
        ax.set_title(f'3D Spiral av Proteinsekvens: {navn}')

        plt.show()
        

class Utfør_funksjoner(Dna, Proteinhåndtering, Proteinvisualisering):
    
    def __init__(self, seq, dna_to_rna, oversett_rna_til_amino, find_protein, proteiner):
        self.seq = seq
        #self.dna_to_rna = dna_to_rna
        #self.oversett_rna_til_amino = oversett_rna_til_amino
        #self.find_protein = find_protein
        self.proteiner = proteiner
        super().__init__(dna_to_rna, oversett_rna_til_amino, find_protein)
        self.dna_to_rna = dna_to_rna
        self.oversett_rna_til_amino = oversett_rna_til_amino
        self.find_protein = find_protein

    def analyser_protein(self):
        # Rydder. Fjerner linjeskift 
        seq = seq.replace("\n", "").replace("\r", "")
        
        # utfører oversett_rna_til_aminoelsen
        rna = self.dna_to_rna(seq)
        protein_sequence = self.oversett_rna_til_amino(rna)
        
        # finner proteinene
        de_funnet = self.find_protein(protein_sequence)
        print("Funnet proteiner:", de_funnet)
        # sjekker om noen er funnet
        if not de_funnet:
            print("Ingen proteiner ble funnet som samsvarer med sekvensen.")
        else:
            # Lar folk (Rasmus og oss) velge hvilket protein vi ønkser å se
            print("Velg et protein fra listen:")
            for i, navn in enumerate(de_funnet):
                print(f"{i + 1}. {navn}")
            
            valg = int(input("\nSkriv nummeret til proteinet du vil visualisere: ")) - 1
            protein_navn = de_funnet[valg]
            protein_sekvens =self.proteiner[protein_navn]
            
            return protein_navn, protein_sekvens