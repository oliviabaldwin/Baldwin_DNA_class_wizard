#!/usr/bin/env python
# coding: utf-8

standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}

# # Class called seq
#
# ## This class has three methods:
# ### 1. info -this should print the name, type, organism, and sequence of the instance
# ### 2. length -this should count the length of the sequence string
# ### 3. fasta_out -this should write the name, organism, type, and sequence as a fasta file.


class seq:
    def __init__(self, name, organism, sequence, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    def info(self):
        print(self.name)
        print(self.organism)
        print(self.sequence)
        print(self.type)

    def length(self):
        length = len(self.sequence)
        print(length)

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


# # Child class called protein.
# ## Overwrite the parent class seq function fasta_out to include the protein size in the first line of the fasta file


class protein(seq):
    def __init__(self, name, organism, sequence, type, size):
        super().__init__(name, organism, sequence, type)
        self.size = size

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            self.size
            + ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()

    def mol_weight(self):
        mw = sum(aa_mol_weights[x] for x in self.sequence)
        print(mw)


# Child class called nucleotide.
# ## New method called gc_content that calculates the percent of letters that are G or C and then prints the gc content percentage


class nucleotide(seq):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def gc_content(self):
        counter = self.sequence.count("G") + self.sequence.count("C")
        percentage = (counter / len(self.sequence)) * 100
        print(percentage)


# # "Grandchild classes" DNA and RNA, which will be child classes of nucleotide.

# ## DNA Class and Methods:
# ### transribe - changes DNA to RNA
# ### six_frames - prints 3 reading frames for coding strand and 3 reading frames for the reverse complement strand
# ### reverse_complement - prints the reverse complement of the DNA strand


class DNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def transcribe(self):
        RNA = self.sequence.replace("T", "U")
        print(RNA)

    def six_frames(self):
        print(self.sequence)
        print(self.sequence[1:])
        print(self.sequence[2:])
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        reverse_comp = "".join(
            complement.get(base, base) for base in reversed(self.sequence)
        )
        print(reverse_comp)
        print(reverse_comp[1:])
        print(reverse_comp[2:])

    def reverse_complement(self):
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        rc = "".join(complement.get(base, base) for base in reversed(self.sequence))
        print(rc)


# ## RNA class and methods
# ### start - finds the start codon 'AUG'
# ### translate - translates the RNA to amino acid sequence


class RNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def start(self):
        print(self.sequence.find("AUG"))

    def translate(self):
        start = self.sequence.find("AUG")
        coding_seq = self.sequence[start : len(self.sequence)]
        aa_seq = ""
        for i in range(0, len(coding_seq), 3):
            codon = coding_seq[i : i + 3]
            aa_seq += standard_code.get(codon, "-")
        print(aa_seq)


# # Part B - Testing New Functions

# ##DNA Test
uidA = DNA(
    name="uidA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    organism="Bacteria",
    type="DNA",
)

uidA.fasta_out()
uidA.six_frames()
uidA.reverse_complement()
uidA.transcribe()

# ##RNA Test
uidA_RNA = RNA(
    name="uidA_RNA",
    sequence="CGCAUGUUACGUCCUGUAGAAACCCCAACCCGUGAAAUCAAAAAA",
    organism="Bacteria",
    type="RNA",
)

uidA_RNA.fasta_out()
uidA_RNA.translate()

# ##Protein Test
uidA_protein = protein(
    name="aidA_protein",
    sequence="MLRPVETPTREIKK",
    organism="Bacteria",
    type="protein",
    size="100",
)

uidA_protein.fasta_out()
uidA_protein.mol_weight()
