RNA_SEQUENCE = "AUGGCCAUUUAA"

CODONS = {
    "UUU UUC": "Phe", "UUA UUG CUU CUC CUA CUG": "Leu",
    "UCU UCC UCA UCG AGU AGC": "Ser", "UAU UAC": "Tyr",
    "UAA UAG UGA": "Stop", "UGU UGC": "Cys", "UGG": "Trp",
    "CCU CCC CCA CCG": "Pro", "CAU CAC": "His", "CAA CAG": "Gln",
    "CGU CGC CGA CGG AGA AGG": "Arg",
    "AUU AUC AUA": "Ile", "AUG": "Met",
    "ACU ACC ACA ACG": "Thr", "AAU AAC": "Asn", "AAA AAG": "Lys",
    "GUU GUC GUA GUG": "Val", "GCU GCC GCA GCG": "Ala",
    "GAU GAC": "Asp", "GAA GAG": "Glu",
    "GGU GGC GGA GGG": "Gly",
}

# Convert the grouped codons above into one simple lookup dictionary.
CODON_TABLE = {
    codon: amino_acid
    for codon_group, amino_acid in CODONS.items()
    for codon in codon_group.split()
}

print("RNA sequence:", RNA_SEQUENCE)
print("\nCodon translations:")

protein = []
for i in range(0, len(RNA_SEQUENCE) - 2, 3):
    codon = RNA_SEQUENCE[i:i + 3]
    amino_acid = CODON_TABLE[codon]
    print(f"{codon} -> {amino_acid}")

    if amino_acid == "Stop":
        break
    protein.append(amino_acid)

if len(RNA_SEQUENCE) % 3:
    print("Incomplete codon at the end:", RNA_SEQUENCE[len(RNA_SEQUENCE) // 3 * 3:])

print("\nProtein:", "-".join(protein))
