from collections import Counter
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
from ollama import chat


FASTA_FILE = "covid.fna"
LLM_NAME = "llama3.2"

BASES = "ACGT"

ALL_TRIPLETS = [
    "".join(triplet)
    for triplet in product(BASES, repeat=3)
]


CODON_GROUPS = {
    "TTT TTC": "Phenylalanine",
    "TTA TTG CTT CTC CTA CTG": "Leucine",
    "ATT ATC ATA": "Isoleucine",
    "ATG": "Methionine",
    "GTT GTC GTA GTG": "Valine",
    "TCT TCC TCA TCG AGT AGC": "Serine",
    "CCT CCC CCA CCG": "Proline",
    "ACT ACC ACA ACG": "Threonine",
    "GCT GCC GCA GCG": "Alanine",
    "TAT TAC": "Tyrosine",
    "TAA TAG TGA": "Stop",
    "CAT CAC": "Histidine",
    "CAA CAG": "Glutamine",
    "AAT AAC": "Asparagine",
    "AAA AAG": "Lysine",
    "GAT GAC": "Aspartic acid",
    "GAA GAG": "Glutamic acid",
    "TGT TGC": "Cysteine",
    "TGG": "Tryptophan",
    "CGT CGC CGA CGG AGA AGG": "Arginine",
    "GGT GGC GGA GGG": "Glycine",
}

CODON_TABLE = {
    triplet: amino_acid
    for group, amino_acid in CODON_GROUPS.items()
    for triplet in group.split()
}


# READ FASTA

def read_fasta(filename):
    sequences = []
    current = []

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
            else:
                current.append(line)

    if current:
        sequences.append("".join(current).upper())

    return sequences


# COUNT TRIPLETS

def count_triplets(sequences):
    counts = Counter()
    total = 0

    # Sliding window moves one nucleotide at a time.
    for sequence in sequences:
        for i in range(len(sequence) - 2):
            triplet = sequence[i:i + 3]

            # Ignore ambiguous bases such as N.
            if all(base in BASES for base in triplet):
                counts[triplet] += 1
                total += 1

    # Include every possible DNA triplet, even if its count is zero.
    for triplet in ALL_TRIPLETS:
        counts.setdefault(triplet, 0)

    return counts, total


# AI ANALYSIS

def ask_ai_about_foods(top_10, amino_acid_counts):
    triplet_information = "\n".join(
        f"{triplet}: {frequency} occurrences -> "
        f"{CODON_TABLE[triplet]}"
        for triplet, frequency in top_10
    )

    amino_acid_information = "\n".join(
        f"{amino_acid}: {frequency} occurrences"
        for amino_acid, frequency in amino_acid_counts.most_common()
    )

    prompt = f"""
I am doing a bioinformatics project analyzing a viral DNA sequence.

TOP 10 MOST FREQUENT OVERLAPPING DNA TRIPLETS:
{triplet_information}

AMINO ACIDS ASSOCIATED WITH THESE TRIPLETS:
{amino_acid_information}

Analyze the nutritional information for educational purposes.

Organize your answer into these sections:

1. AMINO ACIDS IDENTIFIED
List the amino acids represented by the top 10 triplets.
Exclude stop signals because they are not amino acids.

2. FOODS TO LIMIT OR AVOID (FOR COMPARISON)
Give examples of foods that tend to be relatively rich in the
listed amino acids per typical serving.

Explain which amino acids these foods provide and why their
amino-acid content is relatively high.

Do not claim that these foods are harmful or should medically
be avoided because of my genome analysis.

3. FOODS TO CHOOSE
Recommend practical foods that tend to provide relatively low
amounts of several of the listed amino acids per typical serving.

Include specific examples, such as fruits, vegetables, and other
appropriate foods. Explain which listed amino acids they tend
to provide in lower amounts compared with protein-rich foods.

4. BEST FOOD COMPARISONS
Provide a table with these columns:
- Food
- Typical serving
- Relative amounts of the listed amino acids
- Explanation

If exact nutritional data is unavailable, use qualitative terms
rather than inventing numbers.

5. MEAL IDEAS
Give three example meal or snack ideas using foods that are
relatively low in the listed amino acids. Keep the suggestions
balanced and nutritionally adequate.

6. SCIENTIFIC LIMITATIONS
My triplets were counted using a sliding window moving one
nucleotide at a time. They are not necessarily protein-coding
codons, so their frequencies do not establish the amino-acid
composition of viral proteins.

These results do not demonstrate that eating foods low in these
amino acids affects COVID-19 or any other viral infection.

Do not recommend eliminating essential amino acids or following
a restrictive diet. Explain that all amino acids have important
biological roles and that fruits and vegetables can also contain
some amino acids.

Use clear language suitable for a bioinformatics student.
"""

    print("\\nAsking the local AI model...\\n")

    response = chat(
        model=LLM_NAME,
        messages=[
            {
                "role": "user",
                "content": prompt
            }
        ]
    )

    return response["message"]["content"]


# MAIN PROGRAM

def main():
    if not Path(FASTA_FILE).exists():
        print("Genome file not found:", FASTA_FILE)
        return

    sequences = read_fasta(FASTA_FILE)

    if not sequences:
        print("No sequences found in the FASTA file.")
        return

    counts, total = count_triplets(sequences)

    top_10 = sorted(
        counts.items(),
        key=lambda item: (-item[1], item[0])
    )[:10]

    print("Total valid windows:", total)
    print("\nTop 10 most frequent DNA triplets:")
    print(f"{'Triplet':<10}{'Frequency':<12}{'Amino acid'}")

    # Count the occurrences associated with each amino acid
    # among the selected top 10 triplets.
    amino_acid_counts = Counter()

    for triplet, frequency in top_10:
        amino_acid = CODON_TABLE[triplet]

        print(f"{triplet:<10}{frequency:<12}{amino_acid}")

        if amino_acid != "Stop":
            amino_acid_counts[amino_acid] += frequency

    # Plot the top 10 triplets.
    labels = [triplet for triplet, frequency in top_10]
    values = [frequency for triplet, frequency in top_10]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, values)

    plt.title("Top 10 Most Frequent Genome Triplets")
    plt.xlabel("DNA triplet")
    plt.ylabel("Frequency")
    plt.bar_label(bars, padding=3)
    plt.tight_layout()

    plt.savefig("top_10_triplet_frequencies.png", dpi=150)
    plt.show()

    # Ask the AI for an educational food comparison.
    try:
        ai_answer = ask_ai_about_foods(
            top_10,
            amino_acid_counts
        )

        print("\n" + "=" * 60)
        print("AI NUTRITION ANALYSIS")
        print("=" * 60)
        print(ai_answer)

    except Exception as error:
        print("\nCould not contact the local AI model.")
        print("Check that Ollama is running and the model is installed.")
        print("Error:", error)


if __name__ == "__main__":
    main()
