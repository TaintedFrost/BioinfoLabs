from collections import Counter
from itertools import product
from pathlib import Path
import matplotlib.pyplot as plt

#FASTA_FILE = "GCF_001503615.1_ViralProj307951_genomic.fna"
FASTA_FILE = "covid.fna"


BASES = "ACGT"
ALL_CODONS = [
    "".join(codon)
    for codon in product(BASES, repeat=3)
]


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


def main():
    if not Path(FASTA_FILE).exists():
        print("Genome file not found:", FASTA_FILE)
        return

    sequences = read_fasta(FASTA_FILE)
    counts = Counter()
    total = 0

    # Sliding window: 3 letters, moving forward by 1.
    for sequence in sequences:
        for i in range(len(sequence) - 2):
            triplet = sequence[i:i + 3]

            # Ignore ambiguous nucleotides such as N.
            if all(base in BASES for base in triplet):
                counts[triplet] += 1
                total += 1

    # Ensure all 64 possible triplets are represented.
    for codon in ALL_CODONS:
        counts.setdefault(codon, 0)

    top_10 = sorted(
        counts.items(),
        key=lambda item: (-item[1], item[0])
    )[:10]

    print("Total valid windows:", total)
    print("\nTop 10 most frequent triplets:")
    print("Triplet    Frequency")

    for codon, frequency in top_10:
        print(f"{codon}       {frequency}")

    # Plot the 10 most frequent triplets.
    labels = [codon for codon, frequency in top_10]
    values = [frequency for codon, frequency in top_10]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, values)

    plt.title("Top 10 Most Frequent Genome Triplets")
    plt.xlabel("DNA triplet")
    plt.ylabel("Frequency")
    plt.bar_label(bars, padding=3)
    plt.tight_layout()

    plt.savefig("top_10_triplet_frequencies.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
