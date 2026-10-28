from pathlib import Path
from collections import defaultdict

# Configuration
FASTA_FILE = "covid.fna"
MIN_SEQUENCE_LENGTH = 1000
MAX_SEQUENCE_LENGTH = 3000

MIN_REPEAT_LENGTH = 6
MAX_REPEAT_LENGTH = 10

MIN_OCCURRENCES = 3
MAX_OCCURRENCES = 5


def read_fasta(filename):
    """Read DNA from a FASTA file."""
    sequence = []

    with open(filename, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if line and not line.startswith(">"):
                sequence.append(line.upper())

    dna = "".join(sequence)

    # Keep only valid DNA bases
    dna = "".join(base for base in dna if base in "ACGT")

    return dna


def find_repetitions(dna):
    """Find DNA sequences of 6-10 bases that occur 3-5 times."""
    repetitions = defaultdict(list)

    for length in range(MIN_REPEAT_LENGTH, MAX_REPEAT_LENGTH + 1):
        for i in range(len(dna) - length + 1):
            fragment = dna[i:i + length]
            repetitions[fragment].append(i)

    results = []

    for fragment, positions in repetitions.items():
        occurrences = len(positions)

        if MIN_OCCURRENCES <= occurrences <= MAX_OCCURRENCES:
            results.append({
                "sequence": fragment,
                "length": len(fragment),
                "occurrences": occurrences,
                "positions": positions
            })

    # Show the most frequent qualifying repetitions first
    results.sort(key=lambda item: (
        -item["occurrences"],
        item["length"],
        item["sequence"]
    ))

    return results


def main():
    filename = Path(__file__).parent / FASTA_FILE

    if not filename.exists():
        print(f"Error: {FASTA_FILE} was not found.")
        print("Place it in the same folder as main.py.")
        return

    dna = read_fasta(filename)

    # Use a sequence between 1000 and 3000 bases
    if len(dna) < MIN_SEQUENCE_LENGTH:
        print("Error: The DNA sequence contains fewer than 1000 valid bases.")
        return

    dna = dna[:MAX_SEQUENCE_LENGTH]

    print("DNA REPETITION DETECTOR")
    print("=" * 60)
    print(f"Source file: {FASTA_FILE}")
    print(f"DNA sequence length: {len(dna)} bases")
    print(
        f"Repeated sequence length: {MIN_REPEAT_LENGTH}-"
        f"{MAX_REPEAT_LENGTH} bases"
    )
    print(
        f"Required occurrences: {MIN_OCCURRENCES}-"
        f"{MAX_OCCURRENCES}"
    )

    results = find_repetitions(dna)

    print(f"\nNumber of qualifying repeated sequences: {len(results)}")
    print("-" * 80)

    if not results:
        print("No repetitions matching the specified conditions were found.")
        return

    print(
        f"{'DNA sequence':<15}"
        f"{'Length':<10}"
        f"{'Occurrences':<15}"
        f"Positions (0-based)"
    )
    print("-" * 80)

    for item in results:
        positions = ", ".join(map(str, item["positions"]))

        print(
            f"{item['sequence']:<15}"
            f"{item['length']:<10}"
            f"{item['occurrences']:<15}"
            f"{positions}"
        )


if __name__ == "__main__":
    main()
