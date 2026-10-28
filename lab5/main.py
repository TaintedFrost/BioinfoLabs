import random
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# Configuration
FASTA_FILE = "covid.fna"
NUMBER_OF_SAMPLES = 10
MIN_LENGTH = 100
MAX_LENGTH = 3000
RANDOM_SEED = 42


def read_fasta(filename):
    """Read a DNA sequence from a FASTA file."""
    sequence = []

    with open(filename, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            # Ignore FASTA header lines
            if line and not line.startswith(">"):
                sequence.append(line.upper())

    # Keep only valid DNA nucleotides
    dna = "".join(sequence)
    dna = "".join(base for base in dna if base in "ACGT")

    return dna


def create_samples(dna):
    """Generate 10 random DNA fragments and store them in an array."""
    samples = []

    for i in range(NUMBER_OF_SAMPLES):
        length = random.randint(
            MIN_LENGTH,
            min(MAX_LENGTH, len(dna))
        )

        # Choose a random starting position
        start = random.randint(0, len(dna) - length)

        fragment = dna[start:start + length]

        samples.append({
            "id": i + 1,
            "start": start,
            "length": len(fragment),
            "sequence": fragment
        })

    return samples


def migration_distance(length):
    """
    Calculate simulated migration distance.

    Short fragments travel farther.
    Long fragments travel less far.
    """
    min_size = MIN_LENGTH
    max_size = min(MAX_LENGTH, len(dna_sequence))

    if max_size == min_size:
        return 5.0

    max_distance = 10.0

    distance = max_distance * (
        1 - (length - min_size) / (max_size - min_size)
    )

    return distance


def simulate_gel(samples):
    """Draw a simulated agarose gel electrophoresis result."""
    ladder_sizes = [100, 200, 500, 1000, 1500, 2000, 2500, 3000]

    # Keep ladder sizes within the fragment size range
    ladder_sizes = [
        size for size in ladder_sizes
        if size <= min(MAX_LENGTH, len(dna_sequence))
    ]

    if MIN_LENGTH not in ladder_sizes:
        ladder_sizes.insert(0, MIN_LENGTH)

    number_of_lanes = NUMBER_OF_SAMPLES + 1

    fig, ax = plt.subplots(figsize=(12, 7))

    # Gel background
    ax.add_patch(
        Rectangle(
            (0.5, -0.7),
            number_of_lanes,
            12,
            facecolor="#163b35",
            edgecolor="black",
            linewidth=2
        )
    )

    # Draw lane separators
    for lane in range(number_of_lanes):
        ax.add_patch(
            Rectangle(
                (lane + 0.65, -0.5),
                0.7,
                11,
                facecolor="#24574b",
                edgecolor="#356d60",
                linewidth=0.5
            )
        )

    # DNA ladder in the first lane
    ladder_x = 1

    for size in ladder_sizes:
        y = migration_distance(size)

        ax.plot(
            [ladder_x + 0.1, ladder_x + 0.6],
            [y, y],
            color="#ffe66d",
            linewidth=3
        )

        ax.text(
            ladder_x - 0.05,
            y,
            str(size),
            color="black",
            fontsize=8,
            ha="right",
            va="center"
        )

    # Draw the 10 DNA samples
    for i, sample in enumerate(samples):
        x = i + 2
        y = migration_distance(sample["length"])

        ax.plot(
            [x + 0.1, x + 0.6],
            [y, y],
            color="#7cf5ff",
            linewidth=4,
            solid_capstyle="round"
        )

        # Label each band with its fragment size
        ax.text(
            x + 0.35,
            y - 0.25,
            str(sample["length"]),
            color="white",
            fontsize=8,
            ha="center",
            va="bottom"
        )

    # Gel labels and formatting
    ax.set_xlim(-0.3, number_of_lanes + 1)
    ax.set_ylim(11, -1)
    ax.set_xticks(range(1, number_of_lanes + 1))
    ax.set_xticklabels(
        ["Ladder"] + [f"Sample {i}" for i in range(1, 11)],
        rotation=35,
        ha="right"
    )

    ax.set_ylabel("Migration through gel (increasing downward)")
    ax.set_title("Simulated DNA Gel Electrophoresis")

    ax.set_facecolor("#f4f4f4")
    ax.grid(axis="y", linestyle="--", alpha=0.25)

    plt.tight_layout()
    plt.show()


# Main program
if __name__ == "__main__":
    random.seed(RANDOM_SEED)

    file_path = Path(__file__).parent / FASTA_FILE

    if not file_path.exists():
        print(f"Error: {FASTA_FILE} was not found.")
        print("Place it in the same folder as main.py.")
        raise SystemExit(1)

    dna_sequence = read_fasta(file_path)

    if len(dna_sequence) < MIN_LENGTH:
        print("Error: The FASTA file contains fewer than 100 valid DNA bases.")
        raise SystemExit(1)

    # Use a sequence of at most 3000 bases
    dna_sequence = dna_sequence[:MAX_LENGTH]

    print("DNA GEL ELECTROPHORESIS SIMULATION")
    print("-" * 45)
    print(f"Source file: {FASTA_FILE}")
    print(f"DNA sequence length: {len(dna_sequence)} bases")

    # Store the 10 random fragments in an array
    samples = create_samples(dna_sequence)

    print("\nGenerated DNA fragments:")
    print("-" * 65)
    print(f"{'Sample':<10}{'Start':<12}{'Length (bp)':<15}{'Sequence preview'}")
    print("-" * 65)

    for sample in samples:
        preview = sample["sequence"][:20]

        print(
            f"{sample['id']:<10}"
            f"{sample['start']:<12}"
            f"{sample['length']:<15}"
            f"{preview}..."
        )

    print("\nSimulating DNA migration through the gel...")
    print("Shorter fragments migrate farther than longer fragments.")

    # Display the simulated gel
    simulate_gel(samples)
