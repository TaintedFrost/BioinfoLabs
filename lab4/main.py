import os
import random
from collections import defaultdict

# ============================================================
# PARAMETERS
# ============================================================

FASTA_FILE = "covid.fna"

N_SAMPLES = 20_000
MIN_LEN = 100
MAX_LEN = 150
MIN_OVERLAP = 10

RANDOM_SEED = 42

# ============================================================
# 1. READ ORIGINAL GENOME FROM FASTA
# ============================================================

folder = os.path.dirname(os.path.abspath(__file__))
fasta_path = os.path.join(folder, FASTA_FILE)

if not os.path.exists(fasta_path):
    raise FileNotFoundError(
        f"Could not find {FASTA_FILE}. "
        "Place it in the same folder as main.py."
    )

with open(fasta_path, "r") as source:
    parts = []

    for line in source:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if parts:
                break
            continue

        parts.append(line)

genome = "".join(parts).upper()

# Keep only valid DNA bases
genome = "".join(base for base in genome if base in "ATCG")

if len(genome) < MAX_LEN:
    raise ValueError("The genome is shorter than the maximum sample length.")

genome_path = os.path.join(folder, "genome.txt")

with open(genome_path, "w") as target:
    target.write(genome)

print("=" * 60)
print("GENOME SAMPLING")
print("=" * 60)

print(f"Original genome length: {len(genome):,} nucleotides")

# ============================================================
# 2. GENERATE 20,000 RANDOM SAMPLES
# ============================================================

random.seed(RANDOM_SEED)

samples = []

for _ in range(N_SAMPLES):

    length = random.randint(MIN_LEN, MAX_LEN)

    start = random.randint(0, len(genome) - length)

    sample = genome[start:start + length]

    samples.append(sample)

samples_path = os.path.join(folder, "samples.txt")

# Save each sample on a separate line
with open(samples_path, "w") as target:
    target.write("\n".join(samples))

sample_coverage = sum(len(s) for s in samples) / len(genome)

print(f"Samples generated: {len(samples):,}")
print(f"Sample lengths: {MIN_LEN}-{MAX_LEN} nucleotides")
print(f"Sample coverage: {sample_coverage:.2f}x")
print(f"Samples saved to: samples.txt")

# ============================================================
# 3. LOAD THE SAMPLES
# ============================================================

print("\n" + "=" * 60)
print("GENOME RECONSTRUCTION")
print("=" * 60)

with open(samples_path, "r") as source:
    samples = [
        line.strip()
        for line in source
        if line.strip()
    ]

max_sample_length = max(len(s) for s in samples)

print(f"Loaded samples: {len(samples):,}")

# ============================================================
# 4. BUILD INDEX
# ============================================================

def build_index(sequences):
    """
    Create an index using the first MIN_OVERLAP
    nucleotides of every sequence.
    """

    index = defaultdict(list)

    for i, sequence in enumerate(sequences):
        prefix = sequence[:MIN_OVERLAP]
        index[prefix].append(i)

    return index


# ============================================================
# 5. EXTEND SEED TO THE RIGHT
# ============================================================

def extend_right(seed, sequences, index, used):
    """
    Extend the seed by finding samples whose beginnings
    overlap with the end of the current seed.

    Only overlaps of at least MIN_OVERLAP are accepted.
    """

    while True:

        best_index = None
        best_overlap = 0

        n = len(seed)

        # Search possible overlap positions near the end
        # of the current seed.
        for pos in range(
            max(0, n - max_sample_length),
            n - MIN_OVERLAP + 1
        ):

            # Get the sequence at this position
            key = seed[pos:pos + MIN_OVERLAP]

            # Only inspect samples with matching prefixes
            for i in index.get(key, []):

                if used[i]:
                    continue

                sample = sequences[i]

                overlap = n - pos

                # Check if the sample is already contained
                # inside the current seed.
                if len(sample) <= overlap:

                    if seed[pos:pos + len(sample)] == sample:
                        used[i] = True

                    continue

                # Check the complete overlap
                if sample[:overlap] == seed[pos:]:

                    # Prefer the longest overlap
                    if overlap > best_overlap:

                        best_index = i
                        best_overlap = overlap

        # No matching sample found
        if best_index is None:
            return seed

        # Mark the selected sample as used
        used[best_index] = True

        # Add only the nucleotides that are not
        # already present in the seed
        seed += sequences[best_index][best_overlap:]

# ============================================================
# 6. ASSEMBLE THE GENOME
# ============================================================

# Reverse all samples to allow assembly from the left
reversed_samples = [s[::-1] for s in samples]

# Create indexes for both directions
index_forward = build_index(samples)
index_reverse = build_index(reversed_samples)

# Track which samples have already been used
used = [False] * len(samples)

assembled_strings = []

# Start with the longest unused sample
sorted_indices = sorted(
    range(len(samples)),
    key=lambda i: -len(samples[i])
)

for i in sorted_indices:

    if used[i]:
        continue

    used[i] = True

    seed = samples[i]

    # Extend to the right
    seed = extend_right(
        seed,
        samples,
        index_forward,
        used
    )

    # Reverse the seed and extend to the right again.
    # This effectively extends the original seed to the left.
    seed = extend_right(
        seed[::-1],
        reversed_samples,
        index_reverse,
        used
    )[::-1]

    assembled_strings.append(seed)

    # Show progress
    if len(assembled_strings) % 100 == 0:
        print(
            f"Contigs assembled: {len(assembled_strings):,}"
        )

# Sort assembled sequences from longest to shortest
assembled_strings.sort(key=len, reverse=True)

# Save assembled sequences
strings_path = os.path.join(folder, "strings.txt")

with open(strings_path, "w") as target:
    target.write("\n".join(assembled_strings))

print("\nAssembly finished!")

print(f"Number of assembled strings: {len(assembled_strings):,}")
print(f"Longest assembled string: {len(assembled_strings[0]):,} nt")
print(
    f"Total assembled length: "
    f"{sum(map(len, assembled_strings)):,} nt"
)

print("Assembled sequences saved to: strings.txt")

# ============================================================
# 7. CHECK GENOME INTEGRITY
# ============================================================

print("\n" + "=" * 60)
print("GENOME INTEGRITY CHECK")
print("=" * 60)

# Track which original genome positions are covered
covered = bytearray(len(genome))

wrong_joins = []

for contig in assembled_strings:

    # Find every exact occurrence of this contig
    pos = genome.find(contig)

    if pos == -1:

        # The assembled sequence does not exist
        # in the original genome
        wrong_joins.append(contig)

    else:

        # Mark all occurrences in the original genome
        while pos != -1:

            covered[pos:pos + len(contig)] = (
                b"\x01" * len(contig)
            )

            pos = genome.find(contig, pos + 1)

# Calculate statistics
covered_positions = sum(covered)

reconstruction_percentage = (
    100 * covered_positions / len(genome)
)

correct_contigs = len(assembled_strings) - len(wrong_joins)

print(f"Original genome length:       {len(genome):,} nt")
print(f"Sample coverage:              {sample_coverage:.2f}x")
print(f"Number of samples:            {len(samples):,}")
print(f"Assembled strings:            {len(assembled_strings):,}")
print(f"Correct strings:              {correct_contigs:,}")
print(f"Strings with wrong joins:     {len(wrong_joins):,}")
print(f"Genome reconstructed:         {reconstruction_percentage:.2f}%")

# ============================================================
# 8. FINAL RESULT
# ============================================================

if len(assembled_strings) == 1 and assembled_strings[0] == genome:

    print("\nPERFECT: The original genome was reconstructed exactly!")

elif reconstruction_percentage >= 99:

    print("\nExcellent reconstruction!")

else:

    print("\nThe genome was partially reconstructed.")

print("=" * 60)
