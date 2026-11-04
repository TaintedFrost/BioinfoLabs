import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict

# ==================================================
# CONFIGURATION
# ==================================================

FASTA_FILE = "covid.fna"

# Minimum length of the Terminal Inverted Repeats
MIN_TIR_LENGTH = 10

# Length of the Target Site Duplications
TSD_LENGTH = 4

# True  -> a candidate needs TIRs AND identical TSDs (strict)
# False -> TIRs alone are enough (many more, mostly chance hits)
REQUIRE_TSD = True

# Size of the DNA between the two TIRs
MIN_INTERNAL_LENGTH = 25
MAX_INTERNAL_LENGTH = 1000

# Ignore trivial TIRs such as AAAAAAAAAA or ATATATATAT
MIN_DISTINCT_BASES = 3

# How many candidates to print per list
MAX_PRINTED = 25

COMPLEMENT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


# ==================================================
# FASTA / DNA FUNCTIONS
# ==================================================

def read_fasta(path):

    records = []

    name = None
    chunks = []

    with open(path) as file:

        for line in file:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:
                    records.append((name, "".join(chunks)))

                name = line[1:].strip()
                chunks = []

            else:

                if name is None:
                    name = "unnamed sequence"

                chunks.append(line.upper().replace("U", "T"))

    if name is not None:
        records.append((name, "".join(chunks)))

    return records


def reverse_complement(sequence):
    return ''.join(
        COMPLEMENT.get(base, 'N')
        for base in reversed(sequence)
    )


def is_valid(sequence):
    # Only A, C, G, T (skips N and other ambiguity codes)
    return set(sequence) <= set("ACGT")


# ==================================================
# K-MER INDEX
#
# Maps every k-mer to the (sorted) list of positions
# where it starts. This lets us find the matching
# right TIR instantly instead of scanning everything.
# ==================================================

def build_kmer_index(dna, k):

    index = defaultdict(list)

    for i in range(len(dna) - k + 1):

        kmer = dna[i:i + k]

        if is_valid(kmer):
            index[kmer].append(i)

    return index


# ==================================================
# TIR MEASUREMENT
#
# Counts how many bases from the outside edges the
# two ends stay reverse-complementary.
# ==================================================

def measure_tir(dna, start, end):

    length = 0

    while start + length < end - length:

        left_base = dna[start + length]
        right_base = dna[end - length]

        if COMPLEMENT.get(left_base) != right_base:
            break

        length += 1

    return length


# ==================================================
# DETECT TRANSPOSABLE ELEMENT CANDIDATES
#
# Conditions:
# 1. TIRs (at least MIN_TIR_LENGTH) are reverse
#    complements of each other.
# 2. Internal DNA is MIN..MAX_INTERNAL_LENGTH bases.
# 3. Flanking TSDs are identical (if REQUIRE_TSD).
#
# Coordinates in the result are 1-based and inclusive.
# ==================================================

def detect_elements(dna):

    k = MIN_TIR_LENGTH

    index = build_kmer_index(dna, k)

    min_length = 2 * k + MIN_INTERNAL_LENGTH
    max_length = 2 * k + MAX_INTERNAL_LENGTH

    first_start = TSD_LENGTH if REQUIRE_TSD else 0

    detected = []

    for start in range(first_start, len(dna) - k + 1):

        left_tir = dna[start:start + k]

        if not is_valid(left_tir):
            continue

        if len(set(left_tir)) < MIN_DISTINCT_BASES:
            continue

        expected_right = reverse_complement(left_tir)

        positions = index.get(expected_right)

        if not positions:
            continue

        # Right TIR starts at p, element ends at p + k - 1,
        # so element length = p + k - start
        low = bisect_left(positions, start + min_length - k)
        high = bisect_right(positions, start + max_length - k)

        for p in positions[low:high]:

            end = p + k - 1

            left_tsd = ""
            right_tsd = ""

            if REQUIRE_TSD:

                if end + 1 + TSD_LENGTH > len(dna):
                    continue

                left_tsd = dna[start - TSD_LENGTH:start]
                right_tsd = dna[end + 1:end + 1 + TSD_LENGTH]

                if left_tsd != right_tsd or not is_valid(left_tsd):
                    continue

            else:

                # Keep only the outermost pair of a longer
                # inverted repeat (avoids near-duplicates)
                if (
                    start > 0
                    and end + 1 < len(dna)
                    and COMPLEMENT.get(dna[start - 1]) == dna[end + 1]
                ):
                    continue

            tir_length = measure_tir(dna, start, end)

            detected.append({
                "start": start + 1,
                "end": end + 1,
                "length": end - start + 1,
                "tir_length": tir_length,
                "internal_length": end - start + 1 - 2 * tir_length,
                "left_tir": dna[start:start + tir_length],
                "right_tir": dna[end + 1 - tir_length:end + 1],
                "left_tsd": left_tsd,
                "right_tsd": right_tsd
            })

    return detected


# ==================================================
# CHECK OVERLAPPING CANDIDATES
# ==================================================

def check_overlaps(elements):

    overlaps = []

    ordered = sorted(elements, key=lambda e: e["start"])

    for i in range(len(ordered)):

        for j in range(i + 1, len(ordered)):

            a = ordered[i]
            b = ordered[j]

            if b["start"] > a["end"]:
                break

            overlap_end = min(a["end"], b["end"])

            overlaps.append({
                "first": (a["start"], a["end"]),
                "second": (b["start"], b["end"]),
                "bases": overlap_end - b["start"] + 1
            })

    return overlaps


# ==================================================
# EXPECTED NUMBER OF CHANCE HITS
#
# Rough estimate for a uniformly random sequence.
# Helps judging whether the detections mean anything.
# ==================================================

def expected_by_chance(dna_length):

    sizes = MAX_INTERNAL_LENGTH - MIN_INTERNAL_LENGTH + 1

    expected = dna_length * sizes * (4 ** -MIN_TIR_LENGTH)

    if REQUIRE_TSD:
        expected *= 4 ** -TSD_LENGTH

    return expected


# ==================================================
# ANALYSE ONE SEQUENCE
# ==================================================

def analyse(name, dna):

    print("\n" + "=" * 70)
    print(name)
    print("=" * 70)

    valid_bases = sum(dna.count(b) for b in "ACGT")

    print("\nSequence length:", len(dna), "bases")
    print("Unknown bases (N etc.):", len(dna) - valid_bases)

    detected = detect_elements(dna)

    print("\nSETTINGS:")
    print("Minimum TIR length:", MIN_TIR_LENGTH)
    print("TSD required:", REQUIRE_TSD,
          f"(length {TSD_LENGTH})" if REQUIRE_TSD else "")
    print("Internal length:",
          f"{MIN_INTERNAL_LENGTH}-{MAX_INTERNAL_LENGTH} bases")

    print("\nDETECTED TRANSPOSABLE ELEMENT CANDIDATES:")

    if not detected:
        print("None found.")

    for i, e in enumerate(detected[:MAX_PRINTED], 1):

        print(
            f"\nCandidate {i}: Start = {e['start']}, "
            f"End = {e['end']}, "
            f"Length = {e['length']} bases"
        )

        print("TIR length:", e["tir_length"])
        print("Left TIR: ", e["left_tir"])
        print("Right TIR:", e["right_tir"])

        if REQUIRE_TSD:
            print("Left TSD: ", e["left_tsd"])
            print("Right TSD:", e["right_tsd"])

    if len(detected) > MAX_PRINTED:
        print(
            f"\n... {len(detected) - MAX_PRINTED} more "
            f"candidates not shown (raise MAX_PRINTED)."
        )

    overlaps = check_overlaps(detected)

    print("\nOVERLAP ANALYSIS:")

    if overlaps:

        for o in overlaps[:MAX_PRINTED]:

            print(
                f"{o['first'][0]}-{o['first'][1]} and "
                f"{o['second'][0]}-{o['second'][1]} overlap by "
                f"{o['bases']} bases."
            )

        if len(overlaps) > MAX_PRINTED:
            print(f"... {len(overlaps) - MAX_PRINTED} more overlaps.")

    else:

        print("No overlapping candidates.")

    print("\nSUMMARY:")
    print("Candidates found:", len(detected))
    print("Overlapping pairs:", len(overlaps))
    print(
        "Expected by pure chance (random DNA): ~"
        f"{expected_by_chance(len(dna)):.3f}"
    )

    return detected


# ==================================================
# MAIN PROGRAM
# ==================================================

def main():

    path = sys.argv[1] if len(sys.argv) > 1 else FASTA_FILE

    print("=" * 70)
    print("TRANSPOSABLE ELEMENT SEARCH IN A REAL GENOME")
    print("TIR + TSD DETECTION")
    print("=" * 70)

    try:
        records = read_fasta(path)
    except FileNotFoundError:
        print(f"\nFile not found: {path}")
        print("Put it next to main2.py or run: python main2.py <file.fna>")
        return

    if not records:
        print(f"\nNo sequences found in {path}")
        return

    print(f"\nFile: {path}")
    print("Sequences read:", len(records))

    for name, dna in records:
        analyse(name, dna)


if __name__ == "__main__":
    main()