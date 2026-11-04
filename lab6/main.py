
import random

# ==================================================
# CONFIGURATION
# ==================================================

DNA_LENGTH = 300

TIR_LENGTH = 8
TSD_LENGTH = 4

MIN_INTERNAL_LENGTH = 25
MAX_INTERNAL_LENGTH = 70

random.seed(42)

COMPLEMENT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


# ==================================================
# DNA FUNCTIONS
# ==================================================

def generate_dna(length):
    return ''.join(random.choice("ATCG") for _ in range(length))


def reverse_complement(sequence):
    return ''.join(
        COMPLEMENT[base]
        for base in reversed(sequence)
    )


# ==================================================
# EXPERIMENT DEFINITIONS
#
# Coordinates are 0-based and inclusive.
#
# Format: (start, end)
# ==================================================

EXPERIMENTS = {

    "1. No overlap": [
        (20, 69),
        (85, 134),
        (150, 199),
        (215, 264)
    ],

    "2. Partial overlap (2 elements)": [
        (30, 85),
        (70, 125),
        (160, 209),
        (230, 279)
    ],

    "3. Nested overlap (3 elements)": [
        (40, 125),
        (50, 115),
        (60, 105),
        (200, 249)
    ],

    "4. Shared starting position": [
        (30, 79),
        (30, 94),
        (145, 194),
        (225, 274)
    ],

    "5. Shared ending position": [
        (80, 129),
        (65, 129),
        (160, 209),
        (225, 274)
    ],

    "6. Three-way partial overlap": [
        (25, 84),
        (45, 104),
        (65, 124),
        (190, 239)
    ]
}


# ==================================================
# DNA CONSTRAINT SYSTEM
#
# Ensures overlapping TIRs and TSDs are compatible.
#
# Relation 0 = identical bases
# Relation 1 = complementary bases
# ==================================================

class DNAConstraints:

    def __init__(self, length):

        self.parent = list(range(length))
        self.parity = [0] * length
        self.length = length

    def find(self, x):

        if self.parent[x] != x:

            old_parent = self.parent[x]

            root = self.find(old_parent)

            self.parity[x] ^= self.parity[old_parent]

            self.parent[x] = root

        return self.parent[x]

    def get_parity(self, x):

        self.find(x)

        return self.parity[x]

    def add_constraint(self, a, b, relation):

        root_a = self.find(a)
        parity_a = self.parity[a]

        root_b = self.find(b)
        parity_b = self.parity[b]

        if root_a == root_b:

            if (parity_a ^ parity_b) != relation:

                return False

            return True

        self.parent[root_a] = root_b

        self.parity[root_a] = (
            parity_a ^ parity_b ^ relation
        )

        return True

    def generate_sequence(self):

        # Assign a random nucleotide to each independent group
        root_bases = {}

        for i in range(self.length):

            root = self.find(i)

            if root not in root_bases:

                root_bases[root] = random.choice("ATCG")

        # Construct the final sequence
        sequence = []

        for i in range(self.length):

            root = self.find(i)

            base = root_bases[root]

            if self.parity[i] == 1:

                base = COMPLEMENT[base]

            sequence.append(base)

        return ''.join(sequence)


# ==================================================
# BUILD DNA WITH TIRs AND TSDs
# ==================================================

def create_dna(element_specs):

    constraints = DNAConstraints(DNA_LENGTH)

    original = []

    for start, end in element_specs:

        element_length = end - start + 1

        internal_length = (
            element_length - 2 * TIR_LENGTH
        )

        if not (
            MIN_INTERNAL_LENGTH
            <= internal_length
            <= MAX_INTERNAL_LENGTH
        ):

            raise ValueError(
                f"Invalid element length: {element_length}"
            )

        if (
            start < TSD_LENGTH
            or end + TSD_LENGTH >= DNA_LENGTH
        ):

            raise ValueError(
                "Not enough space for target site duplications!"
            )

        # ------------------------------------------
        # TSD constraints
        #
        # Left and right TSDs must be identical.
        # ------------------------------------------

        for i in range(TSD_LENGTH):

            left_position = start - TSD_LENGTH + i

            right_position = end + 1 + i

            if not constraints.add_constraint(
                left_position,
                right_position,
                0
            ):

                raise ValueError(
                    "Incompatible TSD arrangement!"
                )

        # ------------------------------------------
        # TIR constraints
        #
        # Left and right TIRs must be reverse
        # complements of each other.
        # ------------------------------------------

        for i in range(TIR_LENGTH):

            left_position = start + i

            right_position = end - i

            if not constraints.add_constraint(
                left_position,
                right_position,
                1
            ):

                raise ValueError(
                    "Incompatible TIR arrangement!"
                )

        original.append({
            "start": start + 1,
            "end": end + 1,
            "length": element_length,
            "internal_length": internal_length
        })

    # Generate DNA satisfying all constraints
    dna = constraints.generate_sequence()

    return dna, original


# ==================================================
# DETECT TRANSPOSABLE ELEMENTS
#
# Conditions:
# 1. TIRs must be reverse complements.
# 2. Internal DNA must be 25-70 bases.
# 3. Flanking TSDs must be identical.
# ==================================================

def detect_elements(dna):

    detected = []

    for start in range(TSD_LENGTH, len(dna)):

        # Left TIR
        left_tir = dna[
            start:start + TIR_LENGTH
        ]

        if len(left_tir) != TIR_LENGTH:
            continue

        expected_right = reverse_complement(left_tir)

        # Left TSD
        left_tsd = dna[
            start - TSD_LENGTH:start
        ]

        # Search possible internal lengths
        for internal_length in range(
            MIN_INTERNAL_LENGTH,
            MAX_INTERNAL_LENGTH + 1
        ):

            right_start = (
                start + TIR_LENGTH + internal_length
            )

            right_end = right_start + TIR_LENGTH

            if right_end + TSD_LENGTH > len(dna):
                continue

            right_tir = dna[
                right_start:right_end
            ]

            right_tsd = dna[
                right_end:right_end + TSD_LENGTH
            ]

            # Verify both terminal structures
            if (
                right_tir == expected_right
                and left_tsd == right_tsd
            ):

                detected.append((
                    start + 1,
                    right_end
                ))

    return sorted(set(detected))


# ==================================================
# CHECK OVERLAPPING ELEMENTS
# ==================================================

def check_overlaps(original):

    overlaps = []

    for i in range(len(original)):

        for j in range(i + 1, len(original)):

            a = original[i]
            b = original[j]

            overlap_start = max(
                a["start"],
                b["start"]
            )

            overlap_end = min(
                a["end"],
                b["end"]
            )

            if overlap_start <= overlap_end:

                overlaps.append({
                    "first": i + 1,
                    "second": j + 1,
                    "bases": overlap_end - overlap_start + 1
                })

    return overlaps


# ==================================================
# RUN ONE EXPERIMENT
# ==================================================

def run_experiment(name, specs):

    print("\n" + "=" * 70)
    print(name)
    print("=" * 70)

    dna, original = create_dna(specs)

    detected = detect_elements(dna)

    print("\nDNA sequence:")
    print(dna)

    print("\nDNA length:", len(dna), "bases")

    # ----------------------------------------------
    # Original elements
    # ----------------------------------------------

    print("\nORIGINAL TRANSPOSABLE ELEMENTS:")

    for i, element in enumerate(original, 1):

        print(
            f"TE {i}: Start = {element['start']}, "
            f"End = {element['end']}, "
            f"Length = {element['length']} bases"
        )

    # ----------------------------------------------
    # TIR and TSD structures
    # ----------------------------------------------

    print("\nSEQUENCE STRUCTURE:")

    for i, element in enumerate(original, 1):

        start = element["start"] - 1
        end = element["end"]

        left_tsd = dna[
            start - TSD_LENGTH:start
        ]

        left_tir = dna[
            start:start + TIR_LENGTH
        ]

        right_tir = dna[
            end - TIR_LENGTH:end
        ]

        right_tsd = dna[
            end:end + TSD_LENGTH
        ]

        print(f"\nTE {i}:")

        print("Left TSD: ", left_tsd)
        print("Left TIR: ", left_tir)
        print("Right TIR:", right_tir)
        print("Right TSD:", right_tsd)

    # ----------------------------------------------
    # Detection results
    # ----------------------------------------------

    print("\nDETECTED TRANSPOSABLE ELEMENTS:")

    for i, (start, end) in enumerate(detected, 1):

        print(
            f"Detection {i}: Start = {start}, "
            f"End = {end}, "
            f"Length = {end - start + 1} bases"
        )

    # ----------------------------------------------
    # Verification
    # ----------------------------------------------

    original_positions = {
        (e["start"], e["end"])
        for e in original
    }

    detected_positions = set(detected)

    true_positives = (
        original_positions & detected_positions
    )

    false_positives = (
        detected_positions - original_positions
    )

    missed = (
        original_positions - detected_positions
    )

    print("\nVERIFICATION:")

    for start, end in sorted(true_positives):

        print(f"CORRECT: TE at positions {start}-{end}")

    for start, end in sorted(missed):

        print(f"MISSED: TE at positions {start}-{end}")

    for start, end in sorted(false_positives):

        print(f"FALSE POSITIVE: {start}-{end}")

    # ----------------------------------------------
    # Overlap analysis
    # ----------------------------------------------

    overlaps = check_overlaps(original)

    print("\nOVERLAP ANALYSIS:")

    if overlaps:

        for overlap in overlaps:

            print(
                f"TE {overlap['first']} and "
                f"TE {overlap['second']} overlap by "
                f"{overlap['bases']} bases."
            )

    else:

        print("No overlapping elements.")

    # ----------------------------------------------
    # Summary
    # ----------------------------------------------

    print("\nEXPERIMENT SUMMARY:")

    print("Original elements:", len(original))
    print("Total detections:", len(detected))
    print("Correct detections:", len(true_positives))
    print("False positives:", len(false_positives))
    print("Missed elements:", len(missed))

    return {
        "name": name,
        "original": len(original),
        "detected": len(detected),
        "correct": len(true_positives),
        "false": len(false_positives),
        "missed": len(missed)
    }


# ==================================================
# MAIN PROGRAM
# ==================================================

def main():

    results = []

    print("=" * 70)
    print("TRANSPOSABLE ELEMENT DETECTION")
    print("TIR + TSD + OVERLAP SIMULATION")
    print("=" * 70)

    print("\nRunning six experiments...")

    for name, specs in EXPERIMENTS.items():

        result = run_experiment(name, specs)

        results.append(result)

    # ----------------------------------------------
    # Final comparison
    # ----------------------------------------------

    print("\n\n")

    print("=" * 85)
    print("FINAL COMPARISON OF ALL EXPERIMENTS")
    print("=" * 85)

    print(
        f"{'Experiment':<38}"
        f"{'Original':<11}"
        f"{'Detected':<11}"
        f"{'Correct':<10}"
        f"{'False +':<10}"
        f"{'Missed':<10}"
    )

    print("-" * 85)

    for result in results:

        print(
            f"{result['name']:<38}"
            f"{result['original']:<11}"
            f"{result['detected']:<11}"
            f"{result['correct']:<10}"
            f"{result['false']:<10}"
            f"{result['missed']:<10}"
        )

    print("-" * 85)

    print("\nAll six experiments completed.")


if __name__ == "__main__":
    main()
