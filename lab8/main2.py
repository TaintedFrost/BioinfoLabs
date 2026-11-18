#!/usr/bin/env python3
"""
Global sequence alignment (Needleman-Wunsch) of two DNA sequences.

Scoring (change with the options): match +1, mismatch -1, gap -1.

Usage
  python global_alignment.py                                  # the two lab sequences
  python global_alignment.py --s1 ACGT --s2 AGT               # your own sequences
  python global_alignment.py --match 2 --mismatch -1 --gap -2 # other scoring
  python global_alignment.py --max-alignments 5               # show up to 5 optimal alignments
"""
import argparse

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"


def score_matrix(a, b, match, mismatch, gap):
    """Fill the Needleman-Wunsch score matrix F, size (len(a)+1) x (len(b)+1)."""
    n, m = len(a), len(b)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if a[i - 1] == b[j - 1] else mismatch
            F[i][j] = max(F[i - 1][j - 1] + s,   # diagonal: align a[i-1] with b[j-1]
                          F[i - 1][j] + gap,     # up:       a[i-1] against a gap
                          F[i][j - 1] + gap)     # left:     b[j-1] against a gap
    return F


def traceback(F, a, b, match, mismatch, gap, limit=1):
    """Walk back from the bottom-right corner; return up to `limit` optimal alignments."""
    results = []

    def walk(i, j, top, bottom):
        if len(results) >= limit:
            return
        if i == 0 and j == 0:
            results.append((top[::-1], bottom[::-1]))
            return
        if i > 0 and j > 0:
            s = match if a[i - 1] == b[j - 1] else mismatch
            if F[i][j] == F[i - 1][j - 1] + s:
                walk(i - 1, j - 1, top + a[i - 1], bottom + b[j - 1])
        if i > 0 and F[i][j] == F[i - 1][j] + gap:
            walk(i - 1, j, top + a[i - 1], bottom + "-")
        if j > 0 and F[i][j] == F[i][j - 1] + gap:
            walk(i, j - 1, top + "-", bottom + b[j - 1])

    walk(len(a), len(b), "", "")
    return results


def global_align(a, b, match=1, mismatch=-1, gap=-1, limit=1):
    """Return (optimal score, list of optimal alignments, score matrix)."""
    F = score_matrix(a, b, match, mismatch, gap)
    return F[len(a)][len(b)], traceback(F, a, b, match, mismatch, gap, limit), F


def print_matrix(F, a, b):
    print("Score matrix:")
    print("      -  " + "  ".join(f"{c:>3}" for c in b))
    for i, row in enumerate(F):
        label = "-" if i == 0 else a[i - 1]
        print(f"  {label}  " + "  ".join(f"{v:>3}" for v in row))


def show_alignment(top, bottom):
    mid = "".join("|" if x == y else (" " if "-" in (x, y) else ".") for x, y in zip(top, bottom))
    print(f"  s1: {top}")
    print(f"      {mid}")
    print(f"  s2: {bottom}")
    matches = mid.count("|")
    mism = mid.count(".")
    gaps = mid.count(" ")
    print(f"  length {len(top)}, matches {matches}, mismatches {mism}, gaps {gaps}, "
          f"identity {100 * matches / len(top):.1f} %")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s1", default=S1)
    ap.add_argument("--s2", default=S2)
    ap.add_argument("--match", type=int, default=1)
    ap.add_argument("--mismatch", type=int, default=-1)
    ap.add_argument("--gap", type=int, default=-1)
    ap.add_argument("--max-alignments", type=int, default=1,
                    help="how many optimal alignments to print (default 1)")
    args = ap.parse_args()

    a, b = args.s1.upper(), args.s2.upper()
    score, alignments, F = global_align(a, b, args.match, args.mismatch, args.gap, args.max_alignments)

    print(f"s1 = {a}  ({len(a)} nt)")
    print(f"s2 = {b}  ({len(b)} nt)")
    print(f"scoring: match {args.match:+d}, mismatch {args.mismatch:+d}, gap {args.gap:+d}\n")
    print_matrix(F, a, b)
    print(f"\nOptimal global alignment score: {score}\n")
    for k, (top, bottom) in enumerate(alignments, 1):
        print(f"Alignment {k}:")
        show_alignment(top, bottom)
        print()


if __name__ == "__main__":
    main()