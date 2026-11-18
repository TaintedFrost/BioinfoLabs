#!/usr/bin/env python3
"""
Sliding-window DNA patterns for gene-promoter sequences (PromKappa-style).

Information content follows the formula from the board:

    s = {x_1, ..., x_|s|}
    sigma(s) = 100 - ( sum_{v=1}^{|s|-1} [ sum_{i=1}^{|s|-v} f(x_i, x_{i+v}) / (|s|-v) * 100 ] / (|s|-1) )
    f(x_i, x_{i+v}) = 1 if x_i == x_{i+v} else 0

i.e. for every shift v the percentage of positions where the sequence equals
its own copy shifted by v is computed, those percentages are averaged over all
shifts (that average is the Kappa index of coincidence, the curve PromKappa
draws) and sigma is 100 minus that average.

For every window of length L the program computes
  * C+G %                       cg_content
  * CpG observed/expected       cpg_obs_exp
  * Kappa IC (coincidence) %    kappa_coincidence
  * sigma = 100 - Kappa IC      information_content   <- default information measure
  * Shannon information (bits)  shannon_information   (optional)
then plots the pattern, computes its centers of weight, and plots the
per-window centers on a second chart.

Usage
  python promoter_patterns.py                        # built-in test sequence S, L=30
  python promoter_patterns.py --fasta promoter.fa    # your own promoter (FASTA / plain text)
  python promoter_patterns.py --ic coincidence       # plot Kappa IC (what PromKappa draws)
  python promoter_patterns.py --ic shannon           # Shannon information instead
  python promoter_patterns.py --fasta p.fa --compare reference.csv
"""
import argparse
import csv
import math
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# S = CGGACTGATCTATCT (A x 27) CGTAGCATCTATCGATCTATCTAGCGATCTA
TEST_SEQ = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA"


# ----------------------------------------------------------------- input ---
def read_sequence(path):
    """Read a FASTA / plain-text file and return an upper-case ACGT string."""
    with open(path) as fh:
        lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith(">")]
    seq = "".join(lines).upper().replace("U", "T")
    return "".join(b for b in seq if b in "ACGT")


def sliding_windows(seq, L, step=1):
    """Yield (start_index_0based, window) for every window of length L."""
    for i in range(0, len(seq) - L + 1, step):
        yield i, seq[i:i + L]


# ------------------------------------------------------- per-window values ---
def cg_content(window):
    """C+G percentage of a window."""
    return 100.0 * (window.count("C") + window.count("G")) / len(window)


def cpg_obs_exp(window):
    """CpG observed/expected ratio = (#CG * L) / (#C * #G); 0 if undefined."""
    c, g = window.count("C"), window.count("G")
    return window.count("CG") * len(window) / (c * g) if c and g else 0.0


def f(a, b):
    """Board's f: 1 if the two bases are equal, else 0."""
    return 1 if a == b else 0


def kappa_coincidence(window):
    """
    Kappa index of coincidence (percent): for each shift v = 1..n-1 the share of
    positions with x_i == x_{i+v}, averaged over all shifts.
    """
    n = len(window)
    total = 0.0
    for v in range(1, n):
        total += 100.0 * sum(f(window[i], window[i + v]) for i in range(n - v)) / (n - v)
    return total / (n - 1)


def information_content(window):
    """Board formula: sigma(s) = 100 - Kappa index of coincidence."""
    return 100.0 - kappa_coincidence(window)


def shannon_information(window, background=None):
    """Total information (bits) of the window composition vs. the background (uniform by default)."""
    n = len(window)
    bg = background or {b: 0.25 for b in "ACGT"}
    per_base = sum((v / n) * math.log2((v / n) / bg[b]) for b, v in Counter(window).items())
    return per_base * n


# ------------------------------------------------------------ centers -------
def center_of_weight(positions, weights):
    """Center of weight: sum(x * w) / sum(w)."""
    total = sum(weights)
    return sum(x * w for x, w in zip(positions, weights)) / total if total else float("nan")


def gc_mass_center(window, start):
    """Center of weight of ONE window: mean (1-based, sequence) position of its G and C bases."""
    pos = [start + i + 1 for i, b in enumerate(window) if b in "CG"]
    return sum(pos) / len(pos) if pos else float("nan")


# ----------------------------------------------------------------- pattern ---
def build_pattern(seq, L, step=1):
    rows = []
    for start, w in sliding_windows(seq, L, step):
        kc = kappa_coincidence(w)
        rows.append({
            "start0": start,                         # 0-based start (as PromKappa prints it)
            "start": start + 1,                      # 1-based first base
            "center": start + (L + 1) / 2,           # 1-based window middle
            "cg": cg_content(w),
            "cpg_oe": cpg_obs_exp(w),
            "sigma": 100.0 - kc,
            "kappa_ic": kc,
            "shannon_ic": shannon_information(w),
            "gc_center": gc_mass_center(w, start),
        })
    return rows


def composition(seq):
    n = len(seq)
    return {b: 100.0 * seq.count(b) / n for b in "ACGT"}


def self_check():
    """Values PromKappa shows for the test sequence S."""
    assert len(TEST_SEQ) == 73
    assert abs(cg_content(TEST_SEQ) - 27.3972602739726) < 1e-9
    comp = composition(TEST_SEQ)
    assert [int(comp[b]) for b in "ATCG"] == [52, 20, 16, 10]   # PromKappa truncates, it does not round
    # PromKappa labels the window 11b-41b with Kappa IC 73.69
    assert abs(kappa_coincidence(TEST_SEQ[11:41]) - 73.69) < 0.05


# ------------------------------------------------------------------- plots ---
MEASURES = {
    "sigma": ("sigma", "sigma = 100 - Kappa IC (information content)", "%"),
    "coincidence": ("kappa_ic", "Kappa IC (coincidence) %", "%"),
    "shannon": ("shannon_ic", "Shannon information (bits)", " bits"),
}


def plot_pattern(rows, L, col, label, out):
    x = [r["center"] for r in rows]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.plot(x, [r["cg"] for r in rows], color="tab:red", marker="o", ms=3, label="C+G %")
    if col == "shannon_ic":
        ax2 = ax.twinx()
        ax2.plot(x, [r[col] for r in rows], color="black", marker="o", ms=3, label=label)
        ax2.set_ylabel("bits")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
        ax.set_ylabel("C+G %")
    else:
        ax.plot(x, [r[col] for r in rows], color="black", marker="o", ms=3, label=label)
        ax.set_ylabel("percent")
        ax.legend(fontsize=8, loc="upper left")
    ax.set_xlabel("Window center (position in sequence)")
    ax.set_title(f"DNA pattern, sliding window L={L}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_centers(rows, cw_cg, cw_ic, short_label, out):
    pts = [(r["center"], r["gc_center"]) for r in rows if not math.isnan(r["gc_center"])]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.scatter(*zip(*pts), s=18, color="tab:blue", label="center of weight of each window (G/C bases)")
    ax.axhline(cw_cg, color="tab:red", ls="--", label=f"pattern center, C+G% = {cw_cg:.2f}")
    ax.axhline(cw_ic, color="black", ls=":", label=f"pattern center, {short_label} = {cw_ic:.2f}")
    ax.set_xlabel("Window center (position in sequence)")
    ax.set_ylabel("Center of weight (position)")
    ax.set_title("Centers of the patterns")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------- compare ---
def compare_with_reference(rows, ref_csv, col):
    """Compare with a reference CSV that has columns cg and <col> (or ic)."""
    with open(ref_csv) as fh:
        ref = list(csv.DictReader(fh))
    if len(ref) != len(rows):
        print(f"[compare] different number of windows: ours={len(rows)} reference={len(ref)}")
        return False
    rcol = col if col in ref[0] else "ic"
    d_cg = max(abs(float(a["cg"]) - b["cg"]) for a, b in zip(ref, rows))
    d_ic = max(abs(float(a[rcol]) - b[col]) for a, b in zip(ref, rows))
    print(f"[compare] max |diff|  C+G%: {d_cg:.4g}   {col}: {d_ic:.4g}")
    same = d_cg < 0.01 and d_ic < 0.01
    print("[compare] patterns are", "THE SAME" if same else "DIFFERENT")
    return same


# -------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", help="promoter sequence file (FASTA or plain text)")
    ap.add_argument("-L", "--window", type=int, default=30)
    ap.add_argument("--step", type=int, default=1)
    ap.add_argument("--ic", choices=list(MEASURES), default="sigma",
                    help="information measure used for the pattern (default: sigma, the board formula)")
    ap.add_argument("--prefix", default="pattern", help="output file prefix")
    ap.add_argument("--compare", help="reference CSV with columns cg and sigma/kappa_ic/ic")
    args = ap.parse_args()

    self_check()
    seq = read_sequence(args.fasta) if args.fasta else TEST_SEQ
    L = args.window
    if len(seq) < L:
        sys.exit(f"Sequence ({len(seq)} nt) is shorter than the window ({L}).")

    col, label, unit = MEASURES[args.ic]
    short = {"sigma": "sigma", "kappa_ic": "Kappa IC", "shannon_ic": "Shannon Ic"}[col]

    rows = build_pattern(seq, L, args.step)
    xs = [r["center"] for r in rows]
    cw_cg = center_of_weight(xs, [r["cg"] for r in rows])
    cw_ic = center_of_weight(xs, [r[col] for r in rows])

    comp = {b: int(v) for b, v in composition(seq).items()}
    vals = [r[col] for r in rows]
    print(f"sequence length: {len(seq)} nt, windows: {len(rows)}, L={L}, step={args.step}")
    print(f"CG = {cg_content(seq)} %")
    print("A = {A} %  T = {T} %  C = {C} %  G = {G} %   (truncated, as PromKappa shows them)".format(**comp))
    print(f"C+G% first / min / max window : {rows[0]['cg']:.2f} / {min(r['cg'] for r in rows):.2f} / {max(r['cg'] for r in rows):.2f}")
    print(f"{short} first / min / max window: {vals[0]:.2f} / {min(vals):.2f} / {max(vals):.2f}{unit}")
    print(f"center of weight, C+G% pattern : {cw_cg:.2f}")
    print(f"center of weight, {short} pattern: {cw_ic:.2f}")

    with open(f"{args.prefix}.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    plot_pattern(rows, L, col, label, f"{args.prefix}_pattern.png")
    plot_centers(rows, cw_cg, cw_ic, short, f"{args.prefix}_centers.png")
    print(f"saved {args.prefix}.csv, {args.prefix}_pattern.png, {args.prefix}_centers.png")

    if args.compare:
        compare_with_reference(rows, args.compare, col)


if __name__ == "__main__":
    main()