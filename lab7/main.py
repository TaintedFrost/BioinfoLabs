import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------
MIN_LEN, MAX_LEN = 1000, 3000          # allowed length of the DNA to digest
COMPLEMENT = str.maketrans("ACGT", "TGCA")

LADDER = (100, 200, 300, 500, 750, 1000, 1500, 2000, 3000)   # marker lane (bp)
GEL_TOP, GEL_BOTTOM = 3000, 50         # size range shown on the gel (bp)

# Input 1: recognized sequences (NAME=SITE, ^ = where the top strand is cut).
# EcoRI, BamHI and HindIII are the ones from the board. DraI and SspI are
# AT-rich sites: the SARS-CoV-2 genome is AT-rich, so they are much more
# likely to cut it than GC-rich sites such as SmaI (CCC^GGG).
DEFAULT_ENZYMES = (
    "EcoRI=G^AATTC",
    "BamHI=G^GATCC",
    "HindIII=A^AGCTT",
    #"DraI=TTT^AAA",
    #"SspI=AAT^ATT",
)


# --------------------------------------------------------------------------
# Enzymes
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Enzyme:
    name: str
    site: str   # recognized sequence, 5'->3' on the top strand
    cut: int    # the top strand is cut AFTER this many bases of the site

    @property
    def cleavage(self) -> str:
        """Cleavage site written with ^ (e.g. G^AATTC)."""
        return self.site[: self.cut] + "^" + self.site[self.cut:]

    @property
    def bottom_cut(self) -> int:
        """Where the bottom strand is cut (top-strand coordinates): mirror image
        of the top cut, because the sites are palindromic."""
        return len(self.site) - self.cut

    @property
    def ends(self) -> str:
        if self.cut < self.bottom_cut:
            return f"5' overhang ({self.site[self.cut:self.bottom_cut]})"
        if self.cut > self.bottom_cut:
            return f"3' overhang ({self.site[self.bottom_cut:self.cut]})"
        return "blunt"


def parse_enzyme(spec: str) -> Enzyme:
    """'EcoRI=G^AATTC' -> Enzyme('EcoRI', 'GAATTC', 1)."""
    name, sep, marked = spec.partition("=")
    name, marked = name.strip(), marked.strip().upper()
    if not sep or not name or marked.count("^") != 1:
        raise ValueError(f"Bad enzyme '{spec}'. Use NAME=SITE with one ^, e.g. EcoRI=G^AATTC")
    cut = marked.index("^")
    site = marked.replace("^", "")
    if not re.fullmatch(r"[ACGT]+", site):
        raise ValueError(f"{name}: the recognized sequence must only use A, C, G, T.")
    if not 0 < cut < len(site):
        raise ValueError(f"{name}: the ^ must be inside the recognized sequence.")
    if site != reverse_complement(site):
        raise ValueError(f"{name}: only palindromic recognition sites are supported.")
    return Enzyme(name, site, cut)


# --------------------------------------------------------------------------
# Sequence input (pure helpers)
# --------------------------------------------------------------------------
def complement(seq: str) -> str:
    """Complementary strand, written left to right (so it reads 3'->5')."""
    return seq.translate(COMPLEMENT)


def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]


def locate(path_str: str) -> Path:
    """Look for the file in the current folder, then next to main.py."""
    path = Path(path_str)
    for candidate in (path, Path(__file__).resolve().parent / path):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"'{path_str}' not found. Put it in the same folder as main.py or use --fasta PATH."
    )


# A line of DNA: only nucleotide letters (A, C, G, T and the IUPAC codes such as N).
DNA_LINE = re.compile(r"[ACGTNRYKMSWBDHV]+", re.IGNORECASE)


def read_fasta(path: Path) -> tuple[str, str]:
    """Read a DNA file -> (name, sequence in upper case).

    The '>' header line is optional, so these all work:
      * a normal FASTA file ('>MN908947.3 ...' then the sequence)
      * a file where the header line was deleted (only the sequence)
      * a file where just the '>' was deleted (header text left on line 1)
    Only the first record is used."""
    name, parts = path.name, []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():   # utf-8-sig drops a BOM
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if parts:
                break                       # start of a 2nd record: stop
            name = line[1:].strip() or name
        elif DNA_LINE.fullmatch(line):
            parts.append(line)
        elif not parts:
            name = line                     # header text whose '>' was erased
        else:
            raise ValueError(f"{path.name}: unexpected non-DNA line inside the sequence: '{line[:40]}'")
    if not parts:
        raise ValueError(f"{path.name}: no DNA sequence found.")
    return name, "".join(parts).upper()


def select_window(seq: str, start: int, length: int) -> str:
    """Nucleotides start .. start+length-1 (1-based) of the full sequence."""
    if start < 1:
        raise ValueError("--start must be 1 or greater.")
    window = seq[start - 1: start - 1 + length]
    if not MIN_LEN <= len(window) <= MAX_LEN:
        raise ValueError(
            f"The selected DNA has {len(window)} nt; the exercise needs between "
            f"{MIN_LEN} and {MAX_LEN}. Change --start / --length."
        )
    return window


# --------------------------------------------------------------------------
# Digestion
# --------------------------------------------------------------------------
def find_sites(seq: str, enzyme: Enzyme) -> list[int]:
    """0-based start of every occurrence of the recognized sequence
    (lookahead, so overlapping occurrences are found too)."""
    return [m.start() for m in re.finditer(f"(?={enzyme.site})", seq)]


def top_cuts(seq: str, enzyme: Enzyme) -> list[int]:
    """Cleavage positions: number of nucleotides to the left of each cut."""
    return [start + enzyme.cut for start in find_sites(seq, enzyme)]


def fragments(seq_len: int, cuts: list[int]) -> list[tuple[int, int]]:
    """(first nt, last nt), 1-based, of every piece of the linear DNA."""
    bounds = [0, *sorted(set(cuts)), seq_len]
    return [(a + 1, b) for a, b in zip(bounds, bounds[1:])]


def lengths(frags: list[tuple[int, int]]) -> list[int]:
    return [last - first + 1 for first, last in frags]


def cut_diagram(seq: str, enzyme: Enzyme, start: int, flank: int = 4) -> str:
    """Draw one cut like on the board: the two fragments with their sticky ends."""
    comp = complement(seq)
    w0 = max(0, start - flank)
    w1 = min(len(seq), start + len(enzyme.site) + flank)
    top_cut, bot_cut = start + enzyme.cut, start + enzyme.bottom_cut
    left_width = max(top_cut, bot_cut) - w0
    right_from = min(top_cut, bot_cut)

    left_top = seq[w0:top_cut].ljust(left_width)
    left_bot = comp[w0:bot_cut].ljust(left_width)
    right_top = " " * (top_cut - right_from) + seq[top_cut:w1]
    right_bot = " " * (bot_cut - right_from) + comp[bot_cut:w1]
    return (f"      5'-...{left_top}   {right_top}...-3'\n"
            f"      3'-...{left_bot}   {right_bot}...-5'")


# --------------------------------------------------------------------------
# Gel electrophoresis simulation
# --------------------------------------------------------------------------
# Small DNA runs faster: the distance travelled is roughly proportional to
# -log10(size). So the vertical position of a band is linear in log10(size).
def band_row(size: int, rows: int) -> int:
    size = min(max(size, GEL_BOTTOM), GEL_TOP)
    span = math.log10(GEL_TOP) - math.log10(GEL_BOTTOM)
    return round((math.log10(GEL_TOP) - math.log10(size)) / span * (rows - 1))


def render_gel(lanes: list[tuple[str, list[int]]], rows: int = 28) -> str:
    width = 9
    ladder_rows = {band_row(s, rows): s for s in LADDER}
    lane_rows = [{band_row(s, rows) for s in sizes} for _, sizes in lanes]

    out = ["bp      " + "".join(name[:8].center(width) for name, _ in lanes),
           "        " + "".join("[___]".center(width) for _ in lanes)]
    for r in range(rows):
        axis = f"{ladder_rows[r]:>5} -" if r in ladder_rows else "      |"
        cells = "".join(("=====" if r in present else "").center(width) for present in lane_rows)
        out.append(axis + " " + cells)
    out.append("      +" + "-" * (width * len(lanes) + 1))
    return "\n".join(out)


def save_gel_png(lanes: list[tuple[str, list[int]]], path: str) -> bool:
    """Draw the gel as a picture. Returns False if matplotlib is not installed."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    fig, ax = plt.subplots(figsize=(1.0 + 0.9 * len(lanes), 6))
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")
    for i, (_, sizes) in enumerate(lanes):
        total = sum(sizes)
        for size in sizes:
            shade = 0.35 + 0.65 * math.sqrt(size / total)     # more DNA -> brighter band
            ax.hlines(min(max(size, GEL_BOTTOM), GEL_TOP), i - 0.35, i + 0.35,
                      colors=[(shade, shade, shade)], linewidth=4)
    ax.set_yscale("log")
    ax.set_ylim(GEL_BOTTOM, GEL_TOP * 1.15)
    ax.set_xlim(-0.7, len(lanes) - 0.3)
    ax.set_yticks(LADDER)
    ax.set_yticklabels([str(s) for s in LADDER])
    ax.minorticks_off()
    ax.xaxis.tick_top()
    ax.set_xticks(range(len(lanes)))
    ax.set_xticklabels([name for name, _ in lanes], rotation=45, ha="left")
    ax.tick_params(colors="white")
    ax.set_ylabel("fragment size (bp)", color="white")
    for spine in ax.spines.values():
        spine.set_color("white")
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor="black")
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
def report(header: str, seq: str, start: int, enzymes: list[Enzyme], png: str | None) -> None:
    n = len(seq)
    rule = "=" * 78
    print(rule)
    print("RESTRICTION DIGEST AND GEL ELECTROPHORESIS (linear DNA)")
    print(rule)
    print(f"Source       : {header}")
    print(f"DNA digested : nt {start}-{start + n - 1} of the file ({n} nucleotides)")
    print("Positions    : numbered 1..{} inside the digested DNA".format(n))
    other = n - sum(seq.count(b) for b in "ACGT")
    if other:
        print(f"Warning      : {other} non-ACGT character(s) (e.g. N); no site can match there.")
    print()

    results = {e.name: top_cuts(seq, e) for e in enzymes}

    # Output 1: number of cleavages for each enzyme
    print("1) NUMBER OF CLEAVAGES")
    print(f"{'Enzyme':<10}{'Recognized seq.':<17}{'Cleavage site':<15}{'Ends':<22}Cleavages")
    print("-" * 78)
    for e in enzymes:
        print(f"{e.name:<10}{e.site:<17}{e.cleavage:<15}{e.ends:<22}{len(results[e.name])}")
    print()

    # Output 2: cleavage positions and fragment lengths
    print("2) CLEAVAGE POSITIONS AND FRAGMENT LENGTHS")
    print("-" * 78)
    for e in enzymes:
        cuts = results[e.name]
        frags = fragments(n, cuts)
        print(f"{e.name} ({e.cleavage}): {len(cuts)} cleavage(s)")
        if cuts:
            print(f"   cut between nt {' | '.join(f'{c}/{c + 1}' for c in cuts)}")
        for i, (a, b) in enumerate(frags, 1):
            print(f"   fragment {i}: nt {a:>5} - {b:<5}  {b - a + 1:>5} bp")
        if cuts:
            print("   first cut (both strands, sticky/blunt ends):")
            print(cut_diagram(seq, e, find_sites(seq, e)[0]))
        else:
            print("   no recognition site: the DNA is not cut.")
        print()

    all_cuts = sorted({c for cuts in results.values() for c in cuts})
    all_frags = fragments(n, all_cuts)
    print(f"All {len(enzymes)} enzymes together: {len(all_cuts)} cleavage(s)")
    for i, (a, b) in enumerate(all_frags, 1):
        print(f"   fragment {i}: nt {a:>5} - {b:<5}  {b - a + 1:>5} bp")
    print()

    # Output 3: gel electrophoresis
    lanes = [("Ladder", list(LADDER)), ("Uncut", [n])]
    lanes += [(e.name, lengths(fragments(n, results[e.name]))) for e in enzymes]
    lanes += [("All", lengths(all_frags))]

    print("3) GEL ELECTROPHORESIS SIMULATION")
    print("-" * 78)
    print(f"{len(enzymes)} enzymes -> {len(lanes)} lanes: ladder, uncut DNA, one lane per enzyme, all together")
    print("(top = wells / large fragments, bottom = small fragments; log scale)\n")
    print(render_gel(lanes))
    tiny = sorted({s for _, sizes in lanes[1:] for s in sizes if s < GEL_BOTTOM})
    if tiny:
        print(f"\nFragments below {GEL_BOTTOM} bp ({', '.join(map(str, tiny))} bp) run off the "
              f"bottom of the gel and are drawn on the last row.")
    if png:
        saved = save_gel_png(lanes, png)
        print(f"\nGel picture saved to {png}" if saved
              else "\n(Install matplotlib to also get a picture of the gel: pip install matplotlib)")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Restriction digest of covid.fna + gel simulation.")
    p.add_argument("--fasta", default="covid.fna", help="FASTA file with the DNA (default: covid.fna)")
    p.add_argument("--start", type=int, default=1, help="first nucleotide of the window (default: 1)")
    p.add_argument("--length", type=int, default=3000,
                   help=f"window length, {MIN_LEN}-{MAX_LEN} (default: 3000)")
    p.add_argument("--enzyme", action="append", metavar="NAME=SITE",
                   help="enzyme as NAME=SITE with ^ at the cut, e.g. EcoRI=G^AATTC "
                        "(repeat for several; default: 5 built-in enzymes)")
    p.add_argument("--no-png", action="store_true", help="do not save gel.png")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    try:
        enzymes = [parse_enzyme(s) for s in (args.enzyme or DEFAULT_ENZYMES)]
        header, genome = read_fasta(locate(args.fasta))
        seq = select_window(genome, args.start, args.length)
    except (ValueError, OSError) as err:
        sys.exit(f"Error: {err}")
    report(header, seq, args.start, enzymes, None if args.no_png else "gel.png")


if __name__ == "__main__":
    main()