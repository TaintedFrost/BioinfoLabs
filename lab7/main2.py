import argparse
import math
import re
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------
COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Whole genomes are 19-30 kb, so the ladder goes much higher than in main.py.
# Only the sizes up to the longest DNA are drawn (see ladder_for).
LADDER = (100, 250, 500, 1000, 2000, 3000, 5000, 8000, 12000, 20000, 30000)
GEL_BOTTOM = 50          # fragments smaller than this run off the gel
LANE_W = 5               # width of one lane in the text gel
BAND = "==="
SHORT_LEN = 3000       # length kept when the answer is "y"

DEFAULT_ENZYMES = (
    "EcoRI=G^AATTC",
    "BamHI=G^GATCC",
    "HindIII=A^AGCTT",
    #"DraI=TTT^AAA",
    #"SspI=AAT^ATT",
)

# Lane colours in gel2.png: genome A cyan, genome B amber, ladder white.
COLORS = {"A": (0.55, 0.85, 1.0), "B": (1.0, 0.75, 0.35), "": (1.0, 1.0, 1.0)}


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
        return self.site[: self.cut] + "^" + self.site[self.cut:]

    @property
    def ends(self) -> str:
        bottom_cut = len(self.site) - self.cut      # palindromic site: mirror image
        if self.cut < bottom_cut:
            return f"5' overhang ({self.site[self.cut:bottom_cut]})"
        if self.cut > bottom_cut:
            return f"3' overhang ({self.site[bottom_cut:self.cut]})"
        return "blunt"


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


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
# Reading the genomes
# --------------------------------------------------------------------------
DNA_LINE = re.compile(r"[ACGTNRYKMSWBDHV]+", re.IGNORECASE)


@dataclass(frozen=True)
class Genome:
    tag: str       # "A" or "B": the label on the gel
    name: str      # file name without extension
    header: str    # FASTA header (or the file name if there is none)
    seq: str       # the DNA that is digested (maybe only the first SHORT_LEN nt)
    full_len: int  # length of the whole sequence in the file


def locate(path_str: str) -> Path:
    """Look for the file in the current folder, then next to main2.py."""
    path = Path(path_str)
    for candidate in (path, Path(__file__).resolve().parent / path):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"'{path_str}' not found. Put it in the same folder as main2.py or use --fasta1/--fasta2."
    )


def read_fasta(path: Path) -> tuple[str, str]:
    """Read a DNA file -> (header, sequence in upper case). The '>' line is optional:
    a plain sequence file, or a header whose '>' was erased, also works.
    Only the first record is used."""
    name, parts = path.name, []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if parts:
                break
            name = line[1:].strip() or name
        elif DNA_LINE.fullmatch(line):
            parts.append(line)
        elif not parts:
            name = line
        else:
            raise ValueError(f"{path.name}: unexpected non-DNA line: '{line[:40]}'")
    if not parts:
        raise ValueError(f"{path.name}: no DNA sequence found.")
    return name, "".join(parts).upper()


def load_genome(tag: str, path_str: str, limit: int | None = None) -> Genome:
    """Read a genome; with `limit`, keep only its first `limit` nucleotides."""
    path = locate(path_str)
    header, seq = read_fasta(path)
    return Genome(tag, path.stem, header, seq[:limit] if limit else seq, len(seq))


# --------------------------------------------------------------------------
# Digestion
# --------------------------------------------------------------------------
def top_cuts(seq: str, enzyme: Enzyme) -> list[int]:
    """Cleavage positions: number of nucleotides to the left of each cut
    (overlapping occurrences of the site are found too)."""
    return [m.start() + enzyme.cut for m in re.finditer(f"(?={enzyme.site})", seq)]


def fragment_lengths(seq_len: int, cuts: list[int]) -> list[int]:
    bounds = [0, *sorted(set(cuts)), seq_len]
    return [b - a for a, b in zip(bounds, bounds[1:])]


def all_cuts(seq: str, enzymes: list[Enzyme]) -> list[int]:
    return sorted({c for e in enzymes for c in top_cuts(seq, e)})


# --------------------------------------------------------------------------
# Gel: lanes are grouped (a group = the two lanes A and B of one digest)
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Lane:
    tag: str                    # "A", "B" or "" for the ladder
    sizes: tuple[int, ...]      # fragment sizes in bp


Group = tuple[str, list[Lane]]


def ladder_for(longest: int) -> tuple[int, ...]:
    """Ladder bands up to the longest DNA, so a 3000 nt gel is not squashed
    by a 30000 bp marker."""
    return tuple(s for s in LADDER if s <= longest) or (LADDER[0],)


def build_groups(genomes: list[Genome], enzymes: list[Enzyme]) -> list[Group]:
    def pair(cuts_of) -> list[Lane]:
        return [Lane(g.tag, tuple(fragment_lengths(len(g.seq), cuts_of(g)))) for g in genomes]

    ladder = ladder_for(max(len(g.seq) for g in genomes))
    groups: list[Group] = [("Lad", [Lane("", ladder)]), ("Uncut", pair(lambda g: []))]
    groups += [(e.name, pair(lambda g, e=e: top_cuts(g.seq, e))) for e in enzymes]
    groups.append(("All", pair(lambda g: all_cuts(g.seq, enzymes))))
    return groups


def gel_top(groups: list[Group]) -> int:
    return max(max(lane.sizes) for _, lanes in groups for lane in lanes)


def band_row(size: int, rows: int, top: int) -> int:
    """Small DNA runs farther: the row is linear in log10(size)."""
    size = min(max(size, GEL_BOTTOM), top)
    span = math.log10(top) - math.log10(GEL_BOTTOM)
    return round((math.log10(top) - math.log10(size)) / span * (rows - 1))


def render_gel(groups: list[Group], rows: int = 40) -> str:
    top = gel_top(groups)
    ladder = groups[0][1][0].sizes
    ladder_rows = {band_row(s, rows, top): s for s in ladder}
    margin = " " * 8

    header_names = margin + " ".join(name.center(LANE_W * len(lanes)) for name, lanes in groups)
    header_tags = margin + " ".join("".join(l.tag.center(LANE_W) for l in lanes) for _, lanes in groups)
    wells = margin + " ".join("".join("[_]".center(LANE_W) for _ in lanes) for _, lanes in groups)

    out = [header_names, header_tags, wells]
    for r in range(rows):
        axis = f"{ladder_rows[r]:>5} -" if r in ladder_rows else "      |"
        cells = " ".join(
            "".join((BAND if r in {band_row(s, rows, top) for s in l.sizes} else "").center(LANE_W)
                    for l in lanes)
            for _, lanes in groups
        )
        out.append(axis + " " + cells)
    width = sum(LANE_W * len(lanes) for _, lanes in groups) + len(groups)
    out.append("      +" + "-" * (width + 1))
    return "\n".join(out)


def save_gel_png(groups: list[Group], genomes: list[Genome], path: str) -> bool:
    """Draw the gel as a picture. Returns False if matplotlib is not installed."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        return False

    top = gel_top(groups)
    n_lanes = sum(len(lanes) for _, lanes in groups)
    fig, ax = plt.subplots(figsize=(1.5 + 0.55 * n_lanes, 7))
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    x, ticks, labels = 0.0, [], []
    for name, lanes in groups:
        first = x
        for lane in lanes:
            total = sum(lane.sizes)
            r, g, b = COLORS[lane.tag]
            for size in lane.sizes:
                shade = 0.35 + 0.65 * math.sqrt(size / total)      # more DNA -> brighter
                ax.hlines(min(max(size, GEL_BOTTOM), top), x - 0.4, x + 0.4,
                          colors=[(r * shade, g * shade, b * shade)], linewidth=3.5)
            ticks.append(x)
            labels.append(lane.tag)
            x += 1.0
        ax.annotate(name, xy=((first + x - 1) / 2, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(0, 6), textcoords="offset points", ha="center", va="bottom",
                    color="white", fontsize=9, rotation=45)
        x += 0.6                                                     # gap between groups

    ax.set_yscale("log")
    ax.set_ylim(GEL_BOTTOM, top * 1.15)
    ax.set_xlim(-0.8, x - 0.8)
    ladder = list(groups[0][1][0].sizes)
    ax.set_yticks(ladder)
    ax.set_yticklabels([str(s) for s in ladder])
    ax.minorticks_off()
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.tick_params(colors="white")
    ax.set_ylabel("fragment size (bp)", color="white")
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.legend(handles=[Patch(color=COLORS[g.tag], label=f"{g.tag} = {g.name}") for g in genomes],
              loc="lower center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False,
              labelcolor="white")
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor="black")
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
def numbers_block(label: str, numbers: list[int]) -> str:
    return textwrap.fill(", ".join(map(str, numbers)), width=78,
                         initial_indent=" " * 6 + label, subsequent_indent=" " * (6 + len(label)))


def report(genomes: list[Genome], enzymes: list[Enzyme], png: str | None) -> None:
    rule = "=" * 78
    print(rule)
    print("TWO GENOMES, THE SAME ENZYMES, ONE GEL (linear DNA)")
    print(rule)
    for g in genomes:
        gc = 100 * (g.seq.count("G") + g.seq.count("C")) / len(g.seq)
        size = (f"{len(g.seq):,} nt" if len(g.seq) == g.full_len
                else f"first {len(g.seq):,} of {g.full_len:,} nt")
        print(f"Genome {g.tag}: {g.name} - {size}, GC {gc:.1f} %")
        other = len(g.seq) - sum(g.seq.count(b) for b in "ACGT")
        if other:
            print(f"          warning: {other} non-ACGT character(s) (e.g. N); no site can match there.")
    print("Positions are numbered 1..N along each genome.\n")

    cuts = {(g.tag, e.name): top_cuts(g.seq, e) for g in genomes for e in enzymes}

    # Output 1: number of cleavages
    print("1) NUMBER OF CLEAVAGES")
    names = "".join(f"{g.tag}: {g.name[:7]:<10}" for g in genomes)
    print(f"{'Enzyme':<10}{'Recognized seq.':<17}{'Cleavage site':<15}{'Ends':<22}{names}")
    print("-" * 78)
    for e in enzymes:
        counts = "".join(f"{len(cuts[(g.tag, e.name)]):<13}" for g in genomes)
        print(f"{e.name:<10}{e.site:<17}{e.cleavage:<15}{e.ends:<22}{counts}")
    print()

    # Output 2: cleavage positions and fragment lengths
    print("2) CLEAVAGE POSITIONS AND FRAGMENT LENGTHS")
    print("-" * 78)
    for e in enzymes:
        print(f"{e.name} ({e.cleavage})")
        for g in genomes:
            c = cuts[(g.tag, e.name)]
            print(f"   {g.tag} {g.name}: {len(c)} cleavage(s), {len(c) + 1} fragment(s)")
            if c:
                print(numbers_block("cut after nt : ", c))
                print(numbers_block("fragments bp : ", fragment_lengths(len(g.seq), c)))
        print()
    print(f"All {len(enzymes)} enzymes together")
    for g in genomes:
        c = all_cuts(g.seq, enzymes)
        print(f"   {g.tag} {g.name}: {len(c)} cleavage(s), {len(c) + 1} fragment(s)")
        if c:
            print(numbers_block("fragments bp : ", fragment_lengths(len(g.seq), c)))
    print()

    # Output 3: the gel
    groups = build_groups(genomes, enzymes)
    n_lanes = sum(len(lanes) for _, lanes in groups)
    print("3) GEL ELECTROPHORESIS SIMULATION (both genomes on the same gel)")
    print("-" * 78)
    print(f"{len(enzymes)} enzymes -> {n_lanes} lanes: ladder, uncut A and B, a pair of lanes (A, B) "
          f"per enzyme, all together")
    print("A = " + genomes[0].name + "   B = " + genomes[1].name)
    print("(top = wells / large fragments, bottom = small fragments; log scale)\n")
    print(render_gel(groups))

    tiny = sorted({s for _, lanes in groups[1:] for l in lanes for s in l.sizes if s < GEL_BOTTOM})
    if tiny:
        shown = ", ".join(map(str, tiny[:12])) + (", ..." if len(tiny) > 12 else "")
        print(f"\nFragments smaller than {GEL_BOTTOM} bp ({shown}) run off the bottom of the gel "
              f"and are drawn on the last row.")
    if png:
        print(f"\nGel picture saved to {png}" if save_gel_png(groups, genomes, png)
              else "\n(Install matplotlib to also get a picture: pip install matplotlib)")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def ask_short() -> bool:
    """y/n question at the start. True = first SHORT_LEN nt, False = full genomes."""
    prompt = (f"Keep only the first {SHORT_LEN} nucleotides of each file? "
              f"(y = first {SHORT_LEN}, n = full genomes): ")
    while True:
        try:
            answer = input(prompt).strip().lower()
        except EOFError:                      # no keyboard (e.g. an output panel)
            print("\nNo keyboard input available here, so the full genomes are used. "
                  "Run it in a terminal, or use --short / --full.")
            return False
        if answer in ("y", "yes"):
            return True
        if answer in ("n", "no"):
            return False
        print("Please type y or n.")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Digest two genomes and compare them on one gel.")
    p.add_argument("--fasta1", default="covid.fna", help="genome A (default: covid.fna)")
    p.add_argument("--fasta2", default="ebola.fna", help="genome B (default: ebola.fna)")
    p.add_argument("--enzyme", action="append", metavar="NAME=SITE",
                   help="enzyme as NAME=SITE with ^ at the cut, e.g. EcoRI=G^AATTC "
                        "(repeat for several; default: 5 built-in enzymes)")
    size = p.add_mutually_exclusive_group()
    size.add_argument("--short", action="store_true",
                      help=f"use only the first {SHORT_LEN} nt of each file (skips the question)")
    size.add_argument("--full", action="store_true",
                      help="digest the full genomes (skips the question)")
    p.add_argument("--png", default="gel2.png", help="picture of the gel (default: gel2.png)")
    p.add_argument("--no-png", action="store_true", help="do not save the picture")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    short = args.short or (not args.full and ask_short())
    limit = SHORT_LEN if short else None
    try:
        enzymes = [parse_enzyme(s) for s in (args.enzyme or DEFAULT_ENZYMES)]
        genomes = [load_genome("A", args.fasta1, limit), load_genome("B", args.fasta2, limit)]
    except (ValueError, OSError) as err:
        sys.exit(f"Error: {err}")
    report(genomes, enzymes, None if args.no_png else args.png)


if __name__ == "__main__":
    main()