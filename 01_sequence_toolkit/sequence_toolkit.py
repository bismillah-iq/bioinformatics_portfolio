"""
Sequence Analysis Toolkit
-------------------------
Given a GenBank accession (or a local FASTA/GenBank file), this script:

  1. Reports basic stats (length, base counts, GC content)
  2. Computes the reverse complement, transcript (RNA) and a translation
  3. Finds Open Reading Frames (ORFs) on all six reading frames
  4. Calculates codon usage for the longest ORF
  5. Saves a GC-content sliding-window plot and a codon-usage plot

Usage examples:
    python sequence_toolkit.py NM_000546            # fetch human TP53 mRNA from NCBI
    python sequence_toolkit.py my_gene.fasta        # use a local file
    python sequence_toolkit.py NM_000546 --min-aa 50 --window 100
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # save figures to files; no display window needed
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

# NCBI asks for a real email address so they can contact you if needed.
Entrez.email = "xeeeeva1030@gmail.com"


# --------------------------------------------------------------------------
# Loading data
# --------------------------------------------------------------------------
def load_record(source):
    """Load a sequence record from a local file or fetch it from NCBI."""
    path = Path(source)
    if path.exists():
        fmt = "genbank" if path.suffix.lower() in {".gb", ".gbk"} else "fasta"
        return SeqIO.read(path, fmt)

    print(f"Fetching {source} from NCBI...")
    with Entrez.efetch(db="nucleotide", id=source, rettype="gb", retmode="text") as handle:
        return SeqIO.read(handle, "genbank")


# --------------------------------------------------------------------------
# Core analysis functions
# --------------------------------------------------------------------------
def base_counts(seq):
    counts = Counter(str(seq).upper())
    return {base: counts.get(base, 0) for base in "ACGT"}


def gc_sliding_window(seq, window=100, step=10):
    """Return (positions, gc_percent) for a sliding window along the sequence."""
    positions, values = [], []
    for start in range(0, len(seq) - window + 1, step):
        positions.append(start + window // 2)
        values.append(gc_fraction(seq[start:start + window]) * 100)
    return positions, values


def find_orfs(seq, min_aa=30):
    """
    Find ORFs (ATG ... stop) on all six reading frames.
    Coordinates are 0-based positions on the strand being read
    (for '-' strand hits, positions refer to the reverse complement).
    Only ORFs that end with a stop codon are reported.
    """
    orfs = []
    for strand, s in (("+", seq), ("-", seq.reverse_complement())):
        for frame in range(3):
            sub = s[frame:]
            sub = sub[: len(sub) // 3 * 3]
            protein = str(sub.translate())
            chunks = protein.split("*")
            pos = 0  # amino-acid position of the current chunk
            for i, chunk in enumerate(chunks):
                has_stop = i < len(chunks) - 1
                start_m = chunk.find("M")
                if has_stop and start_m != -1 and len(chunk) - start_m >= min_aa:
                    aa_start = pos + start_m
                    nt_start = frame + aa_start * 3
                    nt_end = frame + (pos + len(chunk)) * 3  # end of last codon before stop
                    orfs.append({
                        "strand": strand,
                        "frame": frame + 1,
                        "nt_start": nt_start + 1,  # 1-based for readability
                        "nt_end": nt_end,
                        "length_aa": len(chunk) - start_m,
                        "protein": chunk[start_m:],
                        "cds": str(s[nt_start:nt_end]),
                    })
                pos += len(chunk) + 1  # +1 for the stop codon
    return sorted(orfs, key=lambda o: o["length_aa"], reverse=True)


def codon_usage(cds):
    """Count codons in a coding sequence (string, length multiple of 3)."""
    codons = [cds[i:i + 3] for i in range(0, len(cds) - len(cds) % 3, 3)]
    return Counter(codons)


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------
def plot_gc(positions, values, overall_gc, title, outfile):
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(positions, values, linewidth=1.5)
    ax.axhline(overall_gc, linestyle="--", color="gray", label=f"Overall GC: {overall_gc:.1f}%")
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("GC content (%)")
    ax.set_title(f"GC content sliding window - {title}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def plot_codon_usage(usage, title, outfile, top_n=30):
    items = usage.most_common(top_n)
    codons = [c for c, _ in items]
    counts = [n for _, n in items]
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(codons, counts)
    ax.set_ylabel("Count")
    ax.set_title(f"Top {len(items)} codons in longest ORF - {title}")
    plt.setp(ax.get_xticklabels(), rotation=90)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Simple DNA sequence analysis toolkit")
    parser.add_argument("source", help="GenBank accession (e.g. NM_000546) or path to a FASTA/GenBank file")
    parser.add_argument("--min-aa", type=int, default=30, help="Minimum ORF length in amino acids (default 30)")
    parser.add_argument("--window", type=int, default=100, help="GC sliding window size in bp (default 100)")
    parser.add_argument("--out", default="results", help="Output folder (default: results)")
    args = parser.parse_args()

    record = load_record(args.source)
    seq = record.seq.upper()
    out = Path(args.out)
    out.mkdir(exist_ok=True)

    # 1. Basic stats
    gc = gc_fraction(seq) * 100
    print(f"\nRecord:  {record.id} - {record.description}")
    print(f"Length:  {len(seq)} bp")
    print(f"Counts:  {base_counts(seq)}")
    print(f"GC:      {gc:.2f}%")

    # 2. Reverse complement, RNA, translation
    print(f"\nReverse complement (first 60 bp): {seq.reverse_complement()[:60]}")
    print(f"Transcript RNA (first 60 nt):      {seq.transcribe()[:60]}")
    trimmed = seq[: len(seq) // 3 * 3]
    print(f"Frame 1 translation (first 30 aa): {trimmed.translate()[:30]}")

    # 3. ORFs
    orfs = find_orfs(seq, min_aa=args.min_aa)
    print(f"\nFound {len(orfs)} ORFs of at least {args.min_aa} aa")
    with open(out / "orfs.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["strand", "frame", "nt_start", "nt_end", "length_aa", "protein"])
        for o in orfs:
            writer.writerow([o["strand"], o["frame"], o["nt_start"], o["nt_end"], o["length_aa"], o["protein"]])
    for o in orfs[:5]:
        print(f"  strand {o['strand']} frame {o['frame']}: {o['nt_start']}-{o['nt_end']} ({o['length_aa']} aa)")

    # 4. Plots
    window = min(args.window, len(seq))
    positions, values = gc_sliding_window(seq, window=window, step=max(1, window // 10))
    plot_gc(positions, values, gc, record.id, out / "gc_content.png")

    if orfs:
        usage = codon_usage(orfs[0]["cds"])
        plot_codon_usage(usage, record.id, out / "codon_usage.png")
        print(f"\nLongest ORF protein (first 60 aa): {orfs[0]['protein'][:60]}")

    print(f"\nSaved results to ./{out}/ (orfs.csv, gc_content.png, codon_usage.png)")


if __name__ == "__main__":
    main()