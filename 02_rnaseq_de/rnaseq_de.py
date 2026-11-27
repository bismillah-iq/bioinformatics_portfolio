"""
RNA-seq Differential Expression Pipeline
----------------------------------------
Takes a raw count matrix (genes x samples) and a sample metadata table and runs:

  1. Low-count gene filtering
  2. QC: library sizes
  3. Normalization + PCA (sample-level overview)
  4. Differential expression with PyDESeq2 (Python implementation of DESeq2)
  5. Volcano plot, MA plot, and a heatmap of the top DE genes
  6. Saves the full results table (CSV)

Usage:
    python rnaseq_de.py --demo
    python rnaseq_de.py --counts counts.csv --meta meta.csv \
        --condition-col condition --treated treated --control untreated
    # paired design (e.g. tumor and normal from the same patient) + gene names:
    python rnaseq_de.py --counts counts.tsv.gz --meta meta.csv \
        --treated tumor --control normal --pair-col patient --annotate

Input formats:
    counts.csv : first column = gene IDs/names, one column per sample, raw integer counts
    meta.csv   : first column = sample names (must match counts columns), plus a
                 column holding the condition of each sample
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

PALETTE = {"up": "#d1495b", "down": "#2e86ab", "ns": "#b8b8b8"}


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------
def make_demo_data(n_genes=6000, n_per_group=4, frac_de=0.10, seed=7):
    """
    Simulated RNA-seq counts (negative binomial) with a known set of truly
    differentially expressed genes. For testing the pipeline only - NOT real data.
    """
    rng = np.random.default_rng(seed)
    base_mean = rng.lognormal(mean=5, sigma=1.5, size=n_genes)
    de_idx = rng.choice(n_genes, int(n_genes * frac_de), replace=False)
    lfc = np.zeros(n_genes)
    lfc[de_idx] = rng.choice([-1, 1], len(de_idx)) * rng.uniform(0.8, 3.0, len(de_idx))

    samples, conds, cols = [], [], {}
    lib_factors = rng.uniform(0.7, 1.3, 2 * n_per_group)
    dispersion = 0.1 + 1.0 / np.sqrt(base_mean)
    for i in range(2 * n_per_group):
        cond = "control" if i < n_per_group else "treated"
        mu = base_mean * lib_factors[i] * (2 ** lfc if cond == "treated" else 1)
        r = 1.0 / dispersion
        counts = rng.negative_binomial(r, r / (r + mu))
        name = f"{cond}_{i % n_per_group + 1}"
        cols[name] = counts
        samples.append(name)
        conds.append(cond)

    counts_df = pd.DataFrame(cols, index=[f"Gene{i + 1:05d}" for i in range(n_genes)])
    meta_df = pd.DataFrame({"condition": conds}, index=samples)
    truth = pd.Series(lfc != 0, index=counts_df.index)
    return counts_df, meta_df, truth


def load_data(counts_path, meta_path):
    def _sep(path):
        name = str(path).lower()
        return "\t" if (".tsv" in name or name.endswith((".txt", ".txt.gz"))) else ","

    counts = pd.read_csv(counts_path, sep=_sep(counts_path), index_col=0)  # .gz is read directly
    meta = pd.read_csv(meta_path, sep=_sep(meta_path), index_col=0)

    # keep only numeric count columns that are listed in the metadata
    counts = counts.select_dtypes(include=[np.number])
    shared = [s for s in meta.index if s in counts.columns]
    if len(shared) < 4:
        raise SystemExit(
            "Fewer than 4 sample names match between counts columns and metadata rows. "
            "Check that the metadata's first column holds the sample names."
        )
    counts = counts[shared].round().astype(int)
    meta = meta.loc[shared]
    return counts, meta



# --------------------------------------------------------------------------
# Gene annotation (optional)
# --------------------------------------------------------------------------
def annotate_symbols(gene_ids):
    """
    Map Ensembl gene IDs to gene symbols using the MyGene.info web service.
    Needs an internet connection. Returns a dict {ensembl_id: symbol}.
    """
    import mygene

    mg = mygene.MyGeneInfo()
    hits = mg.querymany(list(gene_ids), scopes="ensembl.gene", fields="symbol",
                        species="human", verbose=False)
    mapping = {}
    for h in hits:
        if h.get("symbol") and h["query"] not in mapping:
            mapping[h["query"]] = h["symbol"]
    return mapping


def add_labels(res, mapping=None):
    """Add a 'label' column: gene symbol when known, otherwise the original ID."""
    res = res.copy()
    if mapping:
        res["symbol"] = [mapping.get(g, "") for g in res.index]
        labels = [sym if sym else g for g, sym in zip(res.index, res["symbol"])]
    else:
        labels = list(res.index)
    # make labels unique (two IDs can share a symbol)
    seen = {}
    unique = []
    for lab, gid in zip(labels, res.index):
        seen[lab] = seen.get(lab, 0) + 1
        unique.append(lab if seen[lab] == 1 else f"{lab} ({gid})")
    res["label"] = unique
    return res

# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------
def plot_library_sizes(counts, meta, cond_col, outfile):
    sizes = counts.sum(axis=0) / 1e6
    conds = meta[cond_col].astype("category")
    colors = plt.cm.Set2(conds.cat.codes)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(sizes.index, sizes.values, color=colors)
    ax.set_ylabel("Total counts (millions)")
    ax.set_title("Library sizes per sample")
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right")
    handles = [plt.Rectangle((0, 0), 1, 1, color=plt.cm.Set2(i)) for i in range(len(conds.cat.categories))]
    ax.legend(handles, conds.cat.categories)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def plot_pca(log_norm, meta, cond_col, outfile):
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(log_norm.T.values)
    var = pca.explained_variance_ratio_ * 100
    conds = meta[cond_col].astype("category")
    fig, ax = plt.subplots(figsize=(6, 5))
    for i, c in enumerate(conds.cat.categories):
        mask = (conds == c).values
        ax.scatter(pcs[mask, 0], pcs[mask, 1], s=90, color=plt.cm.Set2(i), label=c, edgecolor="k")
    for name, (x, y) in zip(log_norm.columns, pcs):
        ax.annotate(name, (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% variance)")
    ax.set_title("PCA of normalized expression")
    ax.legend(title=cond_col)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def plot_volcano(res, lfc_cut, padj_cut, label_top, outfile, highlight=None):
    d = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    d["neglogp"] = -np.log10(d["padj"].clip(lower=1e-300))
    d["status"] = "ns"
    d.loc[(d["padj"] < padj_cut) & (d["log2FoldChange"] >= lfc_cut), "status"] = "up"
    d.loc[(d["padj"] < padj_cut) & (d["log2FoldChange"] <= -lfc_cut), "status"] = "down"

    fig, ax = plt.subplots(figsize=(7, 6))
    for status, color in PALETTE.items():
        s = d[d["status"] == status]
        ax.scatter(s["log2FoldChange"], s["neglogp"], s=10, alpha=0.7, color=color,
                   label=f"{status} ({len(s)})", linewidths=0)
    ax.axvline(lfc_cut, ls="--", color="gray", lw=0.8)
    ax.axvline(-lfc_cut, ls="--", color="gray", lw=0.8)
    ax.axhline(-np.log10(padj_cut), ls="--", color="gray", lw=0.8)

    if highlight:
        wanted = {g.strip().upper() for g in highlight}
        names = d["symbol"] if "symbol" in d else d.get("label", pd.Series(d.index, index=d.index))
        top = d[names.astype(str).str.upper().isin(wanted)]
        missing = wanted - set(names.astype(str).str.upper())
        if missing:
            print(f"  highlight genes not found: {sorted(missing)}")
    else:
        top = d[d["status"] != "ns"].sort_values("padj").head(label_top)
    texts = []
    for _, row in top.iterrows():
        color = PALETTE["up"] if row["log2FoldChange"] > 0 else PALETTE["down"]
        texts.append(ax.text(row["log2FoldChange"], row["neglogp"], row.get("label", row.name),
                             fontsize=8, fontweight="bold", color=color))
    try:
        from adjustText import adjust_text

        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
    except ImportError:
        pass  # labels may overlap a little without adjustText

    # keep the view readable if a few low-count genes have extreme fold changes
    lim = 9
    outside = int((d["log2FoldChange"].abs() > lim).sum())
    ax.set_xlim(-lim, lim)
    if outside:
        ax.text(0.99, 0.01, f"{outside} genes with |log2FC| > {lim} not shown", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7, color="gray")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Volcano plot")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)
    return d["status"].value_counts().to_dict()


def plot_ma(res, padj_cut, outfile):
    d = res.dropna(subset=["padj", "log2FoldChange", "baseMean"])
    sig = d["padj"] < padj_cut
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(d.loc[~sig, "baseMean"], d.loc[~sig, "log2FoldChange"], s=8, color=PALETTE["ns"], linewidths=0)
    ax.scatter(d.loc[sig, "baseMean"], d.loc[sig, "log2FoldChange"], s=10, color=PALETTE["up"], linewidths=0,
               label=f"padj < {padj_cut}")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Mean of normalized counts")
    ax.set_ylabel("log2 fold change")
    ax.set_title("MA plot")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def plot_heatmap(log_norm, res, meta, cond_col, padj_cut, n_genes, outfile):
    import seaborn as sns

    sig = res.dropna(subset=["padj"]).query("padj < @padj_cut").sort_values("padj").head(n_genes)
    if sig.empty:
        print("No significant genes - skipping heatmap.")
        return
    mat = log_norm.loc[sig.index].copy()
    mat.index = sig["label"] if "label" in sig else sig.index
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)
    order = meta.sort_values(cond_col).index
    z = z[order]
    groups = list(meta[cond_col].astype(str).unique())
    lut = {c: plt.cm.Set2(i) for i, c in enumerate(sorted(groups))}
    col_colors = [lut[str(meta.loc[s, cond_col])] for s in order]
    g = sns.clustermap(z, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, col_cluster=False,
                       col_colors=col_colors, figsize=(7, 9),
                       yticklabels=len(z) <= 40, cbar_kws={"label": "z-score (log expression)"})
    g.fig.suptitle(f"Top {len(z)} differentially expressed genes", y=1.01)
    g.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(g.fig)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="RNA-seq differential expression with PyDESeq2")
    p.add_argument("--counts", help="Raw counts table (genes x samples), CSV/TSV")
    p.add_argument("--meta", help="Sample metadata table, CSV/TSV")
    p.add_argument("--condition-col", default="condition", help="Metadata column with the groups")
    p.add_argument("--treated", default="treated", help="Group to compare (numerator)")
    p.add_argument("--control", default="control", help="Reference group (denominator)")
    p.add_argument("--pair-col", default=None,
                   help="Metadata column identifying matched samples (e.g. patient) for a paired design")
    p.add_argument("--annotate", action="store_true",
                   help="Map Ensembl IDs to gene symbols via MyGene.info (needs internet)")
    p.add_argument("--highlight", default=None,
                   help="Comma-separated gene symbols to label on the volcano plot (e.g. SFTPC,COL11A1)")
    p.add_argument("--min-count", type=int, default=10, help="Keep genes with at least this many total counts")
    p.add_argument("--padj", type=float, default=0.05, help="Adjusted p-value cutoff")
    p.add_argument("--lfc", type=float, default=1.0, help="|log2 fold change| cutoff for the volcano plot")
    p.add_argument("--demo", action="store_true", help="Run on simulated data to test the pipeline")
    p.add_argument("--out", default="results", help="Output folder")
    a = p.parse_args()

    out = Path(a.out)
    out.mkdir(exist_ok=True)
    truth = None

    if a.demo:
        print("Running on SIMULATED data (demo mode).")
        counts, meta, truth = make_demo_data()
        a.condition_col, a.treated, a.control = "condition", "treated", "control"
    else:
        if not (a.counts and a.meta):
            raise SystemExit("Provide --counts and --meta, or use --demo.")
        counts, meta = load_data(a.counts, a.meta)

    cond = a.condition_col
    for g in (a.treated, a.control):
        if g not in set(meta[cond]):
            raise SystemExit(f"Group '{g}' not found in column '{cond}'. Found: {sorted(set(meta[cond]))}")
    keep = meta[cond].isin([a.treated, a.control])
    meta = meta[keep]
    counts = counts[meta.index]
    print(f"Samples: {meta[cond].value_counts().to_dict()}")
    print(f"Genes before filtering: {counts.shape[0]}")

    # 1. filter low-count genes
    counts = counts[counts.sum(axis=1) >= a.min_count]
    print(f"Genes after filtering (total count >= {a.min_count}): {counts.shape[0]}")
    if truth is not None:
        truth = truth.loc[counts.index]

    # 2. QC
    plot_library_sizes(counts, meta, cond, out / "01_library_sizes.png")

    # 3. DESeq2 model
    if a.pair_col:
        if a.pair_col not in meta.columns:
            raise SystemExit(f"--pair-col '{a.pair_col}' not found in metadata columns: {list(meta.columns)}")
        design_cols = [a.pair_col, cond]
        design = f"~{a.pair_col} + {cond}"
        print(f"Paired design: {design}")
    else:
        design_cols = [cond]
        design = f"~{cond}"
    dds = DeseqDataSet(counts=counts.T, metadata=meta[design_cols].astype(str), design=design, quiet=True)
    dds.deseq2()
    normed = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names).T
    log_norm = np.log2(normed + 1)

    plot_pca(log_norm, meta, cond, out / "02_pca.png")

    # 4. DE test
    stats = DeseqStats(dds, contrast=[cond, a.treated, a.control], alpha=a.padj, quiet=True)
    stats.summary()
    res = stats.results_df.sort_values("padj")
    mapping = None
    if a.annotate:
        try:
            print("Looking up gene symbols (MyGene.info)...")
            mapping = annotate_symbols(res.index)
            print(f"  mapped {len(mapping)} of {len(res)} genes")
        except Exception as e:  # no internet, service down, package missing...
            print(f"  gene annotation failed ({type(e).__name__}: {e}); using original IDs")
    res = add_labels(res, mapping)
    res.to_csv(out / "de_results.csv")

    # 5. plots
    highlight = a.highlight.split(",") if a.highlight else None
    status = plot_volcano(res, a.lfc, a.padj, 10, out / "03_volcano.png", highlight)
    plot_ma(res, a.padj, out / "04_ma_plot.png")
    plot_heatmap(log_norm, res, meta, cond, a.padj, 40, out / "05_heatmap_top_genes.png")

    n_sig = int((res["padj"] < a.padj).sum())
    print(f"\n{a.treated} vs {a.control}")
    print(f"Significant genes (padj < {a.padj}): {n_sig}")
    print(f"  Up   (padj<{a.padj} & log2FC >= {a.lfc}):  {status.get('up', 0)}")
    print(f"  Down (padj<{a.padj} & log2FC <= -{a.lfc}): {status.get('down', 0)}")
    print("\nTop 10 genes:")
    print(res[["baseMean", "log2FoldChange", "padj"]].set_index(res["label"]).head(10).round(4).to_string())

    if truth is not None:
        called = (res["padj"] < a.padj).reindex(truth.index).fillna(False)
        tp = int((called & truth).sum())
        print(f"\n[demo check] {tp} of {int(truth.sum())} truly changed genes recovered; "
              f"{int((called & ~truth).sum())} false positives")

    print(f"\nSaved results to ./{out}/")


if __name__ == "__main__":
    main()
