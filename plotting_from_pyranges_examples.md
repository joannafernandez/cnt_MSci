# 📊 Plotting peak distributions (PyRanges → seaborn)
*Copy-paste ready examples for CUT&Tag / CUT&RUN*

This mini-tutorial shows how to go from a `PyRanges` object → a tidy pandas DataFrame → common seaborn plots:
- peak length distributions
- peaks per chromosome
- distance to nearest TSS
- promoter vs gene body vs IGR fractions
- replicate / condition comparisons

---

## ✅ Setup

    import pyranges as pr
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

Recommended defaults:

    sns.set_context("talk")

Load peaks:

    peaks = pr.read_bed("peaks.bed").sort()

Optional: if you have multiple conditions, keep a label:

    peaks = peaks.assign("condition", "WT_siSCR")

---

## 1) Peak length distribution

Compute lengths:

    df = peaks.df.copy()
    df["length"] = df["End"] - df["Start"]

Histogram:

    plt.figure()
    sns.histplot(data=df, x="length", bins=60)
    plt.xlabel("Peak length (bp)")
    plt.ylabel("Count")
    plt.title("Peak length distribution")
    plt.tight_layout()

Log-scaled x (often nicer for peak lengths):

    plt.figure()
    sns.histplot(data=df, x="length", bins=60)
    plt.xscale("log")
    plt.xlabel("Peak length (bp, log scale)")
    plt.ylabel("Count")
    plt.title("Peak length distribution (log x)")
    plt.tight_layout()

Violin (good for comparing conditions):

    plt.figure()
    sns.violinplot(data=df, x="condition", y="length", cut=0)
    plt.yscale("log")
    plt.ylabel("Peak length (bp, log scale)")
    plt.xlabel("")
    plt.title("Peak length by condition")
    plt.tight_layout()

---

## 2) Peaks per chromosome (bar plot)

Counts:

    chr_counts = (
        peaks.df.groupby("Chromosome")
        .size()
        .reset_index(name="n_peaks")
        .sort_values("n_peaks", ascending=False)
    )

Plot:

    plt.figure(figsize=(10,4))
    sns.barplot(data=chr_counts, x="Chromosome", y="n_peaks")
    plt.xticks(rotation=90)
    plt.ylabel("Number of peaks")
    plt.title("Peaks per chromosome")
    plt.tight_layout()

Tip: if you want autosomes only:

    autosomes = [f"chr{i}" for i in range(1, 23)]
    peaks_auto = peaks[peaks.Chromosome.isin(autosomes)]

---

## 3) Peak width vs signal (if signal exists)

If your BED has a signal column (often column 4). Adapt index if needed:

    df = peaks.df.copy()
    df["length"] = df["End"] - df["Start"]
    df["signal"] = pd.to_numeric(df.iloc[:, 3], errors="coerce")

Scatter (downsample if huge):

    df_small = df.sample(n=min(len(df), 20000), random_state=1)

    plt.figure()
    sns.scatterplot(data=df_small, x="length", y="signal", alpha=0.3, linewidth=0)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Peak length (bp)")
    plt.ylabel("Signal")
    plt.title("Peak length vs signal")
    plt.tight_layout()

---

## 4) Distance to nearest TSS (distribution)

### Build TSS PyRanges from GTF

Load genes (GTF) and keep genes:

    genes = pr.read_gtf("genes.gtf")
    genes = genes[genes.Feature == "gene"].sort()

Create a 1bp TSS interval per gene.
This keeps strand-aware TSS:

    gdf = genes.df.copy()

    tss_start = np.where(gdf["Strand"] == "+", gdf["Start"], gdf["End"] - 1)
    tss_end   = tss_start + 1

    tss_df = pd.DataFrame({
        "Chromosome": gdf["Chromosome"].values,
        "Start": tss_start,
        "End": tss_end,
        "gene_id": gdf.get("gene_id", pd.Series([None]*len(gdf))).values,
        "Strand": gdf["Strand"].values,
    })

    tss = pr.PyRanges(tss_df).sort()

Compute nearest TSS:

    nearest = peaks.nearest(tss)

Get distances:

    nd = nearest.df.copy()
    nd["distance_to_tss"] = nd["Distance"].astype(int)

Plot:

    plt.figure()
    sns.histplot(data=nd, x="distance_to_tss", bins=80)
    plt.xlabel("Distance to nearest TSS (bp)")
    plt.ylabel("Count")
    plt.title("Distance to nearest TSS")
    plt.tight_layout()

Often you only care about within ±100 kb:

    plt.figure()
    sns.histplot(data=nd[nd["distance_to_tss"].abs() <= 100_000],
                 x="distance_to_tss", bins=80)
    plt.xlabel("Distance to nearest TSS (bp) within ±100 kb")
    plt.ylabel("Count")
    plt.title("Distance to nearest TSS (zoomed)")
    plt.tight_layout()
