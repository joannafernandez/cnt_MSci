# 📈 pyBigWig Tutorial  
*A practical guide to reading, summarizing, and comparing bigWig files in Python*

This tutorial introduces `pyBigWig`, a Python library for reading bigWig signal tracks.

Focus:

- opening bigWigs
- inspecting chromosomes
- extracting signal
- binning
- plotting
- comparing tracks
- region-based quantification

This is especially useful for CUT&Tag / CUT&RUN workflows.

---

# 📦 Installation

    pip install pybigwig

or

    conda install -c bioconda pybigwig

As always, install directly in terminal, never in the IDE console

---

# 📥 Imports

    import pyBigWig
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

---

# 1️⃣ Opening a bigWig

    bw = pyBigWig.open("sample.bigWig")

Always close when finished:

    bw.close()

---

# 2️⃣ Inspect chromosomes

List chromosomes and lengths:

    bw.chroms()

Example output:

    {'chr1': 248956422, 'chr2': 242193529, ...}

Check if chr7 exists:

    "chr7" in bw.chroms()

Get chromosome length:

    bw.chroms()["chr7"]

---

# 3️⃣ Extract raw signal from a region

Values are returned per-base (or per internal block):

    vals = bw.values("chr7", 1000000, 1010000, numpy=True)

This returns a NumPy array.

NaNs indicate missing signal.

Drop NaNs:

    vals = vals[~np.isnan(vals)]

Mean signal in region:

    vals.mean()

---

# 4️⃣ Turn a region into a DataFrame

Example: extract chr7 in 50 bp bins.

    def bigwig_to_df(bigwig_path, chrom="chr7", bin_size=50):

        bw = pyBigWig.open(bigwig_path)
        chrom_len = bw.chroms()[chrom]

        starts = list(range(0, chrom_len, bin_size))
        ends = [min(s + bin_size, chrom_len) for s in starts]

        values = bw.values(chrom, 0, chrom_len, numpy=True)

        rows = []

        for s, e in zip(starts, ends):
            v = values[s:e]
            v = v[~np.isnan(v)]
            score = v.mean() if len(v) > 0 else 0
            rows.append([chrom, s + bin_size//2, score])

        bw.close()

        return pd.DataFrame(rows, columns=["Chr","Pos","score"])

Usage:

    df = bigwig_to_df("sample.bigWig")

---

# 5️⃣ Simple plotting

    plt.figure(figsize=(10,3))
    plt.plot(df["Pos"], df["score"], linewidth=0.5)
    plt.xlabel("Position")
    plt.ylabel("Signal")
    plt.tight_layout()

---

# 6️⃣ Smoothing signal

Rolling mean (e.g. 5 kb window if bin=50 bp):

    df["smooth"] = df["score"].rolling(window=100, center=True).mean()

Plot:

    plt.figure(figsize=(10,3))
    plt.plot(df["Pos"], df["smooth"])
    plt.tight_layout()

---

# 7️⃣ Comparing two bigWigs (local difference)

Load both:

    df1 = bigwig_to_df("WT.bigWig")
    df2 = bigwig_to_df("KO.bigWig")

Assuming identical bins:

    diff = df1.copy()
    diff["delta"] = df1["score"] - df2["score"]
    diff["smooth_delta"] = diff["delta"].rolling(window=100, center=True).mean()

Plot:

    plt.figure(figsize=(10,3))
    plt.plot(diff["Pos"], diff["smooth_delta"])
    plt.axhline(0, color="black", linewidth=0.5)
    plt.tight_layout()

Interpretation:

Positive = WT enriched  
Negative = KO enriched  

---

# 8️⃣ Region-based quantification

Define locus:

    locus_start = 120000000
    locus_end   = 121000000

Extract region mean:

    def region_mean(df, start, end):
        sub = df[(df["Pos"] >= start) & (df["Pos"] <= end)]
        return sub["score"].dropna().mean()

Example:

    wt_val = region_mean(df1, locus_start, locus_end)
    ko_val = region_mean(df2, locus_start, locus_end)

Difference:

    wt_val - ko_val

---

# 9️⃣ Replicates → statistics

Suppose:

    WT_R1.bigWig
    WT_R2.bigWig
    WT_R3.bigWig

Load:

    dfs = [
        bigwig_to_df("WT_R1.bigWig"),
        bigwig_to_df("WT_R2.bigWig"),
        bigwig_to_df("WT_R3.bigWig"),
    ]

Per-replicate region means:

    vals = [region_mean(d, locus_start, locus_end) for d in dfs]

Mean ± SD:

    np.mean(vals), np.std(vals)

---

# 🔟 Paired statistics (WT vs KO)

After computing:

    wt_vals = [...]
    ko_vals = [...]

Run paired test:

    from scipy.stats import ttest_rel, wilcoxon

    ttest_rel(wt_vals, ko_vals)
    wilcoxon(wt_vals, ko_vals)

---

# 1️⃣1️⃣ Mean ± SD tracks across replicates

Stack replicates:

    stack = pd.DataFrame({
        "Pos": dfs[0]["Pos"],
        "R1": dfs[0]["score"],
        "R2": dfs[1]["score"],
        "R3": dfs[2]["score"],
    })

    stack["mean"] = stack[["R1","R2","R3"]].mean(axis=1)
    stack["sd"]   = stack[["R1","R2","R3"]].std(axis=1)

Plot:

    plt.figure(figsize=(10,3))
    plt.plot(stack["Pos"], stack["mean"])
    plt.fill_between(stack["Pos"],
                     stack["mean"]-stack["sd"],
                     stack["mean"]+stack["sd"],
                     alpha=0.3)
    plt.tight_layout()

---

# 🧠 Important concepts

- bigWig = continuous signal, not peaks
- values represent fragment coverage (after normalization)
- comparisons assume identical binning
- spike-in affects scale but not workflow
- statistics must be done on regions, not single bins

---

# ✅ Typical CUT&Tag usage

- visualize loci
- compute local differences
- summarize signal over LADs / promoters
- generate region matrices
- validate peak changes

---

# 📚 Summary

pyBigWig lets you:

- read bigWigs
- extract signal
- bin genome
- smooth tracks
- compare conditions
- quantify regions
- perform statistics

It pairs naturally with PyRanges (regions) + seaborn/matplotlib (plots).

---

Happy signal wrangling 📈
