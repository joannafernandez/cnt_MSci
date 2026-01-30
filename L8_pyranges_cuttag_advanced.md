# 🧬 Advanced PyRanges for CUT&Tag / CUT&RUN  
*Real-world genomics workflows + teaching exercises*

This tutorial assumes you already know basic PyRanges operations.

Focus:

- Peak annotation (promoter / gene body / IGR)
- Shared peaks across replicates
- Consensus peak sets
- Signal-aware filtering
- Preparing inputs for motif analysis
- Typical CUT&Tag logic

---

# 📥 Setup

    import pyranges as pr
    import pandas as pd

Load datasets:

    peaks = pr.read_bed("peaks.bed")
    genes = pr.read_gtf("genes.gtf")

Extract gene bodies:

    genes = genes[genes.Feature == "gene"]

Sort everything:

    peaks = peaks.sort()
    genes = genes.sort()

---

# 1️⃣ Promoters, gene bodies, IGR

Define promoters (±2 kb from TSS):

    promoters = genes.slack(2000, 0)

Rename for clarity:

    promoters = promoters.assign("region", "promoter")
    genes = genes.assign("region", "gene_body")

---

### Annotate peaks

Peaks in promoters:

    peaks_prom = peaks.join(promoters)

Peaks in gene bodies:

    peaks_gene = peaks.join(genes)

Peaks in IGR:

    peaks_igr = peaks.subtract(pr.concat([promoters, genes]))

---

### Combine annotations

    prom_df = peaks_prom.df.assign(category="promoter")
    gene_df = peaks_gene.df.assign(category="gene_body")
    igr_df  = peaks_igr.df.assign(category="igr")

    annotated = pd.concat([prom_df, gene_df, igr_df])

---

# 2️⃣ Peaks overlapping multiple features

Some peaks span promoter → gene body.

Detect these:

    both = peaks.join(promoters).join(genes)

Label:

    both_df = both.df.assign(category="promoter+gene")

---

# 3️⃣ Shared peaks across replicates

Load replicates:

    r1 = pr.read_bed("rep1.bed")
    r2 = pr.read_bed("rep2.bed")
    r3 = pr.read_bed("rep3.bed")

Sort:

    r1, r2, r3 = r1.sort(), r2.sort(), r3.sort()

Concatenate + cluster:

    all_peaks = pr.concat([r1, r2, r3]).cluster()

Require presence in all 3:

    consensus = all_peaks.df.groupby("Cluster").filter(lambda x: len(x) == 3)

Back to PyRanges:

    consensus = pr.PyRanges(consensus)

---

# 4️⃣ Signal-aware filtering

If your BED has signal column (e.g. column 4):

    peaks = peaks.assign("signal", peaks.df.iloc[:,3].astype(float))

Filter:

    strong = peaks[peaks.signal > 5]

---

# 5️⃣ Summit expansion for motif analysis

Expand summits to 151 bp:

    summits = pr.read_bed("summits.bed")

    summits_151 = summits.slack(75, 76)

Export:

    summits_151.to_bed("summits_151bp.bed")

---

# 6️⃣ Universe peak set across conditions

Combine conditions:

    all_conditions = pr.concat([WT, KO, siSCR, siFM]).cluster()

Create universe:

    universe = all_conditions.merge()

---

# 7️⃣ Boolean presence matrix

Example:

    universe = universe.join(WT, suffix="_WT")
    universe = universe.join(KO, suffix="_KO")

Presence table:

    df = universe.df
    df["WT_present"] = ~df["Start_WT"].isna()
    df["KO_present"] = ~df["Start_KO"].isna()

---

# 8️⃣ Save for HOMER / MEME

    universe.to_bed("universe_peaks.bed")

---
