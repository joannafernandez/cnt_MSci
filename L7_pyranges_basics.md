# 🧬 PyRanges Basics Tutorial  
*A practical guide to genomic interval operations in Python*

PyRanges is a fast, pandas-friendly library for working with genomic intervals (BED/GTF-style data).  
It is built on top of pandas + numpy and provides efficient overlap, join, clustering, and annotation operations.

This tutorial covers the **most important core functions** you’ll use in everyday genomics workflows.

---

## 📦 Installation

    pip install pyranges

or

    conda install -c bioconda pyranges

---

## 📥 Import

    import pyranges as pr
    import pandas as pd

---

# 1️⃣ Creating a PyRanges object

PyRanges objects are created from pandas DataFrames.

Minimum required columns:

- Chromosome  
- Start  
- End  

Example:

    df = pd.DataFrame({
        "Chromosome": ["chr1", "chr1", "chr2"],
        "Start": [100, 500, 200],
        "End": [200, 600, 300],
        "name": ["peak1", "peak2", "peak3"]
    })

    gr = pr.PyRanges(df)
    gr

---

## From BED file

    gr = pr.read_bed("peaks.bed")

From GTF:

    gr = pr.read_gtf("genes.gtf")

---

# 2️⃣ Basic inspection

    gr.df
    len(gr)
    gr.columns
    gr.head()

---

# 3️⃣ Sorting

    gr_sorted = gr.sort()

Always sort before downstream operations.

---

# 4️⃣ Subsetting

By chromosome:

    gr_chr1 = gr[gr.Chromosome == "chr1"]

By length:

    gr_long = gr[(gr.End - gr.Start) > 500]

---

# 5️⃣ Interval length

    gr.lengths()

Add as column:

    gr = gr.assign("length", lambda df: df.End - df.Start)

---

# 6️⃣ Overlaps

Find overlapping intervals:

    overlaps = gr1.overlap(gr2)

Count overlaps:

    gr1.count_overlaps(gr2)

Adds column:

    NumberOverlaps

---

# 7️⃣ Join (annotation)

Equivalent to bedtools intersect.

    joined = peaks.join(genes)

---

# 8️⃣ Nearest feature

    nearest = peaks.nearest(genes)

Adds column:

- Distance

---

# 9️⃣ Subtract

    clean = peaks.subtract(blacklist)

---

# 🔟 Merge overlapping intervals

    merged = peaks.merge()

Preserve metadata:

    merged = peaks.merge(by="sample")

---

# 1️⃣1️⃣ Cluster

    clustered = peaks.cluster()

Adds column:

    Cluster

---

# 1️⃣2️⃣ Coverage

    coverage = peaks.coverage()

Against another PyRanges:

    coverage = peaks.coverage(other)

---

# 1️⃣3️⃣ Concatenate PyRanges

    combined = pr.concat([gr1, gr2, gr3])

---

# 1️⃣4️⃣ Drop duplicate positions

    unique = gr.drop_duplicate_positions()

---

# 1️⃣5️⃣ Export

BED:

    gr.to_bed("output.bed")

CSV:

    gr.df.to_csv("output.csv", index=False)

---

# 1️⃣6️⃣ Convert back to pandas

    df = gr.df

---

# 🧠 Common genomics patterns

Peaks in genes:

    peaks_in_genes = peaks.join(genes)

Peaks in IGR:

    peaks_in_igr = peaks.subtract(genes)

Shared peaks:

    shared = sample1.overlap(sample2)

Shared across multiple samples:

    shared = pr.concat([s1, s2, s3]).cluster()
    shared = shared.df.groupby("Cluster").filter(lambda x: len(x) == 3)

---

# 🚀 Performance tips

- Always sort
- Prefer PyRanges ops over pandas loops
- Use .df only at the end
- Avoid row-wise iteration

---

# 📚 References

- https://pyranges.readthedocs.io  
- https://github.com/pyranges/pyranges  

---

# ✅ Summary

Task | Function  
Overlap | overlap()  
Annotate | join()  
Nearest | nearest()  
Remove | subtract()  
Merge | merge()  
Cluster | cluster()  
Combine | pr.concat()  

---

Happy genome wrangling 🧬
