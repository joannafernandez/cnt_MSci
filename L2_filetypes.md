# NGS Data types and formats
> a key resource with great examples can be found here: https://ngs101.com/the-complete-guide-to-ngs-data-types-and-formats-from-raw-reads-to-analysis-ready-files/

LO1: Understand the generic structure of an NGS workflow

LO2: Ifentify and recognise the various files and data types


---
## The experiment begins with the cellular and molecular biology
<img width="1477" height="409" alt="image" src="https://github.com/user-attachments/assets/a585bce8-421c-4783-ba0f-5310344f4c49" />

## The output from the sequencing experiment is where the bioinformatics pipeline begins
Your previous experience was using .csv files as your input for your analysis 

<img width="318.5" height="272" alt="image" src="https://github.com/user-attachments/assets/bcc1a02a-f083-4743-92c7-8a6d8c23f4d4" />

with this project, you will learn about the other "black box" programs that take raw sequencing data (FastQs) as an input, and how we process this in large workflows 
<img width="1290" height="556" alt="image" src="https://github.com/user-attachments/assets/2f60b907-6432-4876-a751-950dd926ca05" />

---
## Common Genomic File Formats

This section summarises the most common file formats encountered in sequencing-based genomic analyses, their purpose, and how they are typically used.

---

### FASTQ (`.fastq`, `.fq`, `.fastq.gz`)

**Purpose:**  
Stores raw sequencing reads and their associated quality scores.

**Key features:**
- Four lines per read:
  1. Read identifier
  2. Nucleotide sequence
  3. Separator (`+`)
  4. ASCII-encoded quality scores
- Generated directly by the sequencing instrument
- Often compressed (`.gz`)

**Typical use:**
- Input for alignment
- Quality control (read quality, adapter content)

---

### SAM (`.sam`)

**Purpose:**  
Text-based format storing aligned sequencing reads.

**Key features:**
- Human-readable
- Contains alignment position, mapping quality, CIGAR string, flags
- Large file size

**Typical use:**
- Debugging alignments
- Intermediate format before conversion to BAM

---

### BAM (`.bam`)

**Purpose:**  
Binary, compressed version of SAM.

**Key features:**
- Much smaller than SAM
- Efficient for storage and computation
- Requires an index file (`.bai`) for random access

**Typical use:**
- Standard alignment file for downstream analysis
- Input for peak calling, coverage calculation, visualisation

---

### BAI (`.bai`)

**Purpose:**  
Index file for BAM files.

**Key features:**
- Enables fast random access to specific genomic regions
- Required for genome browsers and many analysis tools

**Typical use:**
- Paired with BAM files
- Allows tools to query alignments without reading the entire file

---

### BED (`.bed`)

**Purpose:**  
Defines genomic intervals (regions).

**Key features:**
- Plain text, tab-delimited
- Minimum 3 columns:
  - Chromosome
  - Start (0-based)
  - End (1-based)
- Can include additional annotation columns

**Typical use:**
- Peak regions
- Genomic annotations
- Input for motif analysis or region-based counting

---

### BedGraph (`.bedgraph`)

**Purpose:**  
Stores continuous signal values across genomic intervals.

**Key features:**
- Similar to BED but includes a signal/value column
- Represents coverage or enrichment
- Not indexed by default

**Typical use:**
- Intermediate coverage representation
- Input for BigWig conversion

---

### BigWig (`.bw`, `.bigWig`)

**Purpose:**  
Efficient, indexed binary format for continuous genomic signal.

**Key features:**
- Compressed and indexed
- Supports fast region-based access
- Ideal for genome browser visualisation

**Typical use:**
- Visualisation of coverage or enrichment tracks
- Metaplots and heatmaps
- Sharing processed signal data

---

### GFF / GTF (`.gff`, `.gff3`, `.gtf`)

**Purpose:**  
Annotation format describing genomic features.

**Key features:**
- Structured, tab-delimited format
- Describes genes, transcripts, exons, CDS, etc.
- Includes feature attributes as key–value pairs

**Typical use:**
- Genome annotation
- Assigning reads or peaks to genes
- Required input for counting tools (e.g. featureCounts)

---

## Summary Table

| Format | Type | Stores | Indexed | Typical Role |
|------|-----|--------|--------|-------------|
| FASTQ | Text | Raw reads + quality | No | Sequencer output |
| SAM | Text | Alignments | No | Debugging |
| BAM | Binary | Alignments | Yes (BAI) | Core analysis file |
| BAI | Binary | BAM index | N/A | Random access |
| BED | Text | Genomic regions | No | Peaks / intervals |
| BedGraph | Text | Continuous signal | No | Coverage data |
| BigWig | Binary | Continuous signal | Yes | Visualisation |
| GFF/GTF | Text | Genome annotation | No | Feature annotation |

---

