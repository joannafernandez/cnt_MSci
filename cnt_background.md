## ✂️🏷️ 1. What is CUT&Tag?

**CUT&Tag (Cleavage Under Targets and Tagmentation)** is a targeted chromatin profiling technique used to map protein–DNA interactions and histone modifications genome-wide.  
It combines antibody-directed targeting with **Tn5 transposase–mediated tagmentation**, enabling high-resolution profiling using very low cell numbers.

**Key concept:**  
Instead of fragmenting chromatin globally, CUT&Tag tethers Tn5 directly to a protein of interest and cuts only where that protein is bound.

---

## 🎯 2. Core Principles

- Performed on intact, permeabilised cells or nuclei  
- Antibody binds a chromatin-associated protein or histone modification  
- Protein A–Tn5 fusion binds the antibody  
- Tn5 inserts sequencing adapters in situ at binding sites  
- DNA fragments are PCR-amplified and sequenced  

---

## 🧪🥼 3. Wet Lab Workflow

#### 3.1 Cell Preparation
- Fresh or lightly fixed cells (often unfixed)
- Cells immobilised on Concanavalin A–coated beads
- Mild permeabilisation preserves chromatin structure while allowing antibody access

---

#### 3.2 Primary Antibody Binding
- Incubation with antibody against:
  - Histone modifications (e.g. H3K27me3, H3K4me3)
  - Transcription factors
- Antibody specificity is critical for signal quality

---

#### 3.3 Secondary Antibody Binding (Optional)
- Increases signal strength
- Commonly used for:
  - Weak epitopes
  - Transcription factors

---

#### 3.4 pA-Tn5 Binding
- Protein A–Tn5 transposase pre-loaded with sequencing adapters
- Protein A binds the Fc region of the antibody
- Tn5 is tethered to the chromatin target

---

#### 3.5 Tagmentation
- Addition of Mg²⁺ activates Tn5
- Simultaneous DNA cleavage and adapter insertion
- Occurs only at antibody-bound sites

---

#### 3.6 DNA Release and Cleanup
- DNA released by heat and/or Proteinase K
- Minimal background DNA compared to ChIP-seq
- Cleanup using SPRI beads

---

#### 3.7 Library Amplification
- Indexed PCR amplification
- Fewer cycles required due to high signal-to-noise
- Libraries are sequencing-ready without additional ligation steps

---

## 4. 🖥️ Dry Lab Analysis Workflow

#### 4.1 Quality Control
- Adapter trimming (usually minimal)
- Read quality assessment (FastQC, MultiQC)
- Fragment size distribution inspection

---

#### 4.2 Alignment
- Alignment to reference genome (e.g. hg38, mm10)
- Common aligners:
  - Bowtie2
  - BWA
- High alignment rates expected

---

#### 4.3 Filtering
- Removal of:
  - Low mapping quality reads
  - Mitochondrial reads
- Duplicate removal often unnecessary due to low background

---

#### 4.4 Normalisation
Common approaches:
- Reads per million (RPM)
- Spike-in normalisation (optional)
- Total read count normalisation

---

#### 4.5 Peak Calling
- Narrow marks / transcription factors:
  - SEACR
  - MACS2 (with tuned parameters)
- Broad marks:
  - SEACR (stringent mode)
  - Bin-based enrichment approaches

---

#### 4.6 Visualisation and Downstream Analysis
- BigWig generation for genome browsers
- Heatmaps and metaplots
- Differential binding analysis
- Integration with RNA-seq or other epigenomic datasets

---

## 5. 💪🌽 Strengths and Limitations

#### Strengths
- Very low cell input requirement
- High signal-to-noise
- No sonication
- Rapid workflow
- Minimal background

#### Limitations
- Strong dependence on antibody quality
- Less historical benchmarking than ChIP-seq
- Some transcription factors remain challenging

---

## 6. 📋 Comparison with CUT&RUN and ChIP-seq

| Feature | CUT&Tag | CUT&RUN | ChIP-seq |
|------|--------|--------|---------|
| Fragmentation method | Tn5 tagmentation | MNase cleavage | Sonication |
| Targeting mechanism | Antibody-tethered Tn5 | Antibody-tethered MNase | Immunoprecipitation |
| Cell input | Very low (10²–10⁴) | Low (10³–10⁵) | High (10⁶–10⁷) |
| Signal-to-noise | Very high | High | Moderate–low |
| Background DNA | Minimal | Low | High |
| Library prep complexity | Simple | Moderate | Complex |
| Sequencing depth required | Low | Low–moderate | High |
| Resolution | High | High | Moderate |
| Fixed cell compatibility | Limited | Limited | Excellent |
| Historical validation | Emerging | Established | Gold standard |

---

## 7. Recommended Use Cases

CUT&Tag is well suited for:
- Low-input samples
- High-resolution chromatin profiling
- Rapid experimental turnaround
- Avoiding sonication-induced artefacts
