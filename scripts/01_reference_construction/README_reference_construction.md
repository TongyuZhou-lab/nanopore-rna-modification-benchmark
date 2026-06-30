# Reference Construction Scripts

This folder contains scripts used to construct reference sets for mouse m6A, human m6A, and human m5C benchmark analyses.

## Folder structure

```text
01_reference_construction/
├── m6A_mouse/
├── m6A_human/
└── m5C_human/
```

---

## 1. Mouse m6A

### Scripts

```text
m6A_mouse/build_glori_discrete_truth.py
m6A_mouse/build_mouse_reference_sets.py
```

### Purpose

`build_glori_discrete_truth.py` converts two GLORI mouse replicates into loose and strict GLORI BED files.

`build_mouse_reference_sets.py` combines m6A-CLIP, m6Aboost, and GLORI reference sites to generate the final mouse m6A reference sets.

### Main outputs

```text
strict_consensus_2of3.bed
loose_union_3way.bed
reference_set_counts.tsv
```

### Example

```bash
python m6A_mouse/build_glori_discrete_truth.py \
  --rep1 path/to/GLORI_mouse_rep1.csv.gz \
  --rep2 path/to/GLORI_mouse_rep2.csv.gz \
  --outdir path/to/glori_output \
  --prefix GLORI_mESC
```

```bash
python m6A_mouse/build_mouse_reference_sets.py \
  --truth-clip path/to/m6A_CLIP_reference.bed \
  --truth-m6aboost path/to/m6Aboost_reference.bed \
  --truth-glori02 path/to/GLORI_meanRatio0.2.bed \
  --truth-glori095 path/to/GLORI_meanRatio0.95.bed \
  --outdir path/to/reference_sets
```

---

## 2. Human m6A

### Script

```text
m6A_human/build_human_m6A_truth_DRACH.py
```

### Purpose

This script builds human HEK293T m6A DRACH reference sets from GLORI, m6Aboost, and m6ACE sources.

### Main outputs

```text
truth_sources_DRACH/GLORI.strict.bed
truth_sources_DRACH/GLORI.loose.bed
reference_sets_DRACH/strict_consensus_2of3.bed
reference_sets_DRACH/loose_union_3way.bed
reference_sets_DRACH/reference_set_counts.tsv
```

### Example

```bash
python m6A_human/build_human_m6A_truth_DRACH.py \
  --genome-fa path/to/GRCh38.genome.fa.gz \
  --glori-rep1 path/to/GLORI_human_rep1.csv.gz \
  --glori-rep2 path/to/GLORI_human_rep2.csv.gz \
  --m6aboost-bed path/to/m6Aboost_HEK293T.bed.gz \
  --m6ace-csv path/to/m6ACE-Seq.csv \
  --outdir path/to/human_m6A_reference
```

---

## 3. Human m5C

### Script

```text
m5C_human/build_human_m5c_reference_sets.py
```

### Purpose

This script prepares human HeLa m5C reference sets from three public reference sources and merges them into strict and loose reference sets.

### Main outputs

```text
strict_consensus_2of3.bed
loose_union_3way.bed
reference_set_counts.tsv
```

### Example

```bash
python m5C_human/build_human_m5c_reference_sets.py \
  --gse122260 path/to/GSE122260_Ctrl.GRCh38.bed \
  --gse140995 path/to/GSE140995_transcriptome_wide.GRCh38.bed \
  --gse93749 path/to/GSE93749_HeLa_siCTRL_final.GRCh38.bed \
  --outdir path/to/human_m5C_reference_sets
```

---

## 4. Arabidopsis m6A

Arabidopsis m6A uses a pre-defined binary reference set directly. No multi-source consensus reference construction script is required here.

---

## Notes

All BED files should use 0-based half-open genomic coordinates:

```text
chrom    start    end
```

Large raw sequencing files and intermediate tool outputs are not included in this repository.
