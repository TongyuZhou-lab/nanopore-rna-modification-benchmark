# Scorebed Standardization Scripts

This folder contains scripts for converting tool-specific outputs into genome-level score BED files or score-call tables used in benchmark analyses.

The scripts in this folder only perform scorebed standardization. Fixed-threshold binary calls and downstream evaluation metrics are handled separately.

## Folder structure

```text
02_scorebed_standardization/
├── m5C_human/
│   └── build_human_m5c_scorebeds.py
│
├── m6A_arabidopsis/
│   └── build_arabidopsis_m6A_scorebeds.py
│
├── m6A_human/
│   └── build_human_m6A_scorebeds.py
│
└── m6A_mouse/
    └── build_mouse_m6A_scorebeds.py
```

---

## 1. Human m5C

### Script

```text
m5C_human/build_human_m5c_scorebeds.py
```

### Purpose

This script converts human m5C tool outputs into genome-level score BED files.

For TandemMod, the continuous score is defined as:

```text
p95_ratio = count(probability > 0.95) / total
```

This keeps the TandemMod score definition consistent with the m6A analyses.

### Main outputs

```text
TandemMod_human_m5C.GRCh38.centerC.p95_ratio.score.bed
CHEUI_human_m5C.GRCh38.centerPlus5.centerC.score.bed
SWARM_human_m5C.GRCh38.centerPlus5.centerC.score.bed
human_m5C.GRCh38.scorebed_standardization.summary.tsv
```

### Example command

```bash
python m5C_human/build_human_m5c_scorebeds.py \
  --gtf path/to/gencode.v47.primary_assembly.annotation.gtf \
  --tandemmod path/to/TandemMod_m5C_genome_predictions.tsv \
  --cheui path/to/CHEUI_site_level_m5C_predictions.txt \
  --swarm path/to/SWARM_m5C_predictions.tsv \
  --outdir path/to/human_m5C_scorebeds \
  --assembly-label GRCh38
```

---

## 2. Arabidopsis m6A

### Script

```text
m6A_arabidopsis/build_arabidopsis_m6A_scorebeds.py
```

### Purpose

This script standardizes full-depth and downsampled Arabidopsis m6A tool outputs into genome-level score-call tables.

It supports TandemMod, m6Anet, CHEUI, and SWARM across full-depth and downsampled fractions.

### Coordinate and score conventions

```text
TandemMod    p_0.95 / total
m6Anet       probability_modified; transcript_position mapped with origin=0
CHEUI        centerPlus5
SWARM        transcript position + 4, then genome mapping with origin=0
```

### Main outputs

```text
score_call_tables/*.genome_DRACH.score_call.tsv
positive_call_site_lists/*.primary_eval_positive_call_sites.tsv
qc/standardization_summary.primary_eval_space.tsv
qc/full_depth_sanity_check.primary_eval_space.tsv
```

### Input raw-path table

The `--raw-paths` file should be a TSV file with the following columns:

```text
fraction    tool    path
```

Example:

```text
fraction        tool        path
full_depth      TandemMod   path/to/full_depth/TandemMod/m6A.site_level.tsv
full_depth      m6Anet      path/to/full_depth/m6Anet/data.site_proba.csv
full_depth      CHEUI       path/to/full_depth/CHEUI/site_level_m6A_predictions.txt
full_depth      SWARM       path/to/full_depth/SWARM/pred.tsv
frac0.8_seed11  TandemMod   path/to/frac0.8_seed11/TandemMod/m6A.site_level.tsv
```

### Example command

```bash
python m6A_arabidopsis/build_arabidopsis_m6A_scorebeds.py \
  --label-path path/to/labels.primary_nanopore.full_depth_common_universe.tsv \
  --universe-path path/to/full_depth_four_tool_common_universe.site_ids.tsv \
  --raw-paths path/to/arabidopsis_raw_paths.tsv \
  --gtf path/to/Arabidopsis_thaliana.TAIR10.60.gtf \
  --outdir path/to/arabidopsis_m6A_score_call_tables
```

---

## 3. Human m6A

### Script

```text
m6A_human/build_human_m6A_scorebeds.py
```

### Purpose

This script converts human m6A tool outputs into genome-level score BED files.

### Coordinate and score conventions

```text
TandemMod    p_0.95 / total
m6Anet       txPlus1
CHEUI        centerPlus5
SWARM        centerPlus5
```

### Main outputs

```text
TandemMod_human_DRACH.GRCh38.stdchr.score.bed
m6Anet_human_DRACH.GRCh38.txPlus1.stdchr.score.bed
CHEUI_human_all.GRCh38.centerPlus5.stdchr.score.bed
CHEUI_human_DRACH.GRCh38.centerPlus5.stdchr.score.bed
SWARM_human_all.GRCh38.centerPlus5.stdchr.score.bed
SWARM_human_DRACH.GRCh38.centerPlus5.stdchr.score.bed
human_m6A.GRCh38.scorebed_standardization.summary.tsv
```

### Example command

```bash
python m6A_human/build_human_m6A_scorebeds.py \
  --gtf path/to/gencode.v47.primary_assembly.annotation.gtf \
  --tandemmod path/to/TandemMod_m6A.site_level.tsv \
  --m6anet path/to/m6Anet/data.site_proba.csv \
  --cheui-all path/to/CHEUI/site_level_m6A_predictions.txt \
  --cheui-drach path/to/CHEUI/site_level_m6A_predictions.txt.DRACH.tsv \
  --swarm-all path/to/SWARM_m6A_predictions.tsv \
  --swarm-drach path/to/SWARM_m6A_predictions.tsv.DRACH.tsv \
  --outdir path/to/human_m6A_scorebeds \
  --assembly-label GRCh38
```

---

## 4. Mouse m6A

### Script

```text
m6A_mouse/build_mouse_m6A_scorebeds.py
```

### Purpose

This script converts mouse m6A tool outputs into genome-level score BED files.

### Coordinate and score conventions

```text
TandemMod    p_0.95 / total
m6Anet       txPlus1
CHEUI        centerPlus5
SWARM        shiftPlus5
```

### Main outputs

```text
TandemMod_all.mm39.stdchr.score.bed
m6Anet_all.mm39.txPlus1.stdchr.score.bed
CHEUI_all.mm39.centerPlus5.stdchr.score.bed
SWARM_all.mm39.shiftPlus5.stdchr.score.bed
CHEUI_DRACH_all.mm39.centerPlus5.stdchr.score.bed
SWARM_DRACH_all.mm39.shiftPlus5.stdchr.score.bed
mouse_m6A_scorebeds.list
```

### Example command

```bash
python m6A_mouse/build_mouse_m6A_scorebeds.py \
  --gtf path/to/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf \
  --tandemmod path/to/TandemMod_m6A.site_level.tsv \
  --cheui-all path/to/CHEUI/site_level_m6A_predictions.txt \
  --m6anet path/to/m6Anet/data.site_proba.csv \
  --swarm-all path/to/SWARM_m6A_predictions.tsv \
  --cheui-drach path/to/CHEUI/site_level_m6A_predictions.txt.DRACH.tsv \
  --swarm-drach path/to/SWARM_m6A_predictions.tsv.DRACH.tsv \
  --outdir path/to/mouse_m6A_scorebeds \
  --assembly-label mm39 \
  --swarm-shift 5
```

---

## Notes

All score BED files use 0-based half-open genomic coordinates:

```text
chrom    start    end    score
```

The fourth column is a continuous score, not a fixed-threshold binary call.

Fixed-threshold binary outputs, such as `*.binary_score.bed` or `*.fixed_positive.bed`, are not included in this folder.

Note: Different tools report modification positions using different coordinate conventions, such as transcript positions, k-mer start positions, or genome-level positions. Therefore, tool-specific validated coordinate transformations were applied during scorebed standardization rather than forcing a single offset rule across all tools.

The command-line arguments differ across datasets because the raw tool outputs were organized differently for each benchmark setting. Human and mouse datasets use one input file per tool, whereas the Arabidopsis analysis includes multiple sequencing-depth fractions and therefore uses a raw-path manifest table. These differences only reflect dataset-specific file organization; all scripts standardize tool outputs into genome-level score files or score-call tables.
