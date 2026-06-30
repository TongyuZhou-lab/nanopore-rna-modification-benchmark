# Arabidopsis downsampling analysis

This folder contains two scripts for the Arabidopsis m6A downsampling / read-depth robustness analysis.

## Scripts

### 1. `standardize_arabidopsis_downsample_score_calls.py`

This script standardizes raw outputs from four tools into genome-level score/call tables.

Supported tools:

```text
TandemMod
m6Anet
CHEUI
SWARM
```

It processes full-depth and downsampled results, maps transcript-level outputs to genome coordinates when needed, applies the fixed tool-specific cutoffs, and writes standardized outputs for downstream analysis.

Main outputs:

```text
OUTDIR/
├── score_call_tables/
│   ├── full_depth.TandemMod.genome_DRACH.score_call.tsv
│   ├── frac0.8_seed11.m6Anet.genome_DRACH.score_call.tsv
│   └── ...
├── positive_call_site_lists/
│   ├── full_depth.TandemMod.primary_eval_positive_call_sites.tsv
│   └── ...
└── qc/
    ├── standardization_summary.primary_eval_space.tsv
    └── full_depth_sanity_check.primary_eval_space.tsv   # if expected counts are provided
```

### Required input files

The script needs:

```text
1. label table
2. full-depth common-universe site list
3. raw-path table
4. Arabidopsis GTF file
5. output directory
```

The raw-path table should be a TSV file with three columns:

```text
fraction    tool    path
```

Example:

```text
fraction        tool        path
full_depth      TandemMod   /path/to/full_depth/TandemMod/m6A.site_level.tsv
full_depth      m6Anet      /path/to/full_depth/m6Anet/data.site_proba.csv
full_depth      CHEUI       /path/to/full_depth/CHEUI/site_level_m6A_predictions.txt
full_depth      SWARM       /path/to/full_depth/SWARM/pred.tsv
frac0.8_seed11  TandemMod   /path/to/frac0.8_seed11/TandemMod/m6A.site_level.tsv
frac0.8_seed11  m6Anet      /path/to/frac0.8_seed11/m6Anet/data.site_proba.csv
frac0.8_seed11  CHEUI       /path/to/frac0.8_seed11/CHEUI/site_level_m6A_predictions.txt
frac0.8_seed11  SWARM       /path/to/frac0.8_seed11/SWARM/pred.tsv
```

### Example command

```bash
python standardize_arabidopsis_downsample_score_calls.py \
  --label-path path/to/labels.primary_nanopore.full_depth_common_universe.tsv \
  --universe-path path/to/full_depth_four_tool_common_universe.site_ids.tsv \
  --raw-paths path/to/arabidopsis_downsample_raw_paths.tsv \
  --gtf path/to/Arabidopsis_thaliana.TAIR10.60.gtf \
  --outdir path/to/standardized_full_and_downsample_score_call_tables
```

Optional expected-count sanity check:

```bash
python standardize_arabidopsis_downsample_score_calls.py \
  --label-path path/to/labels.primary_nanopore.full_depth_common_universe.tsv \
  --universe-path path/to/full_depth_four_tool_common_universe.site_ids.tsv \
  --raw-paths path/to/arabidopsis_downsample_raw_paths.tsv \
  --gtf path/to/Arabidopsis_thaliana.TAIR10.60.gtf \
  --expected-full-depth path/to/expected_full_depth_primary.tsv \
  --outdir path/to/standardized_full_and_downsample_score_call_tables
```

The optional expected-count table should contain:

```text
tool    expected_n_primary_eval_positive_calls    expected_TP_primary
TandemMod   52      15
m6Anet      760     251
CHEUI       577     177
SWARM       1597    382
```

Default fixed cutoffs:

```text
TandemMod: p_0.95 > 10 and p_0.95 / total >= 0.2
m6Anet: score >= 0.9
CHEUI: score >= 0.9999
SWARM: probability >= 0.9972 and stoichiometry > 0.1
```

These can be changed using:

```bash
--tandemmod-p95-min 10
--tandemmod-ratio-min 0.2
--m6anet-threshold 0.9
--cheui-threshold 0.9999
--swarm-prob-threshold 0.9972
--swarm-stoich-threshold 0.1
```

---

### 2. `plot_arabidopsis_downsampling.py`

This script generates the final unified-style Arabidopsis downsampling figures from prepared summary tables.

It does not standardize raw tool outputs. It expects already summarized downsampling metric tables.

Main panels:

```text
A. Reference-hit recall
B. Full-depth TP-site retention
C. Recovered reference-hit count
D. Positive-call composition
```

Main outputs:

```text
OUTDIR/
├── full/
├── clean/
└── tables/
```

`full/` figures contain titles, legends, and axis labels.  
`clean/` figures remove titles and legends while keeping axes and tick labels by default.

### Example command

```bash
python plot_arabidopsis_downsampling.py \
  --acd-table path/to/downsample_four_figure_metrics.v3.tsv \
  --b-table path/to/B_fullDepth_TPsite_retention.v3.tsv \
  --outdir path/to/arabidopsis_downsampling_figures
```

Common style options:

```bash
--fig-width 4.0
--fig-height 4.3
--tick-labelsize 12
--axis-labelsize 12
--dpi 300
```

If you want the clean figures to remove x/y axis labels but keep tick numbers:

```bash
python plot_arabidopsis_downsampling.py \
  --acd-table path/to/downsample_four_figure_metrics.v3.tsv \
  --b-table path/to/B_fullDepth_TPsite_retention.v3.tsv \
  --outdir path/to/arabidopsis_downsampling_figures \
  --clean-hide-axis-labels
```

---

## Recommended workflow

First standardize the raw full-depth and downsample outputs:

```bash
python standardize_arabidopsis_downsample_score_calls.py \
  --label-path path/to/labels.primary_nanopore.full_depth_common_universe.tsv \
  --universe-path path/to/full_depth_four_tool_common_universe.site_ids.tsv \
  --raw-paths path/to/arabidopsis_downsample_raw_paths.tsv \
  --gtf path/to/Arabidopsis_thaliana.TAIR10.60.gtf \
  --outdir path/to/standardized_full_and_downsample_score_call_tables
```

Then use the prepared summary tables to draw the final figures:

```bash
python plot_arabidopsis_downsampling.py \
  --acd-table path/to/downsample_four_figure_metrics.v3.tsv \
  --b-table path/to/B_fullDepth_TPsite_retention.v3.tsv \
  --outdir path/to/arabidopsis_downsampling_figures
```

## Notes

`standardize_arabidopsis_downsample_score_calls.py` is the upstream preprocessing script.

`plot_arabidopsis_downsampling.py` is the final plotting script.

For manuscript figures, use the outputs from `plot_arabidopsis_downsampling.py`.
