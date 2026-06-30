# Cross-species m6A ROC/PR plotting

This script generates unified-style ROC and PR figures for cross-species m6A evaluation.

It supports three species:

```text
Mouse
Human
Arabidopsis
```

and four tools:

```text
TandemMod
m6Anet
CHEUI
SWARM
```

## Script

```text
plot_m6a_cross_species_rocp.py
```

## What this script does

The script reads or computes ROC/PR curve points and metric tables, then plots full and clean versions of each figure.

Main outputs:

```text
OUTROOT/
├── full_version/
├── clean_version/
└── tables/
```

The `full_version/` figures contain titles, legends, and axis labels.

The `clean_version/` figures remove titles and legends. By default, clean figures also remove x/y axis labels while keeping axes and tick numbers.

## Input logic

### Mouse

Mouse ROC/PR curve points and metrics are read directly from existing tables:

```text
--mouse-curve
--mouse-metrics
```

### Arabidopsis

Arabidopsis can be handled in either of two ways.

Recommended: provide existing ROC/PR curve and metric tables:

```text
--arab-curve
--arab-metrics
```

Alternatively, provide a labeled input table and let the script recompute curves:

```text
--arab-input
```

If recomputing from `--arab-input`, the input table should contain a label column and score columns for TandemMod, m6Anet, CHEUI, and SWARM.

Default Arabidopsis columns:

```text
label
mod_rate_p95
m6Anet_txShift0
CHEUI_origin1_shiftMinus5
SWARM_txCenterPlus5
```

These can be changed using command-line arguments.

### Human

Human ROC/PR curves are recomputed from:

```text
--human-common-universe
--human-truth-bed
--human-loose-bed
--human-tandemmod-scorebed
--human-m6anet-scorebed
--human-cheui-scorebed
--human-swarm-scorebed
```

The script uses the truth BED as positive sites and removes pseudo-negative sites that fall near the loose-union BED within the selected buffer.

Default human truth settings:

```text
truth set: strict_consensus_2of3
loose-union exclusion buffer: 2 bp
```

## Basic usage

```bash
python plot_m6a_cross_species_rocp.py \
  --outroot path/to/cross_species_m6A_ROCPR_figures \
  --mouse-curve path/to/ROC_PR_curve_points_common.tsv \
  --mouse-metrics path/to/ROC_PR_metrics_common.tsv \
  --arab-curve path/to/ROC_PR_curve_points.primary_nanopore.TandemMod_p95rate.tsv \
  --arab-metrics path/to/ROC_PR_metrics.primary_nanopore.TandemMod_p95rate.tsv \
  --human-common-universe path/to/human_four_tool_common_universe.site_ids.tsv \
  --human-truth-bed path/to/strict_consensus_2of3.bed \
  --human-loose-bed path/to/loose_union_3way.bed \
  --human-tandemmod-scorebed path/to/TandemMod_human_DRACH.GRCh38.stdchr.siteMaxP95.score.bed \
  --human-m6anet-scorebed path/to/m6Anet_human_DRACH.GRCh38.txPlus1.stdchr.score.bed \
  --human-cheui-scorebed path/to/CHEUI_human_DRACH.GRCh38.centerPlus5.stdchr.score.bed \
  --human-swarm-scorebed path/to/SWARM_human_DRACH.GRCh38.centerPlus5.stdchr.score.bed
```

## Example with Arabidopsis recomputation

Use this if you do not already have Arabidopsis ROC/PR curve-point and metric tables.

```bash
python plot_m6a_cross_species_rocp.py \
  --outroot path/to/cross_species_m6A_ROCPR_figures \
  --mouse-curve path/to/ROC_PR_curve_points_common.tsv \
  --mouse-metrics path/to/ROC_PR_metrics_common.tsv \
  --arab-input path/to/common_sites_with_TandemMod_official_call.primary_nanopore.tsv \
  --arab-label-col label \
  --arab-tandemmod-col mod_rate_p95 \
  --arab-m6anet-col m6Anet_txShift0 \
  --arab-cheui-col CHEUI_origin1_shiftMinus5 \
  --arab-swarm-col SWARM_txCenterPlus5 \
  --human-common-universe path/to/human_four_tool_common_universe.site_ids.tsv \
  --human-truth-bed path/to/strict_consensus_2of3.bed \
  --human-loose-bed path/to/loose_union_3way.bed \
  --human-tandemmod-scorebed path/to/TandemMod.score.bed \
  --human-m6anet-scorebed path/to/m6Anet.score.bed \
  --human-cheui-scorebed path/to/CHEUI.score.bed \
  --human-swarm-scorebed path/to/SWARM.score.bed
```

## Main outputs

For each species, the script generates ROC and PR figures.

Example output files:

```text
full_version/mouse_ROC.full.svg
full_version/mouse_PR.full.svg
full_version/human_ROC.full.svg
full_version/human_PR.full.svg
full_version/arabidopsis_ROC.full.svg
full_version/arabidopsis_PR.full.svg

clean_version/mouse_ROC.clean.svg
clean_version/mouse_PR.clean.svg
clean_version/human_ROC.clean.svg
clean_version/human_PR.clean.svg
clean_version/arabidopsis_ROC.clean.svg
clean_version/arabidopsis_PR.clean.svg
```

The script also writes the exact tables used for plotting:

```text
tables/mouse.curve_points.used_for_plot.tsv
tables/mouse.metrics.used_for_plot.tsv
tables/human.curve_points.used_for_plot.tsv
tables/human.metrics.used_for_plot.tsv
tables/arabidopsis.curve_points.used_for_plot.tsv
tables/arabidopsis.metrics.used_for_plot.tsv
tables/cross_species_m6A_ROCPR_metrics.used_for_plot.tsv
tables/cross_species_m6A_ROCPR_source_files.tsv
```

## Common style options

```bash
--fig-width 4.0
--fig-height 4.3
--line-width 1.8
--tick-labelsize 12
--axis-labelsize 9
--legend-size 7.5
--dpi 300
```

The default tick positions are:

```text
0, 0.25, 0.5, 0.75, 1.0
```

You can change them with:

```bash
--tick-values 0,0.25,0.5,0.75,1.0
```

## Clean figure options

By default, clean figures remove axis labels but keep tick numbers.

To keep axis labels in clean figures:

```bash
--clean-keep-axis-labels
```

By default, the top and right spines are removed. To keep the full box around each plot:

```bash
--keep-full-box
```

## Notes

This script is intended for final unified plotting of cross-species m6A ROC/PR results.

It is not responsible for upstream standardization of raw tool outputs. Prepare the mouse, Arabidopsis, and human score/curve inputs before running this script.
