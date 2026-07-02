# Nanopore RNA Modification Benchmark

This repository provides the supplementary document and representative running scripts for our BIBM manuscript:

**Benchmarking Computational Methods for RNA Methylation Detection from Nanopore Direct RNA Sequencing**

## Overview

Nanopore direct RNA sequencing (DRS) enables direct profiling of native RNA molecules and provides a powerful platform for RNA methylation detection. In this study, we benchmarked four representative DRS-based RNA methylation callers, including **m6Anet**, **CHEUI**, **TandemMod**, and **SWARM**, across multiple species, RNA modification types, and analytical settings.

The benchmark includes:

- Cross-species m6A detection in human, mouse, and *Arabidopsis thaliana* datasets;
- Cross-modification evaluation of m5C detection in human data;
- Fixed-threshold and threshold-independent performance evaluation;
- Continuous modification-level agreement with external reference datasets;
- Read-depth robustness analysis using downsampled direct RNA sequencing data;
- Tool-specific parameter sensitivity analysis for SWARM and m6Anet.

This repository is intended to improve the transparency and reproducibility of the analyses reported in the manuscript.

## Repository Structure

```text
.
├── README.md
├── supplementary.pdf
├── .gitignore
└── scripts/
    ├── 01_reference_construction/
    ├── 02_scorebed_standardization/
    ├── 03_main_benchmark/
    ├── 04_downsampling/
    └── 05_swarm_analysis/
