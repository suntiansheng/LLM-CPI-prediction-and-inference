# Code Release Workspace

Companion code release for the paper in `https://arxiv.org/abs/2506.09516`. This repository
contains analysis workflows and utilities for topic modeling and forecasting on
WSJ and Weibo data. Most analyses are written in R and rendered from R Markdown;
a few experiments use Python for neural forecasting.

## Repository Layout
- `WSJ/`: WSJ topic modeling, robustness checks, plotting, and simulations.
- `Weibo/`: Weibo topic modeling, COVID analyses, plotting, and neural forecasts.
- `Joint-TS/`: Joint time-series experiments and logs.
- `weibo-short-version-release/`: Short-release subset of Weibo-derived data.


## Data Files
Data files are stored alongside their analyses. Some scripts expect the working
directory to be the analysis folder. The Weibo release subset is available under
`weibo-short-version-release/`.

## Data Release Note
Our release dataset is limited to Weibo-derived data. The Weibo release data in
this repository span 2019-01-01 to 2023-12-31 (see the `time` column in the CSV
files). The Weibo data range in this repository is shorter than the manuscript,
but results in `weibo-short-version-release/` support the manuscript's findings.
The analytic code is fully
documented in executable R Markdown files and publicly available in this Repository. The underlying
data sets include proprietary Wall Street Journal content and Weibo posts
accessed under platform terms. Because the WSJ terms of use may treat certain
derived data sets (e.g., embeddings) as reproductions, creating potential legal
uncertainty, we do not share the WSJ data. We release only non-expressive
derived data sets from Weibo to support reproducibility.

## Dataset License
The dataset files in `weibo-short-version-release/` are governed by
`DATASET_LICENSE.md` and require contacting the author before publication use.
