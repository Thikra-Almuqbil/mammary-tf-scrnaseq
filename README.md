# Are transcription factors alone enough to define a cell's identity?

A single-cell RNA-seq mini project on the mouse mammary epithelium, done for
the single-cell week of a bioinformatics summer school (Cambridge Bioinformatics
Training / KAUST Academy). Everything here is my own work on a public dataset.

## The short answer

**No, and by a margin I did not expect.**

I clustered the same 23,989 mammary epithelial cells three times, changing only
the genes I gave the algorithm and nothing else:

| Gene set given to the clustering | Genes | Clusters found | Agreement with the full-transcriptome reference (ARI) |
|---|---|---|---|
| All highly variable genes (the reference) | 3,000 | 15 | reference |
| Transcription factors only | 1,346 | 9 | **0.53** |
| Ordinary genes, matched for expression level | 1,346 | 11 | **0.73** |

Transcription factors did **worse** than an equally sparse, equally
hard-to-detect set of ordinary genes. That control arm is what makes the result
mean something: without it, a mediocre TF score could just as easily have meant
"transcription factor mRNAs are too sparsely detected in 10x data to tell".

Looking at where the TF clustering broke down is more interesting than the
number. It merged 6,057 cells into a single cluster, and those cells were the
hormone-sensing cells from **all four developmental stages** plus the luminal
progenitors. So transcription factors kept the **lineage** (hormone-sensing
versus alveolar versus basal) and lost the **developmental state** within that
lineage.

I did not go looking for this answer. I expected transcription factors to do at
least as well as the control, and they did not.

## Why this question

The dataset is Bach et al. (2017), a reference atlas of the mouse mammary gland
across four developmental stages: nulliparous, gestation, lactation and
post-involution. The baseline task for the week was to reproduce their Figures 1
and 2B to 2C, which I did. Then we were each asked to ask something the paper
does not answer.

Transcription factors are the proteins that switch other genes on and off, so
they are the obvious candidate for what makes a cell what it is. That intuition
is widely held and, as far as I could find, not directly tested on this data.
The paper does not test it, because it is an atlas and that is not what atlases
are for.

## What is in this repository

```
data/
  raw/                        the real input data, plus a README explaining
                              which file you have to download yourself
scripts/                      the analysis, 4 R scripts, run in order
results/
  graphs/                     31 figures, all produced by these scripts
  tables/                     the key numbers as CSV
README.md                     this file
REPORT.md                     the full write-up, with the figures embedded,
                              readable without running anything
```

There is no `data/processed/` folder. The scripts do write intermediate objects,
one `saveRDS()` checkpoint per stage, but they total about 18 GB, so there is no
sensible way to put them in a repository. The scripts recreate them, and each
one says at the top what it saves.

## Where the data came from

- **Counts:** Bach et al. (2017), GEO accession
  [GSE106273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106273).
  8 samples, 2 mice per developmental stage, CellRanger version 2 output.
- **Transcription factor list:** AnimalTFDB 4.0, mouse. 1,611 gene symbols.
- **Method:** the Cambridge Bioinformatics Training course material,
  <https://cambiotraining.github.io/single-cell-rnaseq/>. Where I departed from
  it, the script says so and says why.

The 8 count matrices are 12 MB to 31 MB each and GitHub's web uploader rejects
files over 25 MB, so `matrix.mtx.gz` is not in this repository.
`data/raw/README.md` tells you exactly where to get them and where to put them.
Everything else, including the barcode and gene files, is here.

## How to reproduce this

You need R (I used 4.5.2) and these packages:

```r
install.packages(c("tidyverse", "patchwork", "cluster", "mclust"))
install.packages("Seurat")
# Bioconductor
install.packages("BiocManager")
BiocManager::install(c("glmGamPoi"))
# optional, only for part 5 of script 3
install.packages("leidenbase")
```

Then:

1. Download the 8 `matrix.mtx.gz` files as described in `data/raw/README.md`.
2. Open this repository folder as your working directory in R. Every path in
   every script is relative to the repository root, not to `scripts/`.
3. Run the scripts in order. They are numbered, and each reads the checkpoint
   the previous one wrote:

   | Script | What it does | Roughly |
   |---|---|---|
   | `01_import_and_quality_control.R` | read the 8 matrices, calculate QC metrics, compare three threshold schemes, apply the one I chose | 15 min |
   | `02_normalisation_and_dimension_reduction.R` | log-normalisation, SCTransform, PCA, UMAP | 20 to 30 min |
   | `03_clustering_and_cell_type_annotation.R` | remove the non-epithelial cells, tune the clustering, name the 15 cell states. **Figures 1b, 1c, 1d, 2b, 2c** | 70 min |
   | `04_transcription_factor_analysis.R` | **the research question** | 25 min |

   Total about 2 hours 20 if nothing goes wrong.

   Each script is divided into numbered PARTs with a header comment saying what
   that part is for, and each part ends with a `saveRDS()` checkpoint and a
   "safe stopping point" note. You do not have to run a whole script in one
   sitting: every part says which checkpoint to read if you are resuming.

   Two parts are optional and can be skipped entirely. Part 5 of script 3 is
   the Louvain versus Leiden check, which needs an extra package that will not
   install everywhere. Part 8 of script 3 is the ranked marker list, which is
   slow and which the annotation does not depend on.

4. Figures land in `results/graphs/`, tables in `results/tables/`, and
   checkpoints in an `RObjects/` folder that the scripts create and `.gitignore`
   excludes.

You can also just read `REPORT.md`, which has all the figures in it and requires
nothing at all.

## Two things I would tell you before you read the code

**The quality control is per sample, and that is the whole point.** The 8 mice
are in four very different physiological states. A lactating mammary gland is
essentially a milk protein factory, so those cells genuinely have different
library sizes and different gene counts from a nulliparous mouse. One global
threshold would preferentially delete lactation cells for biological reasons and
then bias every stage comparison downstream, including my own research question.
So every threshold is computed inside its own sample. I also compared three
threshold schemes before choosing, in part 3 of script 1, because a MAD-based cutoff is
*relative*: it always removes a slice off the bottom of whatever distribution
you hand it, whether or not anything is actually wrong.

**The scripts are written to match the course practicals, not to be elegant.**
Explicit and repetitive, no custom functions, no clever pipelines. If the
practical wrote five near-identical lines I wrote five lines. That was a
deliberate constraint on this project and it is why the code looks the way it
does.

## Honest limitations

Written out properly in `REPORT.md`, but the short version:

- Cell type and developmental stage are **confounded by design** in this
  dataset. Almost every cluster comes from one stage, so I cannot cleanly
  separate "this is a different cell type" from "this is the same cell type at a
  different stage".
- My clustering split one myoepithelial population in two along a **technical**
  axis: cluster 7 is the standout for dissociation-stress genes. That is a
  limitation of my analysis, not a finding.
- One cluster (c10) is genuinely ambiguous between Procr+ basal and pericyte.
  The original paper had the same ambiguity and resolved it the same way, by
  keeping it and flagging it.
- The negative TF result is a statement about **clustering on TF mRNA in 10x
  data**. It is not a statement about transcription factor protein activity,
  which is what actually does the regulating and which this assay does not
  measure.

## Author

Thikra Alwazzan. Bioinformatics summer school, single-cell week, 2026.
