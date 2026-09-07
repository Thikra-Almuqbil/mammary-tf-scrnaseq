# ==============================================================
# SCRIPT 2 OF 4: NORMALISATION AND DIMENSIONALITY REDUCTION
#
# From a filtered count matrix to a 2D map I can look at.
#
#   PART 1  Log-normalisation, and sctransform with feature selection
#   PART 2  PCA, and choosing how many components to keep
#   PART 3  UMAP, and locating the non-epithelial cells
#
# Why these three belong together: they are all about getting from "counts,
# which are not comparable between cells" to "a small number of dimensions
# that carry the biology". Nothing here removes a cell or makes a claim, it is
# all preparation, and the whole thing runs in one sitting.
#
# HOW TO RUN THIS ONE: from the top, in order. SCTransform in part 1 is the
# slow step, roughly 3 to 8 minutes for 8 layers, and it prints per-layer
# progress. It is working, do not press stop. The checkpoint is saved
# immediately after it, before any of the checks, so the expensive step is
# never lost because a later line errored.
#
# Course reference:
#   Demonstrations/05_NormalisationAndFeatureSelection.R
#   Demonstrations/06_DimensionalityReduction.R (PCA, UMAP)
#   Demonstrations/07_Dataset_Integration.R (split, and JoinLayers)
#
# Expected runtime: about 20 to 30 minutes in total.
# Checkpoints written: RObjects/02a_sct.rds
#                      RObjects/02b_pca.rds
#                      RObjects/02c_umap.rds
# ==============================================================


library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(patchwork)

theme_set(theme_classic())

# Seurat v5 passes objects around internally using the future package, whose
# default size limit is 500 MB. My object is much bigger than that and the
# call fails with a "globals" error without this line. 8 GB is generous.
options(future.globals.maxSize = 8 * 1024^3)

# Set your working directory to the ROOT of this repository before running.

dir.create("results/graphs", recursive = TRUE, showWarnings = FALSE)

qc_seurat_object <- readRDS("RObjects/01c_filtered.rds")
qc_seurat_object



# ==============================================================
# PART 1: NORMALISATION AND FEATURE SELECTION
#
# The problem this solves: two cells can have the same biology and very
# different total counts, just because one was captured more efficiently. Raw
# counts are therefore not comparable between cells, and the difference in
# depth would dominate any distance I calculated. Normalisation removes the
# depth difference so that what is left is biology.
#
# I run two normalisations, deliberately:
#   NormalizeData()  gives the log-normalised RNA assay. I use this whenever
#                    I plot or test one gene at a time, and in script 4.
#   SCTransform()    gives variance-stabilised residuals and picks the 3,000
#                    most variable genes. I use this for PCA and clustering.
# ==============================================================

#### Log-normalisation of the RNA assay ####

# The default shifted-log transformation does three things:
#   1. divide each gene's count by that cell's total UMI count
#   2. multiply by a scale factor, default 10,000, so the units are
#      counts per ten thousand
#   3. take log(x + 1), the +1 so that zeros stay zero instead of going to
#      minus infinity
# (The scale factor is 10,000 rather than a million because at single-cell
# depth a per-million figure would be mostly noise, so the convention is
# per ten thousand.)
qc_seurat_object <- NormalizeData(qc_seurat_object)

# This adds a layer called "data" to the RNA assay, next to "counts"
qc_seurat_object


#### Splitting the RNA assay by sample ####

# Each of my 8 samples was captured on its own 10x run, so in this dataset
# "sample" and "technical batch" are the same thing. Splitting the assay into
# one layer per sample makes SCTransform fit a separate model per sample,
# which is what I want: a lactating sample and a nulliparous sample should
# not share one mean-variance model.
# The course splits by SampleGroup in demonstration 07, but its ARI section
# also compares a SampleName-split version, so both are within scope. I chose
# SampleName because the capture run is the sample, not the stage.
qc_seurat_object[["RNA"]] <- split(qc_seurat_object[["RNA"]],
                                   f = qc_seurat_object$SampleName)

# There should now be 8 counts layers and 8 data layers, 16 in total
qc_seurat_object


#### sctransform ####

# What it does: instead of the log transform, it fits a regularised negative
# binomial model per gene and keeps the Pearson residuals. That removes the
# relationship between a gene's mean and its variance, which the log transform
# leaves behind. The practical effect is that highly expressed genes stop
# dominating the variable-gene list just for being highly expressed.
#
# ONE DELIBERATE DEVIATION FROM THE COURSE.
# The course calls SCTransform with vars.to.regress = "percent.mt". I left it
# out. Regressing a covariate is applied gene by gene across the whole matrix
# and is the most expensive part of the call, and my own diagnostics in script
# 1 showed percent.mt has a median of 0.75 percent and a maximum of 8.2
# percent. There is essentially no variance there to remove, so it was costing
# most of the runtime to accomplish nothing. percent.mt stays in the metadata
# either way, and I TEST this decision explicitly at the end of part 3 by
# plotting percent.mt on the UMAP. If a distinct island had turned out to be a
# mitochondrial hotspot I would have come back and re-run this with the
# regression.
# The course version would be:
#   SCTransform(qc_seurat_object, assay = "RNA",
#               vars.to.regress = "percent.mt", verbose = TRUE)
#
# I am also not regressing out cell cycle, which was on my list of things I
# was allowed to cut for time.

qc_seurat_object <- SCTransform(qc_seurat_object,
                                assay = "RNA",
                                verbose = TRUE)


#### Check it worked, THEN save ####

# I should now see:
#   2 assays, active assay SCT, about 3,000 variable features,
#   SCT layers counts/data/scale.data as well as the 16 RNA layers.
qc_seurat_object

# SCTransform changes the default assay to SCT
DefaultAssay(qc_seurat_object)


#### CHECKPOINT 1 ####

# Saved NOW, before any of the checks below. This object is expensive to
# rebuild and there is no reason to risk it on a typo in a later line.
# Safe stopping point.
saveRDS(qc_seurat_object, "RObjects/02a_sct.rds")


#### The variable features it selected ####

hvgs_sct <- VariableFeatures(qc_seurat_object)

# Should be about 3,000, and definitely not 0
length(hvgs_sct)

# The top 30, ranked by residual variance
head(hvgs_sct, n = 30)


#### Check: did it pick up the genes the paper cares about? ####

# This is my sanity check on the feature selection. These are the genes used
# in the paper's Fig 1d, Fig 2b and Table 1. They are exactly the genes that
# differ between cell types, so if my variable gene selection is working most
# of them should be in the list. If almost none of them were, the selection
# would be picking up noise and I would need to go back.
paper_markers <- c("Krt5", "Krt18",
                   "Esr1", "Prlr", "Pgr", "Cited1",
                   "Aldh1a3", "Kit", "Cd14", "Elf5",
                   "Csn2", "Wap", "Lalba", "Glycam1",
                   "Acta2", "Oxtr", "Krt14", "Pdpn",
                   "Procr", "Zeb2", "Notch3")

# Which of them made it in
intersect(paper_markers, hvgs_sct)

# And which did not. A marker being absent is not automatically a problem: a
# gene expressed in every cell of one large cluster and nowhere else can
# still fall outside the top 3,000 by residual variance.
setdiff(paper_markers, hvgs_sct)


#### Looking ahead: is script 4 going to be feasible? ####

# My research question depends on there being enough transcription factors
# detected in this data to cluster on. I would much rather find out now than
# after building the whole pipeline.

# Confirm the file is where I think it is. If this prints FALSE, fix the path.
file.exists("data/raw/Mus_musculus_TF.txt")

# The AnimalTFDB 4.0 mouse transcription factor list
tf_table <- read_tsv("data/raw/Mus_musculus_TF.txt")

# Unique gene symbols in the list. I got 1,611.
tf_genes <- unique(tf_table$Symbol)
length(tf_genes)

# How many of those are present in my data at all? I got 1,346.
tf_in_data <- intersect(tf_genes, rownames(qc_seurat_object[["RNA"]]))
length(tf_in_data)

# And how many are among the 3,000 variable features?
tf_variable <- intersect(tf_genes, hvgs_sct)
length(tf_variable)

head(tf_variable, n = 40)



# ==============================================================
# PART 2: PCA
#
# Why this comes here: after normalisation I have about 3,000 variable genes,
# so every cell is a point in 3,000-dimensional space. Distances in that many
# dimensions are dominated by noise. PCA finds the directions along which the
# cells actually vary and lets me keep only those, which both denoises the
# data and makes everything downstream tractable.
# ==============================================================

# If you are resuming here:
# qc_seurat_object <- readRDS("RObjects/02a_sct.rds")
# hvgs_sct <- VariableFeatures(qc_seurat_object)
# tf_table <- read_tsv("data/raw/Mus_musculus_TF.txt")
# tf_genes <- unique(tf_table$Symbol)

# PCA has to run on the SCT assay, so I check rather than assume
DefaultAssay(qc_seurat_object)


#### Are transcription factors over-represented among the variable genes? ####

# This is a small piece of my own research question, done here because it is
# cheap and the objects are already loaded. If transcription factors were
# enriched among the genes that drive the variation, that would already be a
# hint about the script 4 answer.

# The pool the 3,000 were chosen FROM. This has to be the SCT assay, not the
# RNA assay, or the test would be against the wrong background.
n_sct_genes <- nrow(qc_seurat_object[["SCT"]])
n_sct_genes

tf_in_sct <- intersect(tf_genes, rownames(qc_seurat_object[["SCT"]]))
length(tf_in_sct)

tf_variable <- intersect(tf_genes, hvgs_sct)
length(tf_variable)

# The two proportions I am comparing
length(tf_in_sct) / n_sct_genes         # TF share of all genes
length(tf_variable) / length(hvgs_sct)  # TF share of the variable genes

# A hypergeometric test of the difference. In words: if I drew 3,000 genes at
# random from the SCT assay, how surprising would it be to get this many
# transcription factors?
#   q = successes observed, minus 1 (see below)
#   m = transcription factors available
#   n = everything else available
#   k = how many genes were drawn
# lower.tail = FALSE gives P(at least this many), and the "minus 1" on q is
# what makes it "at least" rather than "more than".
phyper(q = length(tf_variable) - 1,
       m = length(tf_in_sct),
       n = n_sct_genes - length(tf_in_sct),
       k = length(hvgs_sct),
       lower.tail = FALSE)


#### Running PCA ####

# On the variable features only, not all genes, which is the standard
# approach and what the course does.
qc_seurat_object <- RunPCA(qc_seurat_object,
                           features = VariableFeatures(qc_seurat_object))

Reductions(qc_seurat_object)

# The standard deviation captured by each component, largest first
Stdev(qc_seurat_object, reduction = "pca")


#### Choosing the number of components ####

# Every plot from here on is assigned to a named object and then written to
# disk with ggsave(). I learned this the hard way: printing to the RStudio
# plot pane leaves nothing on disk, so when I wanted a figure for the slides I
# had to re-run the script. Now the figures are made once.
#
# bg = "white" on every ggsave is also something I learned the hard way. On my
# setup some of these came out with a TRANSPARENT background instead of a
# white one, which looks fine in RStudio and fine on a white page, but on any
# dark background the black axis text and titles become invisible. GitHub
# renders READMEs in dark mode by default, so three of my figures were
# unreadable there until I forced the background.
#
# Reading my own elbow plot: there is a clear break after PC4, and the per-PC
# gain drops below about 0.15 from roughly PC20 onwards. I keep 20. This is a
# judgement call, not a rule, and the usual advice is to err on the generous
# side because extra PCs mostly add a little noise whereas too few lose real
# structure.
p_elbow <- ElbowPlot(qc_seurat_object, ndims = 50) +
  labs(title = "PC standard deviation, all cells",
       subtitle = "Clear break after PC4; per-PC gain below ~0.15 from about PC20, so I keep 20")

p_elbow

ggsave("results/graphs/02_elbow_all_cells.png", plot = p_elbow,
       width = 7, height = 5, dpi = 300, bg = "white")


#### What are the first components separating? ####

p_pca_stage <- DimPlot(qc_seurat_object,
                       reduction = "pca",
                       group.by = "SampleGroup") +
  labs(title = "PC1 vs PC2, coloured by developmental stage",
       subtitle = "The four stages occupy different regions of PCA space")

p_pca_stage

ggsave("results/graphs/02_pca_by_stage.png", plot = p_pca_stage,
       width = 8, height = 6, dpi = 300, bg = "white")

# THIS IS THE IMPORTANT ONE. It is the diagnostic behind two decisions: to
# split by sample in part 1, and NOT to run batch correction at all.
# Each stage was sequenced as its own run, so stage and batch are perfectly
# confounded and I cannot separate them by design. What I CAN do is look at
# the two independent mice within each stage. If replicates of one stage sit
# on top of each other while the stages sit apart, then the thing separating
# the stages is not a per-run technical effect, because that would separate
# the replicates too.
p_pca_sample <- DimPlot(qc_seurat_object,
                        reduction = "pca",
                        group.by = "SampleName") +
  labs(title = "Biological replicates co-localise in PCA space",
       subtitle = paste("Two independent mice per stage. Replicate pairs overlap;",
                        "stages separate, so the separation is not a batch effect"))

p_pca_sample

ggsave("results/graphs/02_pca_by_sample.png", plot = p_pca_sample,
       width = 8, height = 6, dpi = 300, bg = "white")

# The two side by side, which is the version I present
p_pca_pair <- p_pca_stage + p_pca_sample +
  plot_annotation(
    title = "Stage separation is biological, not technical",
    subtitle = paste("Left: cells coloured by stage.  Right: the same cells coloured by",
                     "individual mouse. A batch effect would separate the replicates."))

p_pca_pair

ggsave("results/graphs/02_pca_stage_vs_sample.png", plot = p_pca_pair,
       width = 14, height = 6, dpi = 300, bg = "white")

p_pca_23 <- DimPlot(qc_seurat_object,
                    reduction = "pca",
                    dims = c(2, 3),
                    group.by = "SampleGroup") +
  labs(title = "PC2 vs PC3, coloured by developmental stage")

p_pca_23

ggsave("results/graphs/02_pca_pc2_pc3.png", plot = p_pca_23,
       width = 8, height = 6, dpi = 300, bg = "white")

# Which genes load most strongly on the first few components? This is how I
# find out what biological axis each PC represents, rather than treating the
# PCs as anonymous numbers.
print(qc_seurat_object[["pca"]], dims = 1:5, nfeatures = 15)


#### CHECKPOINT 2 ####

# Safe stopping point.
saveRDS(qc_seurat_object, "RObjects/02b_pca.rds")



# ==============================================================
# PART 3: UMAP, AND FINDING THE NON-EPITHELIAL CELLS
#
# I use dims = 1:20, chosen from the elbow plot above.
#
# t-SNE is deliberately NOT run here. The paper's Fig 1b is a t-SNE of 23,184
# cells, that is, AFTER the immune, fibroblast and endothelial cells were
# removed. Running t-SNE now on all 25,802 cells would cost 5 to 10 minutes
# and I would throw the result away. I run t-SNE once, in script 3, on the
# epithelial subset, which is the version that goes in the figure.
# ==============================================================

# If you are resuming here:
# qc_seurat_object <- readRDS("RObjects/02b_pca.rds")


#### Rejoining the RNA layers ####

# The RNA assay is still split into 8 counts layers and 8 data layers from
# part 1. Anything that needs ONE expression matrix, which means FeaturePlot,
# DotPlot and VlnPlot on the RNA assay, errors while it is split. Rejoining is
# the last step of demonstration 07 and it is cheap to split again later if I
# need to.
qc_seurat_object <- JoinLayers(qc_seurat_object, assay = "RNA")

# The RNA assay should now have 2 layers, not 16.
# Careful: printing the object shows the layers of the ACTIVE assay only,
# which is SCT, so I ask for the RNA assay specifically.
Layers(qc_seurat_object[["RNA"]])


#### UMAP ####

# UMAP takes the 20-dimensional PCA space and squashes it into 2 dimensions
# for plotting. It preserves local neighbourhoods reasonably well. Distances
# between far-apart blobs on a UMAP are not meaningful and I do not interpret
# them as such.
qc_seurat_object <- RunUMAP(qc_seurat_object,
                            reduction = "pca",
                            dims = 1:20)

Reductions(qc_seurat_object)

# Save straight away, before the plots
saveRDS(qc_seurat_object, "RObjects/02c_umap.rds")

# I verify the save rather than trusting it. saveRDS on a large Seurat object
# often emits a harmless warning ("no applicable method for 'depth'") that
# comes from object introspection rather than from the write itself, and I did
# not want to assume it was harmless without checking.
file.exists("RObjects/02c_umap.rds")
file.info("RObjects/02c_umap.rds")$size / 1024^2   # size in MB


#### UMAP by stage and by sample ####

p_umap_stage <- DimPlot(qc_seurat_object,
                        reduction = "umap",
                        group.by = "SampleGroup") +
  ggtitle("UMAP coloured by developmental stage (all 25,802 cells)")

p_umap_stage

ggsave("results/graphs/02_umap_by_stage.png",
       plot = p_umap_stage, width = 8, height = 6, dpi = 300, bg = "white")

# One panel per stage, which makes it much easier to see which region belongs
# to which stage than trying to read four overlapping colours
p_umap_split <- DimPlot(qc_seurat_object,
                        reduction = "umap",
                        group.by = "SampleGroup",
                        split.by = "SampleGroup")

p_umap_split

ggsave("results/graphs/02_umap_split_by_stage.png",
       plot = p_umap_split, width = 14, height = 4, dpi = 300, bg = "white")

# Coloured by individual mouse, the replicate check again
p_umap_sample <- DimPlot(qc_seurat_object,
                         reduction = "umap",
                         group.by = "SampleName") +
  ggtitle("UMAP coloured by sample")

p_umap_sample

ggsave("results/graphs/02_umap_by_sample.png",
       plot = p_umap_sample, width = 8, height = 6, dpi = 300, bg = "white")


#### Basal and luminal markers ####

# Krt5 (basal) and Krt18 (luminal) are the two genes in the paper's Fig 1d.
# They should mark opposite ends of the structure, because PC1 came out as a
# basal versus luminal axis.
# I plot from the RNA assay, which is log-normalised, rather than SCT, so the
# values are comparable to the paper's log-transformed normalised counts.
DefaultAssay(qc_seurat_object) <- "RNA"

p_krt <- FeaturePlot(qc_seurat_object,
                     reduction = "umap",
                     features = c("Krt5", "Krt18"))

p_krt

ggsave("results/graphs/02_krt5_krt18.png",
       plot = p_krt, width = 11, height = 5, dpi = 300, bg = "white")


#### Where are the non-epithelial cells? ####

# The paper sorted on EpCAM but still had immune, fibroblast and endothelial
# cells left over, and removed them. I need to find mine. These are the marker
# genes from their Methods:
#   immune       Cd74, Cd52
#   fibroblast   Col1a1, Fn1
#   endothelial  Eng, Emcn, Pecam1
#   epithelial   Epcam, the marker they sorted on
# This is only locating them. The actual removal happens in script 3, after a
# clustering that gives me discrete groups to remove rather than a hand-drawn
# region of a UMAP.
p_contaminants <- FeaturePlot(qc_seurat_object,
                              reduction = "umap",
                              features = c("Epcam", "Cd74",
                                           "Col1a1", "Pecam1"))

p_contaminants

ggsave("results/graphs/02_contaminant_markers.png",
       plot = p_contaminants, width = 10, height = 8, dpi = 300, bg = "white")


#### Did leaving out vars.to.regress cost me anything? ####

# This is the promised test of the decision I made in part 1. If no region of
# the map is a percent.mt hotspot, omitting the regression cost me nothing. If
# a distinct island had lit up, I would have gone back and re-run SCTransform
# with the regression, at the cost of most of an afternoon.
p_qc_metrics <- FeaturePlot(qc_seurat_object,
                            reduction = "umap",
                            features = c("percent.mt", "nCount_RNA"))

p_qc_metrics

ggsave("results/graphs/02_percentmt_ncount.png",
       plot = p_qc_metrics, width = 11, height = 5, dpi = 300, bg = "white")

# Put the default assay back to SCT so that everything downstream behaves.
# Forgetting this line causes confusing errors in the next script.
DefaultAssay(qc_seurat_object) <- "SCT"


#### CHECKPOINT 3: END OF SCRIPT 2 ####

# Saved again, now with the default assay reset. This is the object script 3
# reads.
saveRDS(qc_seurat_object, "RObjects/02c_umap.rds")

list.files("results/graphs")
