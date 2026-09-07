# ==============================================================
# SCRIPT 4 OF 4: THE RESEARCH QUESTION
#
#   Is a mammary epithelial cell's identity readable from its
#   transcription factors alone?
#
#   PART 1  Load the transcription factor list, and measure how detectable
#           transcription factors actually are
#   PART 2  Build the expression-matched control gene set
#   PART 3  Cluster on transcription factors, and on the control set
#   PART 4  The result: adjusted Rand index, and WHERE it broke down
#
# Why this is a question worth asking: transcription factors are the proteins
# that switch other genes on and off, so they are the obvious candidate for
# what makes a cell what it is. That is a widely held intuition. Nobody in the
# Bach et al. paper tests it, because the paper is a reference atlas and this
# is not what atlases are for. So it is a real gap and it is answerable with
# the object script 3 already produced.
#
# THE DESIGN. Three clusterings of the same 23,989 cells:
#
#   REFERENCE  clustering on the 3,000 highly variable genes.
#              Already computed in script 3 (final_cluster / cluster_label).
#
#   TF-ONLY    clustering on all 1,346 AnimalTFDB transcription factors
#              present in the data. NO variance filtering is applied to the
#              TF list, deliberately: my question is whether TF expression
#              recovers identity, not whether the most variable TFs do.
#
#   CONTROL    clustering on the same NUMBER of non-TF genes, sampled to
#              match the TFs' expression distribution.
#
# WHY THE CONTROL ARM IS THE WHOLE EXPERIMENT.
# Transcription factor mRNAs are low abundance, so a mediocre TF result on its
# own is ambiguous: it could mean "transcription factors do not encode cell
# identity", or it could just mean "transcription factors are too sparsely
# detected in 10x data to tell". Those are completely different conclusions
# and I could not distinguish them with two arms.
# The control set is drawn from the same expression bins, so it suffers the
# same dropout. Comparing TF against expression-matched control cancels most
# of the detectability problem and turns an uninterpretable negative into an
# interpretable one. If I only had time for one clever thing in this project,
# this was it.
#
# This script is an EXTENSION beyond the course. The course does not cover
# clustering on a chosen gene set, or building a matched control set. I kept
# the added code as plain as I could: a for loop over 20 bins and nothing else.
#
# Course reference for the pieces that ARE from the course:
#   Demonstrations/06_DimensionalityReduction.R (RunPCA, RunUMAP)
#   Demonstrations/08_Clustering.R (FindNeighbors, FindClusters)
#   Demonstrations/07_Dataset_Integration.R (adjusted Rand index)
#
# Expected runtime: about 25 minutes.
# Checkpoints written: RObjects/04_tf.rds (re-saved after each slow step)
# ==============================================================


library(Seurat)
library(tidyverse)
library(patchwork)

theme_set(theme_classic())

options(future.globals.maxSize = 8 * 1024^3)

# Set your working directory to the ROOT of this repository before running.

dir.create("results/graphs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

epi_seurat_object <- readRDS("RObjects/03f_annotated.rds")

# Everything here uses the log-normalised RNA assay so that all three gene
# sets are treated identically. Using SCT would restrict me to the genes
# SCTransform could scale, which excludes some transcription factors (Spi1,
# for one), and then the TF arm would be handicapped for a reason that has
# nothing to do with my question.
DefaultAssay(epi_seurat_object) <- "RNA"

table(epi_seurat_object$cluster_label)



# ==============================================================
# PART 1: THE TRANSCRIPTION FACTOR LIST, AND HOW DETECTABLE IT IS
# ==============================================================

tf_table <- read_tsv("data/raw/Mus_musculus_TF.txt")

# 1,611 unique symbols in AnimalTFDB 4.0 for mouse
tf_genes_all <- unique(tf_table$Symbol)
length(tf_genes_all)

# Keep the ones actually present in my data. I got 1,346, which is 84 percent
# of the list, so nothing systematic went wrong with symbol matching.
tf_genes <- intersect(tf_genes_all, rownames(epi_seurat_object[["RNA"]]))
length(tf_genes)


#### How detectable are transcription factors, really? ####

# I measure this rather than asserting it, and I need the per-gene means
# anyway to build the control set in part 2.

# The log-normalised expression matrix, genes by cells
rna_data <- GetAssayData(epi_seurat_object, assay = "RNA", layer = "data")

# Mean log-normalised expression per gene, across all 23,989 cells
gene_means <- rowMeans(rna_data)

# The proportion of cells in which each gene is detected at all. This is the
# direct measure of dropout: a gene detected in 3 percent of cells cannot
# contribute much to a distance calculation whatever its biology.
gene_detection <- rowMeans(rna_data > 0)

gene_info <- data.frame(gene = rownames(rna_data),
                        mean_expr = gene_means,
                        detection = gene_detection)

gene_info$is_tf <- gene_info$gene %in% tf_genes

# Transcription factors against everything else. This is the number that
# justifies the control arm existing.
gene_info %>%
  group_by(is_tf) %>%
  summarise(n = n(),
            median_mean_expr = median(mean_expr),
            median_detection = median(detection))



# ==============================================================
# PART 2: BUILDING THE EXPRESSION-MATCHED CONTROL SET
#
# The idea: for every transcription factor, find a non-TF gene expressed at a
# similar level, so the control set as a whole has the same expression profile
# as the TF set and therefore the same dropout problem.
# ==============================================================

# ntile() splits all genes into 20 bins of equal SIZE by mean expression, so
# bin 1 is the lowest-expressed 5 percent and bin 20 the highest.
gene_info$bin <- ntile(gene_info$mean_expr, 20)

# How many transcription factors land in each bin? They are concentrated in
# the low bins, which is the point.
table(gene_info$bin[gene_info$is_tf])

# For each bin, take the same number of non-TF genes.
# I wrote this as a plain for loop over the 20 bins rather than anything
# clever, because I want to be able to explain it line by line.
set.seed(123)

control_genes <- c()

for (b in 1:20) {

  # how many transcription factors are in this expression bin
  n_tf_in_bin <- sum(gene_info$is_tf & gene_info$bin == b)

  # the non-TF genes available to draw from in this bin
  pool <- gene_info$gene[!gene_info$is_tf & gene_info$bin == b]

  # take that many, or all of them if the pool happens to be smaller
  n_to_take <- min(n_tf_in_bin, length(pool))

  control_genes <- c(control_genes, sample(pool, n_to_take))
}

# Should be 1,346, the same as the TF set
length(control_genes)

# THE CHECK THAT THE MATCHING WORKED. If the two sets had very different
# median expression or detection, the control would not be a control and the
# whole comparison would be worthless. They came out close.
gene_info %>%
  mutate(set = case_when(gene %in% tf_genes ~ "transcription factors",
                         gene %in% control_genes ~ "matched control",
                         TRUE ~ "other")) %>%
  filter(set != "other") %>%
  group_by(set) %>%
  summarise(n = n(),
            median_mean_expr = median(mean_expr),
            median_detection = median(detection))


#### The same check, as the figure I would show if challenged ####

# If someone says the comparison was unfair, this is the answer: the two
# density curves sit on top of each other, so the two gene sets are matched
# for exactly the property that would otherwise confound the result.
p_matching <- gene_info %>%
  mutate(set = case_when(gene %in% tf_genes ~ "Transcription factors",
                         gene %in% control_genes ~ "Matched control",
                         TRUE ~ "All other genes")) %>%
  ggplot(aes(x = mean_expr, colour = set)) +
  geom_density(linewidth = 1) +
  scale_x_log10() +
  labs(title = "The control gene set is matched to the transcription factors for expression",
       subtitle = "Both sets experience comparable dropout, so the comparison is not confounded by detectability",
       x = "Mean log-normalised expression (log scale)",
       y = "Density",
       colour = NULL)

p_matching

ggsave("results/graphs/04_expression_matching.png",
       plot = p_matching, width = 9, height = 6, dpi = 300)



# ==============================================================
# PART 3: THE TWO ALTERNATIVE CLUSTERINGS
# ==============================================================

#### Scaling the genes I am about to use ####

# RunPCA reads the scale.data layer, so both gene sets have to be scaled
# first. I scale them in one call rather than two, to save time.
genes_to_scale <- union(tf_genes, control_genes)
length(genes_to_scale)

epi_seurat_object <- ScaleData(epi_seurat_object, features = genes_to_scale)


#### ARM 2: clustering on transcription factors only ####

# PCA on the TF genes alone. reduction.name keeps this separate from the
# reference PCA computed in script 3, so I do not overwrite the thing I am
# comparing against.
epi_seurat_object <- RunPCA(epi_seurat_object,
                            features = tf_genes,
                            reduction.name = "pca_tf",
                            reduction.key = "PCTF_")

# Same k, same dims, same resolution as the reference clustering. This is
# essential: if I changed any of them, a difference in the answer could just
# be the parameters. graph.name keeps this neighbour graph separate too.
epi_seurat_object <- FindNeighbors(epi_seurat_object,
                                   reduction = "pca_tf",
                                   dims = 1:20,
                                   k.param = 25,
                                   graph.name = c("tf_nn", "tf_snn"))

epi_seurat_object <- FindClusters(epi_seurat_object,
                                  graph.name = "tf_snn",
                                  resolution = 0.4,
                                  algorithm = 1,
                                  cluster.name = "tf_cluster")

# How many clusters did the transcription factors alone find? I got 9, against
# 15 for the reference. So TFs are already merging things.
table(epi_seurat_object$tf_cluster)

saveRDS(epi_seurat_object, "RObjects/04_tf.rds")


#### ARM 3: clustering on the expression-matched control genes ####

epi_seurat_object <- RunPCA(epi_seurat_object,
                            features = control_genes,
                            reduction.name = "pca_ctrl",
                            reduction.key = "PCCTRL_")

epi_seurat_object <- FindNeighbors(epi_seurat_object,
                                   reduction = "pca_ctrl",
                                   dims = 1:20,
                                   k.param = 25,
                                   graph.name = c("ctrl_nn", "ctrl_snn"))

epi_seurat_object <- FindClusters(epi_seurat_object,
                                  graph.name = "ctrl_snn",
                                  resolution = 0.4,
                                  algorithm = 1,
                                  cluster.name = "ctrl_cluster")

# I got 11 clusters here
table(epi_seurat_object$ctrl_cluster)

saveRDS(epi_seurat_object, "RObjects/04_tf.rds")



# ==============================================================
# PART 4: THE RESULT
# ==============================================================

#### Adjusted Rand index against the reference clustering ####

# 1 = identical partitions, 0 = agreement no better than chance.
# It is pairwise and label-free, so it does not matter that TF cluster 3 and
# reference cluster 3 are unrelated names.
#
# I use mclust::adjustedRandIndex() rather than the course's
# bluster::pairwiseRand() because bluster is source-only on Windows for
# Bioconductor 3.22 and needs a compiler. The statistic is identical. I call
# it with :: rather than attaching mclust, because attaching mclust masks
# dplyr::count and purrr::map.

# Transcription factors versus the full-transcriptome reference. I got 0.533.
ari_tf <- mclust::adjustedRandIndex(epi_seurat_object$final_cluster,
                                    epi_seurat_object$tf_cluster)
ari_tf

# Expression-matched control versus the same reference. I got 0.728.
ari_ctrl <- mclust::adjustedRandIndex(epi_seurat_object$final_cluster,
                                      epi_seurat_object$ctrl_cluster)
ari_ctrl

# And the two alternatives against each other. I got 0.608.
ari_tf_ctrl <- mclust::adjustedRandIndex(epi_seurat_object$tf_cluster,
                                         epi_seurat_object$ctrl_cluster)
ari_tf_ctrl

# So the transcription factors did WORSE than an equally sparse set of
# ordinary genes. That is the opposite of what I expected when I designed
# this, and it is the finding. I am reporting it as it came out rather than
# looking for a variant of the analysis that agrees with my hypothesis.

ari_results <- data.frame(
  comparison = c("TF-only vs reference",
                 "Expression-matched control vs reference",
                 "TF-only vs control"),
  n_genes = c(length(tf_genes), length(control_genes), NA),
  n_clusters = c(length(unique(epi_seurat_object$tf_cluster)),
                 length(unique(epi_seurat_object$ctrl_cluster)),
                 NA),
  ARI = c(ari_tf, ari_ctrl, ari_tf_ctrl)
)

ari_results

write_csv(ari_results, "results/tables/tf_ari_results.csv")


#### WHERE exactly did the transcription factors fail? ####

# The ARI is one number and one number does not tell me what went wrong.
# This cross-tabulation does. If each annotated cell state maps cleanly onto
# one TF cluster, transcription factors recovered it. If several states
# collapse into one TF cluster, they were not distinguishable by TFs alone.
#
# Reading my own output: TF cluster 0 absorbs 6,057 cells, and they are the
# hormone-sensing cells from ALL FOUR developmental stages plus the luminal
# progenitors. So transcription factors kept the LINEAGE (hormone-sensing
# versus alveolar versus basal) and lost the developmental STATE within that
# lineage. That is a much more specific and more interesting statement than
# "the ARI was low", and it is what I put on the slide.
tf_vs_reference <- table(epi_seurat_object$cluster_label,
                         epi_seurat_object$tf_cluster)

tf_vs_reference

write.csv(as.data.frame.matrix(tf_vs_reference),
          "results/tables/tf_clusters_vs_reference.csv")


#### UMAPs of the two alternative gene spaces ####

epi_seurat_object <- RunUMAP(epi_seurat_object,
                             reduction = "pca_tf",
                             dims = 1:20,
                             reduction.name = "umap_tf")

epi_seurat_object <- RunUMAP(epi_seurat_object,
                             reduction = "pca_ctrl",
                             dims = 1:20,
                             reduction.name = "umap_ctrl")

saveRDS(epi_seurat_object, "RObjects/04_tf.rds")


#### THE MAIN FIGURE ####

# Three embeddings side by side, with every cell keeping the colour of its
# REFERENCE annotation in all three panels. That is what makes the picture
# readable: if a colour stays as one blob, that cell state survived in this
# gene space; if a colour smears across the panel, it did not.
# The question is visual as well as numerical, and I wanted both.

p_ref <- DimPlot(epi_seurat_object,
                 reduction = "tsne",
                 group.by = "cluster_label") +
  NoLegend() +
  ggtitle("Full transcriptome (3,000 variable genes)")

p_tf <- DimPlot(epi_seurat_object,
                reduction = "umap_tf",
                group.by = "cluster_label") +
  NoLegend() +
  ggtitle(paste0("Transcription factors only (", length(tf_genes), " genes)"))

p_ctrl <- DimPlot(epi_seurat_object,
                  reduction = "umap_ctrl",
                  group.by = "cluster_label") +
  ggtitle(paste0("Expression-matched control (", length(control_genes), " genes)"))

p_main <- p_ref + p_tf + p_ctrl +
  plot_annotation(
    title = "Do transcription factors alone encode mammary epithelial cell identity?",
    subtitle = "All cells coloured by their full-transcriptome annotation")

p_main

ggsave("results/graphs/04_tf_vs_control_vs_reference.png",
       plot = p_main, width = 20, height = 6.5, dpi = 300)


#### Final save ####

saveRDS(epi_seurat_object, "RObjects/04_tf.rds")
