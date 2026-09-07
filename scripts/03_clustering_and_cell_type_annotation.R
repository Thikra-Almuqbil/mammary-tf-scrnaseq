# ==============================================================
# SCRIPT 3 OF 4: CLUSTERING AND CELL TYPE ANNOTATION
#
# From a map of all captured cells to fifteen named epithelial cell states.
# This is where the paper's Figures 1 and 2B to 2C get reproduced.
#
#   PART 1  Coarse clustering, only to find the non-epithelial cells
#   PART 2  Remove them, re-embed the epithelium        -> Fig 1b, Fig 1d
#   PART 3  Sweep the clustering parameters, score by silhouette
#   PART 4  Lock in the final clustering                -> Fig 1c
#   PART 5  OPTIONAL: does Leiden agree with Louvain?
#   PART 6  Marker panels, and two clusters I had to think about -> Fig 2b
#   PART 7  Label the fifteen clusters                  -> Fig 2c
#   PART 8  OPTIONAL and SLOW: ranked marker lists
#
# This is the longest script, so every part ends with a checkpoint and a
# "safe stopping point" note. Parts 5 and 8 can both be skipped entirely
# without affecting anything downstream.
#
# Course reference:
#   Demonstrations/08_Clustering.R
#   Demonstrations/09_Cluster_Marker_Genes.R
#   Demonstrations/06_DimensionalityReduction.R (t-SNE, UMAP)
#   Demonstrations/04_Preprocessing_And_QC.R (subset with cells =)
#
# Expected runtime: about 70 minutes for parts 1 to 4 and 6 to 7,
#                   plus 5 to 15 for part 5 and 10 to 20 for part 8.
# Checkpoints written: RObjects/03a_coarse_clusters.rds
#                      RObjects/03b_epithelial_pca.rds
#                      RObjects/03c_epithelial.rds
#                      RObjects/03d_sweep.rds
#                      RObjects/03e_clustered.rds
#                      RObjects/03f_annotated.rds
# ==============================================================


library(Seurat)
library(cluster)    # this is where silhouette() comes from, used in part 3
library(tidyverse)
library(patchwork)

theme_set(theme_classic())

options(future.globals.maxSize = 8 * 1024^3)

# Set your working directory to the ROOT of this repository before running.

dir.create("results/graphs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

qc_seurat_object <- readRDS("RObjects/02c_umap.rds")
qc_seurat_object



# ==============================================================
# PART 1: COARSE CLUSTERING, TO FIND THE NON-EPITHELIAL CELLS
#
# THIS IS NOT MY FINAL CLUSTERING. It is deliberately coarse. I only need
# enough resolution to isolate the immune, fibroblast and endothelial
# populations so I can remove them, following the paper:
#
#   "we flagged clusters that expressed clear markers of non-epithelial cells
#    as 'contaminating cells' and removed them from the downstream analysis.
#    Clusters C16 and C17 were tagged as immune cells based on the expression
#    of Cd74, Cd72 and Cd52, C18 as fibroblasts based on the expression of
#    collagens and fibronectin and C19 as endothelial cells based on the
#    expression of Eng, S1pr1 and Emcn"
#
# The careful clustering, with parameter tuning, happens in part 3 AFTER the
# removal, so that I only pay for that exploration once and so that it is
# tuned on the cells I actually care about.
# ==============================================================

# The RNA layers were already rejoined in script 2, so this is a harmless
# no-op. I left it in so the script also works on an object where RNA is
# still split.
qc_seurat_object <- JoinLayers(qc_seurat_object, assay = "RNA")

# Reminder to myself: printing the object shows the layers of the ACTIVE assay
# only, which is SCT. So "3 layers present: counts, data, scale.data" is about
# SCT, not RNA. This should be just "counts" and "data".
Layers(qc_seurat_object[["RNA"]])


#### The clustering ####

# Step 1: build a shared nearest neighbour graph. Each cell is connected to
# its k closest neighbours in PCA space, and the edge weight depends on how
# many neighbours two cells share. Clustering then means cutting this graph
# into communities.
# I use the same 20 PCs as the UMAP, so the clusters and the picture agree.
qc_seurat_object <- FindNeighbors(qc_seurat_object,
                                  reduction = "pca",
                                  k.param = 20,  # 20 is Seurat's default
                                  dims = 1:20)

Graphs(qc_seurat_object)

# Step 2: cut the graph. resolution controls granularity: higher means more,
# smaller clusters. The default is 0.8 and I use 0.5, lower, because here I
# only want broad populations and I do not want the immune cells split into
# five subtypes I would then have to identify one by one.
# algorithm = 1 is Louvain, Seurat's default, and it is fast.
qc_seurat_object <- FindClusters(qc_seurat_object,
                                 resolution = 0.5,
                                 algorithm = 1,
                                 cluster.name = "coarse_clusters")

# How many clusters and how big. I got 19.
table(qc_seurat_object$coarse_clusters)

# Which stages does each cluster draw from? A cluster made of one stage is
# probably a cell state; a cluster drawing from all four is more likely to be
# a cell type that exists throughout, which is what I expect the immune and
# endothelial cells to look like.
table(qc_seurat_object$coarse_clusters, qc_seurat_object$SampleGroup)

p_coarse <- DimPlot(qc_seurat_object,
                    reduction = "umap",
                    group.by = "coarse_clusters",
                    label = TRUE,
                    label.size = 5) +
  NoLegend() +
  labs(title = "Coarse clustering of all captured cells",
       subtitle = "Louvain, k = 20, resolution 0.5; 19 clusters from 25,802 cells")

p_coarse

ggsave("results/graphs/03_coarse_clusters_umap.png",
       plot = p_coarse, width = 8, height = 6, dpi = 300, bg = "white")


#### Which clusters are not epithelial? ####

# I plot from the RNA assay (log-normalised), not SCT, because several of
# these markers were among the 232 genes SCTransform could not scale across
# all 8 sample layers. Ptprc, Cdh5, Emcn and Eng are all in that list, and on
# the SCT assay they would simply be missing.
DefaultAssay(qc_seurat_object) <- "RNA"

# THE TRAP IN THIS PANEL, which I nearly walked into:
# Acta2 is a myoepithelial marker AS WELL AS a smooth muscle and pericyte
# marker. Myoepithelial cells are EPITHELIAL and are one of the largest
# populations in this dataset (the paper's C14 has 7,741 cells, nearly all
# from the lactation samples). Deleting them because they are Acta2-positive
# would have destroyed the entire lactation arm of my analysis.
# The thing that separates them is epithelial identity, not the contractile
# genes:
#   myoepithelial   Acta2+ AND Epcam+/Krt5+/Krt14+   -> keep
#   pericyte, other Acta2+ but Epcam-/Krt-           -> remove
contaminant_markers <- c(
  # epithelial, these clusters are KEPT
  "Epcam", "Krt8", "Krt18",          # luminal
  "Krt5", "Krt14", "Acta2",          # basal and myoepithelial
  # immune, REMOVE
  "Cd74", "Cd52", "Ptprc", "C1qa", "Lyz2",
  # fibroblast, REMOVE
  "Col1a1", "Col1a2", "Dcn", "Lum", "Fn1",
  # endothelial, REMOVE
  "Pecam1", "Eng", "Emcn", "Cdh5"
)

# Check every one is present before plotting, because DotPlot's error message
# for a missing gene is not obvious
contaminant_markers %in% rownames(qc_seurat_object[["RNA"]])

# The dot plot is the main output of this part. For every cluster it shows
# what fraction of its cells express each marker (dot size) and how strongly
# (colour). Reading across a row tells me what that cluster is. No statistical
# test is involved and none is needed: these are known markers and I am
# asking where they are expressed, not discovering anything.
p_dotplot <- DotPlot(qc_seurat_object,
                     features = contaminant_markers,
                     group.by = "coarse_clusters") +
  scale_colour_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Lineage markers identify immune, fibroblast and endothelial clusters",
       subtitle = "Dot size = percent of cells expressing; colour = mean expression")

p_dotplot

ggsave("results/graphs/03_contaminant_dotplot.png",
       plot = p_dotplot, width = 11, height = 6, dpi = 300, bg = "white")

# The same information as violins, which I find easier to read off when I want
# to be sure about one particular cluster
p_vln <- VlnPlot(qc_seurat_object,
                 features = c("Epcam", "Cd74", "Col1a1", "Pecam1"),
                 group.by = "coarse_clusters",
                 ncol = 2,
                 pt.size = 0)

p_vln

ggsave("results/graphs/03_contaminant_violins.png",
       plot = p_vln, width = 12, height = 8, dpi = 300, bg = "white")

DefaultAssay(qc_seurat_object) <- "SCT"


#### CHECKPOINT 1 ####

# Safe stopping point.
saveRDS(qc_seurat_object, "RObjects/03a_coarse_clusters.rds")



# ==============================================================
# PART 2: REMOVE THE NON-EPITHELIAL CELLS AND RE-EMBED
#
# What I identified from the dot plot above:
#   10  immune       Cd74, Cd52, Ptprc, C1qa, Lyz2 high; Epcam low    677 cells
#   12  fibroblast   Col1a1, Col1a2, Dcn, Lum, Fn1 high               393 cells
#   13  endothelial  Pecam1, Eng, Emcn, Cdh5 high                     365 cells
#   15  fibroblast   Col1a1, Col1a2, Dcn, Lum, Fn1 high               215 cells
#   16  endothelial  Pecam1, Eng, Emcn, Cdh5 high                     163 cells
#                                                            total   1813 cells
#
# TWO CLUSTERS I DECIDED TO KEEP, and why:
#
# Cluster 11 (403 cells) is carried forward but flagged as ambiguous. It is
# Epcam-negative and uniquely positive for the mural markers Des and Cspg4,
# which points at pericyte or vascular smooth muscle. But the paper explicitly
# RETAINED an equivalent cluster ("we retained C15 ... despite being positive
# for 2 out of 4 (Des and Cspg4) pericyte markers") because it matched the
# previously described Procr+ basal population. Procr itself does not settle
# it, because Procr is also high in my endothelial clusters 13 and 16, it is
# an endothelial surface protein (EPCR).
# I keep it because the costs are asymmetric. Keeping a contaminant is
# recoverable: I can label it and exclude it from interpretation. Deleting a
# real cell type is not recoverable, because I would never know it was there.
# I revisit the decision in part 6 with a proper marker panel.
#
# Clusters 0 and 7 (7,803 cells, 7,644 of them lactation) are Krt5+/Krt14+/
# Acta2+ MYOEPITHELIAL cells and must be kept despite low Epcam. Mouse mammary
# basal and myoepithelial cells are genuinely EpCAM-low. This is the mistake I
# was most worried about making.
# ==============================================================

# If you are resuming here:
# qc_seurat_object <- readRDS("RObjects/03a_coarse_clusters.rds")

contaminant_clusters <- c("10", "12", "13", "15", "16")

# I build a vector of the barcodes I want to KEEP and pass that to subset(),
# which is the idiom demonstration 04 uses. The "!" inverts the match, so this
# is every cell that is not in one of those five clusters.
keep_cells <- colnames(qc_seurat_object)[
  !(qc_seurat_object$coarse_clusters %in% contaminant_clusters)
]

# Kept, and removed. I got 23,989 kept and 1,813 removed.
length(keep_cells)
ncol(qc_seurat_object) - length(keep_cells)

epi_seurat_object <- subset(qc_seurat_object, cells = keep_cells)
epi_seurat_object

# The five clusters should now show 0 cells. This is my check that the subset
# did what I meant rather than something subtly different.
table(epi_seurat_object$coarse_clusters)

# Cells per stage after removal. The paper had 23,184 epithelial cells and I
# have 23,989. Mine is higher because I kept cluster 11 and because I did not
# go looking for their doublet cluster C20.
table(epi_seurat_object$SampleGroup)

table(epi_seurat_object$SampleName)


#### Re-running PCA on the cells I kept ####

# Why re-run at all: PC2 and PC5 of the old PCA were the immune and
# endothelial axes. Those axes no longer exist in this object, so the old
# components are not a sensible basis for describing epithelial cells.
#
# Why I do NOT re-run SCTransform: its residuals are computed per gene per
# cell from a model fitted across cells, and removing 7 percent of the cells
# shifts those fits negligibly. The variable feature list still contains some
# immune and stromal genes (C1qa, Cd74, Dcn, Fabp4), but among the retained
# cells those genes now have near-zero variance, so they contribute almost
# nothing to the new components. Re-running SCTransform would have cost
# another 8 minutes to change almost nothing.
DefaultAssay(epi_seurat_object) <- "SCT"

epi_seurat_object <- RunPCA(epi_seurat_object,
                            features = VariableFeatures(epi_seurat_object))

Stdev(epi_seurat_object, reduction = "pca")

p_elbow <- ElbowPlot(epi_seurat_object, ndims = 50) +
  labs(title = "Variance explained by principal components",
       subtitle = "Mammary epithelial cells, n = 23,989")

p_elbow

ggsave("results/graphs/03_elbow_epithelial.png",
       plot = p_elbow, width = 7, height = 5, dpi = 300, bg = "white")

# What are the new components made of? PC4 turned out to carry a dissociation
# stress signature, which becomes relevant in part 6 when one cluster splits
# along exactly that axis.
print(epi_seurat_object[["pca"]], dims = 1:5, nfeatures = 15)

# Save before the slow embedding steps, so a crash in t-SNE does not cost me
# the PCA as well. Safe stopping point.
saveRDS(epi_seurat_object, "RObjects/03b_epithelial_pca.rds")


#### Stop and check before the slow part ####

# I continue with dims = 1:20, the same as before. If the new elbow looked
# materially different from the previous one (which was sd 23.3, 18.6, 17.0,
# 16.6, then a break to 11.9) this is where I would pause and reconsider,
# rather than after spending 20 minutes on the embeddings.


#### t-SNE ####

# perplexity = 50 matches the paper's Methods, so my Fig 1b is comparable to
# theirs. Perplexity is roughly "how many neighbours each cell should try to
# stay near", so it sets how local the picture is.
#
# A trap I hit: demonstration 06 writes `seed = 123`, but the Seurat argument
# is `seed.use`. Written as `seed` it gets passed straight through to Rtsne,
# which has no such argument, so the seed is silently not set and the plot
# changes every run. It has to be seed.use.
epi_seurat_object <- RunTSNE(epi_seurat_object,
                             reduction = "pca",
                             dims = 1:20,
                             perplexity = 50,
                             seed.use = 123)

Reductions(epi_seurat_object)

# Save immediately, t-SNE is the slowest single step in the whole project
saveRDS(epi_seurat_object, "RObjects/03c_epithelial.rds")


#### UMAP ####

epi_seurat_object <- RunUMAP(epi_seurat_object,
                             reduction = "pca",
                             dims = 1:20)

Reductions(epi_seurat_object)


#### My version of the paper's Fig 1b ####

# t-SNE coloured by the four developmental stages. The paper used pink for NP,
# dark green for G, light green for L and purple for PI. I am not matching
# their colours, only the analytical content.
p_fig1b <- DimPlot(epi_seurat_object,
                   reduction = "tsne",
                   group.by = "SampleGroup") +
  labs(title = "Mammary epithelial cells segregate by developmental stage",
       subtitle = "t-SNE of 23,989 cells, perplexity 50, 20 principal components")

p_fig1b

ggsave("results/graphs/03_fig1b_tsne_by_stage.png",
       plot = p_fig1b, width = 8, height = 6, dpi = 300, bg = "white")

p_fig1b_split <- DimPlot(epi_seurat_object,
                         reduction = "tsne",
                         group.by = "SampleGroup",
                         split.by = "SampleGroup") +
  labs(title = "Stage-specific composition of the mammary epithelium",
       subtitle = "Nulliparous, gestation, lactation and post-involution occupy distinct regions")

p_fig1b_split

ggsave("results/graphs/03_fig1b_split.png",
       plot = p_fig1b_split, width = 14, height = 4, dpi = 300, bg = "white")

# The replicate check again, now on the epithelial cells. This is the figure
# I point at when someone asks whether the stage separation is a batch effect.
p_tsne_sample <- DimPlot(epi_seurat_object,
                         reduction = "tsne",
                         group.by = "SampleName") +
  labs(title = "Biological replicates co-localise within each developmental stage",
       subtitle = "Two independent mice per stage; separation between stages is not a batch effect")

p_tsne_sample

ggsave("results/graphs/03_tsne_by_sample.png",
       plot = p_tsne_sample, width = 8, height = 6, dpi = 300, bg = "white")


#### My version of the paper's Fig 1d ####

# Krt5 (basal) and Krt18 (luminal), from the log-normalised RNA assay so the
# values are comparable to the paper's log-transformed normalised counts.
DefaultAssay(epi_seurat_object) <- "RNA"

p_fig1d <- FeaturePlot(epi_seurat_object,
                       reduction = "tsne",
                       features = c("Krt5", "Krt18")) +
  plot_annotation(
    title = "Basal and luminal compartments are mutually exclusive",
    subtitle = "Log-normalised expression of Krt5 (basal) and Krt18 (luminal)")

p_fig1d

ggsave("results/graphs/03_fig1d_krt5_krt18.png",
       plot = p_fig1d, width = 11, height = 5, dpi = 300, bg = "white")

DefaultAssay(epi_seurat_object) <- "SCT"


#### UMAP of the same cells, for comparison ####

# I make this so I can check that the structure I am about to interpret is not
# an artefact of one particular embedding algorithm.
p_umap_stage <- DimPlot(epi_seurat_object,
                        reduction = "umap",
                        group.by = "SampleGroup") +
  labs(title = "UMAP embedding reproduces the stage structure seen by t-SNE",
       subtitle = "23,989 mammary epithelial cells, 20 principal components")

p_umap_stage

ggsave("results/graphs/03_umap_epithelial_by_stage.png",
       plot = p_umap_stage, width = 8, height = 6, dpi = 300, bg = "white")


#### CHECKPOINT 2 ####

# Safe stopping point.
saveRDS(epi_seurat_object, "RObjects/03c_epithelial.rds")



# ==============================================================
# PART 3: CHOOSING THE CLUSTERING PARAMETERS INSTEAD OF GUESSING
#
# Clustering has two knobs, k (how many neighbours in the graph) and
# resolution (how finely to cut it), and the number of clusters I get is
# entirely determined by what I set them to. I did not want to pick a value
# because it gave me the answer I liked, so I swept nine settings and scored
# them with silhouette width.
#
# What silhouette width is, in one paragraph, because I had to look it up:
# for one cell, take its average distance to the other cells in its own
# cluster (call it a) and its average distance to the cells of the nearest
# other cluster (call it b). The silhouette is (b - a) / max(a, b). It is
# near 1 if the cell is much closer to its own cluster than to any other,
# near 0 if it sits on a boundary, and negative if it would fit better
# somewhere else. Averaging over cells gives one number per clustering.
#
# IMPORTANT CAVEAT that I state in the presentation: silhouette width rewards
# compact, well separated clusters, so it systematically prefers COARSER
# solutions. A peak at 15 clusters does not prove 15 is biologically correct.
#
# TWO DEPARTURES FROM THE COURSE, both forced by dataset size:
#
#  1. I compute silhouette on a random subsample of 5,000 cells.
#     Demonstration 08 runs dist() on the full embedding, which is fine for
#     its 2,000-cell demo object (a 16 MB distance matrix). My epithelial
#     object has about 24,000 cells, where the distance matrix is roughly
#     2.3 GB and silhouette() is O(n^2) on top of that.
#
#  2. I sweep with Louvain (algorithm = 1), not Leiden.
#     Demonstration 08's final choice was Leiden. Leiden in Seurat can be
#     slow and memory hungry on a graph this size. I sweep with Louvain, then
#     run Leiden ONCE at the chosen parameters in part 5 and compare.
# ==============================================================

# If you are resuming here:
# epi_seurat_object <- readRDS("RObjects/03c_epithelial.rds")

DefaultAssay(epi_seurat_object) <- "SCT"


#### Sweep 1: k, holding resolution at Seurat's default of 0.8 ####

# Seurat's default k.param is 20 no matter how big the dataset is, and it was
# tuned on much smaller data. With about 24,000 cells a larger k gives a more
# robust graph, because each cell's neighbourhood is estimated from more
# evidence. This loop is copied from demonstration 08's "Testing multiple k".
k_values <- c(15, 20, 25, 30)

for (k in k_values) {
  epi_seurat_object <- FindNeighbors(epi_seurat_object,
                                     reduction = "pca",
                                     k.param = k,
                                     dims = 1:20)

  epi_seurat_object <- FindClusters(epi_seurat_object,
                                    resolution = 0.8,
                                    algorithm = 1,
                                    cluster.name = paste0("Louvain_k", k,
                                                          "_res0.8"))
}

colnames(epi_seurat_object[[]])


#### Sweep 2: resolution, holding k at 25 ####

# k = 25 is a middle value from the sweep above. FindClusters accepts a vector
# of resolutions, so the graph is only built once and then cut five times,
# which is much faster than five separate calls.
epi_seurat_object <- FindNeighbors(epi_seurat_object,
                                   reduction = "pca",
                                   k.param = 25,
                                   dims = 1:20)

epi_seurat_object <- FindClusters(epi_seurat_object,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.2),
                                  algorithm = 1)

# These land in columns named SCT_snn_res.<value>, which is why part 4 copies
# the chosen one into a clearly named column: nothing in the name records that
# these were run at k = 25.
colnames(epi_seurat_object[[]])

# Save before the silhouette calculations, which are the part most likely to
# run out of memory
saveRDS(epi_seurat_object, "RObjects/03d_sweep.rds")


#### How many clusters did each setting give? ####

table(epi_seurat_object$Louvain_k15_res0.8)
table(epi_seurat_object$Louvain_k20_res0.8)
table(epi_seurat_object$Louvain_k25_res0.8)
table(epi_seurat_object$Louvain_k30_res0.8)

table(epi_seurat_object$SCT_snn_res.0.4)
table(epi_seurat_object$SCT_snn_res.0.6)
table(epi_seurat_object$SCT_snn_res.0.8)
table(epi_seurat_object$SCT_snn_res.1)
table(epi_seurat_object$SCT_snn_res.1.2)


#### Silhouette widths on a 5,000-cell subsample ####

# set.seed makes the subsample reproducible, so re-running gives the same
# scores and I am not choosing parameters off a lottery.
set.seed(123)
sil_cells <- sample(colnames(epi_seurat_object), 5000)

# The PCA coordinates for those cells, using the same 20 dimensions the
# clustering used. Scoring in a different space from the one I clustered in
# would not be a fair test.
pca_sub <- Embeddings(epi_seurat_object, reduction = "pca")[sil_cells, 1:20]

# The distance matrix. On 5,000 cells this is about 95 MB, which is fine.
dist_sub <- dist(pca_sub)

# The two helper functions are lifted straight from demonstration 08
calc_mean_sil <- function(clusters, dist_matrix) {
  silhouette_width <- silhouette(clusters, dist_matrix)
  mean(silhouette_width[, "sil_width"])
}

calc_num_clusters <- function(clusters) {
  length(unique(clusters))
}

epi_meta <- epi_seurat_object[[]][sil_cells, ]

cluster_cols <- c("Louvain_k15_res0.8",
                  "Louvain_k20_res0.8",
                  "Louvain_k25_res0.8",
                  "Louvain_k30_res0.8",
                  "SCT_snn_res.0.4",
                  "SCT_snn_res.0.6",
                  "SCT_snn_res.0.8",
                  "SCT_snn_res.1",
                  "SCT_snn_res.1.2")

# For each setting, the number of clusters and the mean silhouette width.
# The as.integer(as.character(...)) is not decoration: the cluster columns are
# factors, and as.integer() on a factor returns the level INDEX, not the
# label, so without the round trip through character the cluster identities
# would be silently wrong.
silhouette_stats <- lapply(cluster_cols, function(this_col) {

  clusters <- epi_meta[, this_col]
  clusters <- as.integer(as.character(clusters))

  tibble(setting = this_col,
         num_clusters = calc_num_clusters(clusters),
         mean_silhouette_width = calc_mean_sil(clusters, dist_sub))
}) %>%
  bind_rows()

# My result: SCT_snn_res.0.4 (that is k = 25, resolution 0.4) scored 0.212,
# and the runner-up scored 0.160. That is about a 30 percent gap, which is
# not a close call, so I did not agonise over it.
print(silhouette_stats, n = Inf)

p_sil <- silhouette_stats %>%
  ggplot(aes(x = num_clusters, y = mean_silhouette_width)) +
  geom_point(size = 3) +
  geom_text(aes(label = setting), hjust = -0.1, size = 3) +
  expand_limits(x = max(silhouette_stats$num_clusters) + 8) +
  labs(x = "Number of clusters",
       y = "Mean silhouette width (5000-cell subsample)",
       title = "Clustering parameter sweep")

p_sil

ggsave("results/graphs/03_silhouette_sweep.png",
       plot = p_sil, width = 9, height = 6, dpi = 300, bg = "white")


#### Two of the options drawn on the t-SNE ####

# Numbers are not enough. I wanted to see what over-clustering looks like on
# this data, and at resolution 1.2 the myoepithelial mass fragments into
# several pieces that the paper reports as one cluster.
p_k25 <- DimPlot(epi_seurat_object,
                 reduction = "tsne",
                 group.by = "SCT_snn_res.0.8",
                 label = TRUE, label.size = 4) +
  NoLegend() + ggtitle("k = 25, resolution 0.8")

p_k25_r12 <- DimPlot(epi_seurat_object,
                     reduction = "tsne",
                     group.by = "SCT_snn_res.1.2",
                     label = TRUE, label.size = 4) +
  NoLegend() + ggtitle("k = 25, resolution 1.2")

p_compare <- p_k25 + p_k25_r12

p_compare

ggsave("results/graphs/03_cluster_options.png",
       plot = p_compare, width = 14, height = 6, dpi = 300, bg = "white")


#### CHECKPOINT 3 ####

# Safe stopping point.
saveRDS(epi_seurat_object, "RObjects/03d_sweep.rds")



# ==============================================================
# PART 4: LOCKING IN THE FINAL CLUSTERING
#
# Chosen from the sweep above: Louvain, k = 25, resolution 0.4.
#
# Why:
#   - highest mean silhouette width of the nine settings tested, 0.212 against
#     0.160 for the runner-up, roughly a 30 percent gap
#   - it gives exactly 15 clusters, which is the number Bach et al. report.
#     This was NOT my selection criterion. I chose on silhouette width and 15
#     fell out of it, so the agreement is independent corroboration rather
#     than the target. I want to be clear about the direction of that, because
#     picking the resolution that reproduces the paper's cluster count and
#     then presenting the match as a result would be circular.
#   - at higher resolutions the myoepithelial mass fragments into four or five
#     pieces, whereas the paper reports it as a single cluster (C14)
#
# Caveats I record rather than hide:
#   - silhouette width prefers coarser solutions, so its peak at 15 does not
#     prove 15 is biologically right.
#   - the paper reached 15 clusters by a two-step route (coarse clustering,
#     then hierarchical subclustering within each), so their clusters are
#     finer grained than any single-pass Louvain can be. If part 6 had found a
#     cluster with mixed markers, resolution 0.6 (17 clusters) was my fallback.
# ==============================================================

# SCT_snn_res.0.4 was produced at k = 25 in part 3, but nothing in that column
# name records the k. I copy it to a self-documenting name so that nothing
# later depends on me remembering which sweep it came from.
epi_seurat_object$final_cluster <- epi_seurat_object$SCT_snn_res.0.4

# 15 clusters, sizes from 6,418 down to 61
table(epi_seurat_object$final_cluster)

# Which stages does each cluster come from? The paper's Table 1 shows clusters
# that are almost entirely single-stage, and mine came out the same way.
table(epi_seurat_object$final_cluster, epi_seurat_object$SampleGroup)

# And per sample. This is a real check, not a formality: a cluster present in
# only ONE of the two mice of a stage would be suspicious, because it would
# suggest a per-animal or per-run artefact rather than a genuine cell state.
table(epi_seurat_object$final_cluster, epi_seurat_object$SampleName)


#### My version of the paper's Fig 1c ####

p_fig1c <- DimPlot(epi_seurat_object,
                   reduction = "tsne",
                   group.by = "final_cluster",
                   label = TRUE,
                   label.size = 5) +
  NoLegend() +
  labs(title = "Fifteen transcriptional states in the mammary epithelium",
       subtitle = "Louvain clustering, k = 25, resolution 0.4; 23,989 cells")

p_fig1c

ggsave("results/graphs/03_fig1c_clusters_tsne.png",
       plot = p_fig1c, width = 8, height = 6, dpi = 300, bg = "white")

# The same clusters on the UMAP. If the cluster boundaries agreed with the
# t-SNE only, I would worry that I was reading structure out of one particular
# embedding rather than out of the data.
p_fig1c_umap <- DimPlot(epi_seurat_object,
                        reduction = "umap",
                        group.by = "final_cluster",
                        label = TRUE,
                        label.size = 5) +
  NoLegend() +
  labs(title = "Cluster assignment is consistent between t-SNE and UMAP",
       subtitle = "Louvain, k = 25, resolution 0.4")

p_fig1c_umap

ggsave("results/graphs/03_clusters_umap.png",
       plot = p_fig1c_umap, width = 8, height = 6, dpi = 300, bg = "white")


#### Cluster composition by stage ####

# This is the demonstration 07 barplot idiom applied to my clusters.
# position = "fill" makes every bar the same height so I am reading
# proportions, not cluster sizes.
# This one figure is the empirical justification for two decisions: not batch
# correcting, and being careful about what "cell state" means here. Almost
# every cluster is one colour, which means cell identity and developmental
# stage are confounded BY DESIGN in this dataset. That is not a flaw in my
# analysis, it is what the experiment is, and it limits what any
# stage-versus-identity claim can mean.
p_composition <- data.frame(
    Cluster = epi_seurat_object$final_cluster,
    Stage = epi_seurat_object$SampleGroup
  ) %>%
  ggplot(aes(x = Cluster)) +
  geom_bar(aes(fill = Stage), position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Nearly every cluster is restricted to a single developmental stage",
       subtitle = "Cell type and developmental stage are confounded by design in this dataset",
       y = "Proportion of cells")

p_composition

ggsave("results/graphs/03_cluster_composition_by_stage.png",
       plot = p_composition, width = 9, height = 6, dpi = 300, bg = "white")


#### CHECKPOINT 4 ####

# Safe stopping point. This is the object part 6 and part 7 read.
saveRDS(epi_seurat_object, "RObjects/03e_clustered.rds")



# ==============================================================
# PART 5: OPTIONAL. DOES THE ANSWER DEPEND ON THE ALGORITHM?
#
# SKIP THIS PART ENTIRELY IF THE PACKAGE WILL NOT INSTALL.
# Reporting Louvain alone is a perfectly defensible cut and it is Seurat's own
# default algorithm. Nothing downstream depends on this part.
#
# PREREQUISITE, run once:
#   install.packages("leidenbase")
# If that reports "only available in source form, may need compilation", stop
# here and move to part 6.
#
# IF FindClusters HANGS FOR MORE THAN ABOUT 10 MINUTES, STOP IT AND SKIP.
# ==============================================================

# This should print TRUE. Note the print(): requireNamespace() returns its
# answer INVISIBLY, so on its own it appears to do nothing at all, which
# confused me for a few minutes.
print(requireNamespace("leidenbase", quietly = TRUE))

# Same k, same dims, same resolution. The ONLY thing that changes is the
# algorithm, otherwise the comparison would not isolate the algorithm.
epi_seurat_object <- FindNeighbors(epi_seurat_object,
                                   reduction = "pca",
                                   k.param = 25,
                                   dims = 1:20)

# algorithm = 4 is Leiden. I set random.seed because Leiden warns without it.
epi_seurat_object <- FindClusters(epi_seurat_object,
                                  resolution = 0.4,
                                  algorithm = 4,
                                  random.seed = 123,
                                  cluster.name = "Leiden_k25_res0.4")

table(epi_seurat_object$Leiden_k25_res0.4)

# The adjusted Rand index compares two partitions of the SAME cells. It works
# pairwise: for every pair of cells it asks whether the two clusterings agree
# about whether that pair belongs together, and then corrects for how much
# agreement you would expect by chance alone.
#   1 = identical partitions
#   0 = agreement no better than chance
# It does not care what the clusters are called, which is exactly what I need
# here, because Louvain's cluster 3 and Leiden's cluster 3 are unrelated
# labels.
#
# A FORCED SUBSTITUTION: demonstration 07 uses bluster::pairwiseRand(mode =
# "index") for this. Bioconductor 3.22 ships bluster for Windows in source
# form only and building it needs Rtools, which I did not have. I use
# mclust::adjustedRandIndex() instead, which computes the identical statistic.
# I call it with :: rather than attaching mclust, because attaching mclust
# masks dplyr::count and purrr::map and breaks later code.
#
# I got 0.876, which I read as: the two algorithms found the same structure
# and mostly disagree about where a few boundary cells go. So my results are
# not an artefact of choosing Louvain.
mclust::adjustedRandIndex(epi_seurat_object$final_cluster,
                          epi_seurat_object$Leiden_k25_res0.4)

# Which Louvain cluster maps onto which Leiden cluster? A near-diagonal table
# means the same structure, numbered differently.
table(epi_seurat_object$final_cluster,
      epi_seurat_object$Leiden_k25_res0.4)

p_louvain <- DimPlot(epi_seurat_object,
                     reduction = "tsne",
                     group.by = "final_cluster",
                     label = TRUE, label.size = 4) +
  NoLegend() +
  ggtitle("Louvain (algorithm 1)")

p_leiden <- DimPlot(epi_seurat_object,
                    reduction = "tsne",
                    group.by = "Leiden_k25_res0.4",
                    label = TRUE, label.size = 4) +
  NoLegend() +
  ggtitle("Leiden (algorithm 4)")

p_algo <- p_louvain + p_leiden +
  plot_annotation(
    title = "Louvain and Leiden recover the same population structure",
    subtitle = "Both at k = 25, resolution 0.4, 20 principal components")

p_algo

ggsave("results/graphs/03_louvain_vs_leiden.png",
       plot = p_algo, width = 14, height = 6, dpi = 300, bg = "white")

# This part only ADDS a column, it does not touch final_cluster, so saving
# back over the checkpoint is safe.
saveRDS(epi_seurat_object, "RObjects/03e_clustered.rds")



# ==============================================================
# PART 6: MARKER PANELS, AND THE TWO CLUSTERS I HAD TO THINK ABOUT
#
# The paper's Table 1 assignments, which are what I annotate against:
#   Hormone sensing progenitor  (Hsp)  Esr1, Prlr, Pgr, S100a6, Cited1
#   Hormone sensing diff.       (Hsd)  Esr1, Prlr, Pgr, S100a6, Cited1
#   Luminal progenitor          (Lp)   Aldh1a3, Cd14, Kit
#   Alveolar differentiated     (Avd)  Wap, Csn2, Glycam1, Lalba
#   Alveolar progenitor         (Avp)  Wap, Csn2, Glycam1, Lalba, Aldh1a3, Cd14, Kit
#   Basal                       (Bsl)  Krt4, Krt14, Pdpn, Etv5, Acta2
#   Myoepithelial               (Myo)  Oxtr, Acta2, Krt4, Krt14
#   Procr+ basal                (Prc)  Procr, Igfbp4, Gng11, Zeb2
#
# This part is fast, about a minute, and it is what actually annotates the
# clusters. The slow ranked-marker version is part 8 and is a bonus.
# ==============================================================

# If you are resuming here:
# epi_seurat_object <- readRDS("RObjects/03e_clustered.rds")

# I work from the RNA assay (log-normalised) for the rest of this script, for
# two reasons:
#  1. some of the paper's markers, Wap and Glycam1 in particular, were among
#     the 232 features SCTransform could not scale across all 8 sample layers,
#     so on the SCT assay they would not be there.
#  2. the SCT route for marker detection fails on this object entirely. See
#     the long note in part 8.
DefaultAssay(epi_seurat_object) <- "RNA"

Layers(epi_seurat_object[["RNA"]])


#### The paper's Table 1 marker panel across my clusters ####

# THIS IS THE PLOT THAT ANNOTATES THE CLUSTERS. For each cluster it shows what
# fraction of its cells express each known marker (dot size) and how strongly
# (colour). Reading down a column tells me which clusters are hormone-sensing,
# which are alveolar, which are basal, and so on.
#
# No statistical test is involved and none is needed. These are published
# markers and I am asking where they are expressed. A test would be answering
# a different question, namely "which genes distinguish this cluster", which
# is part 8 and is a bonus rather than the basis of the annotation.

table1_markers <- c(
  # hormone sensing
  "Esr1", "Prlr", "Pgr", "Cited1", "S100a6", "Areg", "Foxa1", "Ly6a",
  # luminal progenitor
  "Aldh1a3", "Cd14", "Kit", "Elf5",
  # alveolar
  "Wap", "Csn2", "Glycam1", "Lalba", "Olah", "Thrsp",
  # basal and myoepithelial
  "Krt5", "Krt14", "Krt4", "Acta2", "Oxtr", "Pdpn", "Etv5", "Krt15",
  # Procr+ basal
  "Procr", "Igfbp4", "Gng11", "Zeb2", "Notch3",
  # general luminal
  "Krt8", "Krt18", "Epcam"
)

# Anything missing from the RNA assay? Should print character(0).
table1_markers[!(table1_markers %in% rownames(epi_seurat_object[["RNA"]]))]

p_table1 <- DotPlot(epi_seurat_object,
                    features = table1_markers,
                    group.by = "final_cluster") +
  scale_colour_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Canonical marker expression assigns an identity to every cluster",
       subtitle = "Dot size = percent of cells expressing; colour = mean log-normalised expression")

p_table1

ggsave("results/graphs/03_table1_markers_dotplot.png",
       plot = p_table1, width = 14, height = 7, dpi = 300, bg = "white")


#### The genes shown in the paper's Fig 2b ####

p_fig2b <- FeaturePlot(epi_seurat_object,
                       reduction = "tsne",
                       features = c("Esr1", "Prlr", "Csn2",
                                    "Acta2", "Oxtr", "Pdpn",
                                    "Procr", "Zeb2", "Aldh1a3"),
                       ncol = 3) +
  plot_annotation(
    title = "Canonical lineage markers occupy discrete regions of the embedding",
    subtitle = "Log-normalised expression of nine cell-state markers across 23,989 mammary epithelial cells")

p_fig2b

ggsave("results/graphs/03_fig2b_markers.png",
       plot = p_fig2b, width = 13, height = 12, dpi = 300, bg = "white")


#### Settling the cluster 11 question I left open in part 1 ####

# Coarse cluster 11 was carried forward flagged as ambiguous: Epcam-negative,
# uniquely Des+ and Cspg4+, positive for mural markers, but the paper
# explicitly retained an equivalent cluster as Procr+ basal cells.
# First: which final cluster did it become?
table(epi_seurat_object$coarse_clusters, epi_seurat_object$final_cluster)

p_mural <- VlnPlot(epi_seurat_object,
                   features = c("Des", "Cspg4", "Rgs5", "Cox4i2",
                                "Epcam", "Krt14"),
                   group.by = "final_cluster",
                   ncol = 3,
                   pt.size = 0)

p_mural

ggsave("results/graphs/03_mural_markers.png",
       plot = p_mural, width = 14, height = 8, dpi = 300, bg = "white")


#### Is cluster 14 a stress artefact or a doublet cluster? ####

# Cluster 14 (61 cells) is the only one that is not stage-restricted:
# NP 10, G 1, L 20, PI 30. A small population appearing across all four
# stages is the signature of something technical rather than biological,
# because nothing in mouse mammary biology is present at 0.25 percent
# frequency in every physiological state.
# PC4 of the epithelial PCA carried a dissociation stress signature, so:
#   - high Fos, Jun, Egr1, Hspa1a with no coherent cell-type markers means
#     dissociation stress
#   - co-expression of markers from two different lineages, for example both
#     Csn2 and Krt14, means doublets
p_stress <- VlnPlot(epi_seurat_object,
                    features = c("Fos", "Jun", "Egr1",
                                 "Hspa1a", "Csn2", "Krt14"),
                    group.by = "final_cluster",
                    ncol = 3,
                    pt.size = 0)

p_stress

ggsave("results/graphs/03_stress_markers.png",
       plot = p_stress, width = 14, height = 8, dpi = 300, bg = "white")



# ==============================================================
# PART 7: LABELLING THE FIFTEEN CLUSTERS
#
# Identities assigned from the Table 1 dot plot in part 6, cross-referenced
# against the cluster-by-stage table from part 4.
#
#  cluster  n     stage  defining markers                     identity
#  0        6418  L      Krt5/Krt14/Acta2/Oxtr, Epcam-low     Myoepithelial
#  1        3156  PI     Esr1/Prlr/Pgr/Cited1/Areg/Foxa1      Hormone-sensing
#  2        2311  NP     same, strongest in the dataset       Hormone-sensing
#  3        2249  PI     Aldh1a3/Cd14/Kit/Elf5                Luminal progenitor
#  4        2198  G      Elf5/Wap/Csn2/Glycam1/Olah/Thrsp     Alveolar diff.
#  5        1846  G      Elf5/Csn2/Glycam1/Kit                Alveolar
#  6        1394  G      Krt5/Krt14/Acta2/Pdpn/Krt15          Basal
#  7        1272  L      as c0 PLUS Fos/Jun/Egr1/Hspa1a       Myoepithelial, stressed
#  8        1121  NP     Esr1/Prlr/Cited1/Ly6a + Aldh1a3      Hormone-sensing prog.
#  9        932   NP     Krt5/Krt14/Acta2/Pdpn                Basal
#  10       419   L      Procr/Igfbp4/Gng11/Zeb2/Notch3       Procr+ basal
#  11       342   L      Wap/Csn2/Lalba + Aldh1a3/Cd14/Kit    Alveolar progenitor
#  12       159   L      Krt5/Krt14/Acta2/Oxtr + some Csn2    Myoepithelial
#  13       111   G      Esr1/Prlr/Pgr/Cited1                 Hormone-sensing
#  14       61    mixed  Epcam+ with BOTH Csn2 and Krt14      Doublets (likely)
#
# TWO OF THESE NEED A CAVEAT STATED OUT LOUD, and I say both in the talk:
#
#  c10 (Prc). Expresses all four of Table 1's Procr+ basal markers strongly,
#      but is ALSO Epcam-negative and positive for four mural markers (Des,
#      Cspg4, Rgs5, Cox4i2). This is the same ambiguity the paper describes,
#      "we retained C15 ... despite being positive for 2 out of 4 (Des and
#      Cspg4) pericyte markers", and I resolve it the same way, by retaining
#      it and flagging it rather than pretending it is clean.
#
#  c7 (Myo, stressed). Clusters 0 and 7 together are 7,690 cells against the
#      paper's single C14 (Myo) at 7,741. Cluster 7 is the standout for the
#      immediate-early genes Fos, Jun and Egr1 and the heat-shock gene Hspa1a,
#      which is the classic enzymatic-dissociation signature and exactly what
#      PC4 of the epithelial PCA predicted. So my clustering split ONE
#      biological population along a TECHNICAL axis. That is a real limitation
#      of my analysis and I report it as one.
# ==============================================================

#### Applying the labels ####

# case_match() maps each cluster number to a label. It is the dplyr equivalent
# of the recoding block in demonstration 09.
# I keep the cluster NUMBER inside the label, for example "Myo (c0)", so that
# every figure can be traced back to the numbered clustering and nobody has to
# take my word for which cluster became which cell type.
epi_seurat_object$cluster_label <- case_match(
  as.character(epi_seurat_object$final_cluster),
  "0"  ~ "Myo (c0)",
  "1"  ~ "Hsd-PI (c1)",
  "2"  ~ "Hsd-NP (c2)",
  "3"  ~ "Lp-PI (c3)",
  "4"  ~ "Avd-G (c4)",
  "5"  ~ "Avd-G (c5)",
  "6"  ~ "Bsl-G (c6)",
  "7"  ~ "Myo stressed (c7)",
  "8"  ~ "Hsp-NP (c8)",
  "9"  ~ "Bsl-NP (c9)",
  "10" ~ "Prc (c10)",
  "11" ~ "Avp-L (c11)",
  "12" ~ "Myo (c12)",
  "13" ~ "Hsd-G (c13)",
  "14" ~ "Doublets (c14)"
)

# I order the labels biologically rather than alphabetically: hormone-sensing,
# luminal progenitor, alveolar, basal, myoepithelial, Procr+ basal, then the
# artefact cluster. Every figure downstream inherits this order, which is what
# makes the heatmap readable.
label_order <- c("Hsp-NP (c8)", "Hsd-NP (c2)", "Hsd-G (c13)", "Hsd-PI (c1)",
                 "Lp-PI (c3)",
                 "Avp-L (c11)", "Avd-G (c4)", "Avd-G (c5)",
                 "Bsl-NP (c9)", "Bsl-G (c6)",
                 "Myo (c0)", "Myo stressed (c7)", "Myo (c12)",
                 "Prc (c10)",
                 "Doublets (c14)")

epi_seurat_object$cluster_label <- factor(epi_seurat_object$cluster_label,
                                          levels = label_order)

# Set the labels as the active identity, as demonstration 09 does
Idents(epi_seurat_object) <- "cluster_label"

# Every cluster should have a label and none should be NA. An NA here would
# mean I mistyped a cluster number in case_match above, which is easy to do
# and silent.
table(epi_seurat_object$cluster_label)

# Labels against developmental stage. This is my annotated version of the
# paper's Table 1.
table(epi_seurat_object$cluster_label, epi_seurat_object$SampleGroup)

annotation_table <- as.data.frame.matrix(
  table(epi_seurat_object$cluster_label, epi_seurat_object$SampleGroup)
)
annotation_table$total <- rowSums(annotation_table)
annotation_table

write.csv(annotation_table,
          "results/tables/annotation_table_clusters_by_stage.csv")


#### Annotated version of the paper's Fig 1c ####

# repel = TRUE stops the labels sitting on top of each other, which they do
# badly with 15 clusters.
p_annotated <- DimPlot(epi_seurat_object,
                       reduction = "tsne",
                       group.by = "cluster_label",
                       label = TRUE,
                       label.size = 3.5,
                       repel = TRUE) +
  NoLegend() +
  labs(title = "Fifteen annotated cell states in the mouse mammary epithelium",
       subtitle = "23,989 cells; identities assigned from canonical lineage markers (Bach et al., 2017)")

p_annotated

ggsave("results/graphs/03_annotated_clusters_tsne.png",
       plot = p_annotated, width = 9, height = 7, dpi = 300, bg = "white")

# A version with the legend instead, for a slide where in-plot labels are
# too crowded to read from the back of a room
p_annotated_legend <- DimPlot(epi_seurat_object,
                              reduction = "tsne",
                              group.by = "cluster_label") +
  labs(title = "Fifteen annotated cell states in the mouse mammary epithelium",
       subtitle = "23,989 cells across four developmental stages")

p_annotated_legend

ggsave("results/graphs/03_annotated_clusters_legend.png",
       plot = p_annotated_legend, width = 11, height = 7, dpi = 300, bg = "white")


#### The paper's Fig 2c: the marker heatmap ####

# The paper's Fig 2c is a heatmap of key markers across cells, and their
# legend says "for visualisation purposes only 100 randomly selected cells
# were shown for large clusters". I do the same with subset(downsample = 100).
#
# This is a PLOTTING decision, not an analytical one. Without it, cluster 0
# with 6,418 cells would occupy most of the width and the 61-cell cluster
# would be a single invisible line, so the picture would be about cluster
# sizes rather than about marker expression. Nothing is being tested here, so
# downsampling costs no statistical power.

# set.seed makes the random selection reproducible, so the figure does not
# change every time I regenerate it for the slides.
set.seed(123)
epi_downsampled <- subset(epi_seurat_object, downsample = 100)

ncol(epi_downsampled)
table(epi_downsampled$cluster_label)

# DoHeatmap plots the scale.data layer, so the marker genes have to be scaled
# first or the heatmap comes out empty. I scale only these 34 genes rather
# than the whole matrix, which takes seconds instead of minutes.
# Scaling here means centring and dividing by the standard deviation per gene,
# so the colour shows "high or low FOR THIS GENE", which is what makes genes
# with very different absolute levels comparable in one picture.
epi_downsampled <- ScaleData(epi_downsampled, features = table1_markers)

p_fig2c <- DoHeatmap(epi_downsampled,
                     features = table1_markers,
                     group.by = "cluster_label",
                     size = 3,
                     angle = 45) +
  scale_fill_viridis_c() +
  labs(title = "Lineage-restricted marker expression defines fifteen mammary epithelial states",
       subtitle = "34 canonical markers, scaled expression; up to 100 randomly selected cells per cluster")

p_fig2c

ggsave("results/graphs/03_fig2c_marker_heatmap.png",
       plot = p_fig2c, width = 14, height = 9, dpi = 300, bg = "white")


#### Cluster composition by stage, with the labels ####

p_composition_labelled <- data.frame(
    Cluster = epi_seurat_object$cluster_label,
    Stage = epi_seurat_object$SampleGroup
  ) %>%
  ggplot(aes(x = Cluster)) +
  geom_bar(aes(fill = Stage), position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Each cell state is confined to one or two developmental stages",
       subtitle = "Cell identity and developmental stage are confounded by the experimental design",
       y = "Proportion of cells")

p_composition_labelled

ggsave("results/graphs/03_annotated_composition.png",
       plot = p_composition_labelled, width = 10, height = 6, dpi = 300, bg = "white")


#### CHECKPOINT 5: END OF THE REQUIRED PART OF SCRIPT 3 ####

# This is the object script 4 reads. Part 8 below is optional.
saveRDS(epi_seurat_object, "RObjects/03f_annotated.rds")



# ==============================================================
# PART 8: OPTIONAL AND SLOW. RANKED MARKER LISTS.
#
# SKIP THIS IF YOU ARE SHORT OF TIME. Nothing depends on it. The annotation
# in part 7 came from the published marker panel, not from these tests. I
# restructured the script this way after discovering that FindAllMarkers was
# going to take about 11 hours on my machine.
#
# WHY THIS IS SLOW, AND WHY IT IS ON THE RNA ASSAY.
#
# (a) The assay. The course runs PrepSCTFindMarkers() then FindAllMarkers() on
#     the SCT assay. That route FAILS here, and it fails SILENTLY, which cost
#     me an afternoon. SCTransform fitted one model per sample layer.
#     PrepSCTFindMarkers() compares each model's STORED median UMI count and
#     returns early saying "Minimum UMI unchanged", while FindMarkers()
#     compares the OBSERVED median UMIs and then refuses every single test
#     with "multiple models with unequal library sizes". The observed medians
#     stopped matching the stored ones the moment I removed 1,813 cells in
#     part 2. The result is an EMPTY table plus warnings, not an error, so
#     nothing tells you it went wrong except that the answer is zero rows.
#
# (b) The speed. Seurat v5's Wilcoxon test uses the presto package if it is
#     installed, and otherwise loops base R's wilcox.test one gene at a time.
#     presto is not on CRAN, r-universe has no build for R 4.5, and building
#     from source needs Rtools, which I did not have. Measured with the slow
#     fallback: 41 minutes for cluster 0 alone, so roughly 11 hours for all 15.
#
# The three arguments below make it tractable, and none of them costs me
# anything I keep:
#
#   max.cells.per.ident = 200
#     Downsamples each cluster before testing. Markers are large effects and
#     200 cells detects them comfortably. Note that at n = 24,000 the p-values
#     are meaningless anyway: Wilcoxon returns p around 1e-300 for a log2 fold
#     change of 0.1, because with enough cells everything is significant. My
#     selection below is by EFFECT SIZE, not by p-value, which is the right
#     way round. Clusters smaller than 200 (mine are 159, 111 and 61) are
#     unaffected.
#
#   min.pct = 0.25, logfc.threshold = 0.5
#     Skip testing genes that could not pass the effect-size filter I apply
#     afterwards. This is not an approximation, it is declining to compute
#     p-values I was going to throw away.
#
# Expected runtime with these settings: roughly 10 to 20 minutes.
# ==============================================================

set.seed(123)

markers_all <- FindAllMarkers(epi_seurat_object,
                              assay = "RNA",
                              group.by = "final_cluster",
                              max.cells.per.ident = 200,
                              min.pct = 0.25,
                              logfc.threshold = 0.5)

# THE CHECK THAT MATTERS, given how the SCT route failed. If this is 0,
# something went wrong and the warnings above the output are the reason.
nrow(markers_all)

table(markers_all$cluster)

saveRDS(markers_all, "RObjects/03g_markers_all.rds")
write_csv(markers_all, "results/tables/markers_all.csv")


#### Top markers per cluster ####

# pct.1 > 0.70 means the gene is expressed in at least 70 percent of the
# cluster's cells, so it is a marker OF the cluster rather than of a corner of
# it. abs(avg_log2FC) > 1 means at least a two-fold difference.
top_markers_all <- markers_all %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.70, abs(avg_log2FC) > 1) %>%
  slice_max(n = 15, order_by = avg_log2FC) %>%
  summarise(top_genes = paste(gene, collapse = ", "))

print(top_markers_all, n = Inf, width = Inf)

write_csv(top_markers_all, "results/tables/top_markers_per_cluster.csv")

# A less strict version, because some clusters lose everything at pct.1 > 0.7
top_markers_relaxed <- markers_all %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(n = 15, order_by = avg_log2FC) %>%
  summarise(top_genes = paste(gene, collapse = ", "))

print(top_markers_relaxed, n = Inf, width = Inf)

write_csv(top_markers_relaxed, "results/tables/top_markers_relaxed.csv")
