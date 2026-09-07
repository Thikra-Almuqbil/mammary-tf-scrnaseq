# ==============================================================
# SCRIPT 1 OF 4: IMPORT AND QUALITY CONTROL
#
# From 8 CellRanger folders on disk to one filtered Seurat object.
#
#   PART 1  Read the 8 matrices and build the Seurat object
#   PART 2  Calculate the QC metrics and look at them per sample
#   PART 3  Compare three threshold schemes BEFORE cutting anything
#   PART 4  Apply the filter I chose, and report retention per sample
#
# Dataset: Bach et al. (2017) Nat Commun, mouse mammary epithelium,
#          GEO GSE106273, 8 samples across 4 developmental stages.
#
# The order matters here and it is deliberate. Parts 2 and 3 remove nothing.
# I look at the real distributions and cost out three different rules first,
# and only then commit in part 4. Filtering is irreversible, so I did not want
# to do it from a default.
#
# Course reference:
#   Demonstrations/04_Preprocessing_And_QC.R, sections "Multiple samples",
#   "QC: remove undetected genes", "QC: visualising metrics",
#   "QC: filtering low-quality droplets", "Apply filtering across all samples"
#
# Expected runtime: about 15 minutes in total.
# Checkpoints written: RObjects/01a_imported.rds
#                      RObjects/01b_qc_metrics.rds
#                      RObjects/01c_filtered.rds
# ==============================================================


# I load the libraries I need for this script and nothing else, so it is
# obvious later which functions came from where.
library(Seurat)
library(tidyverse)

# One ggplot theme for every plot in the script
theme_set(theme_classic())

# I removed the setwd() line that used to be here, because it had my own
# laptop path in it and would not work on anyone else's machine.
# Set your working directory to the ROOT of this repository before running,
# not to the scripts folder. Every path below is relative to that root.
# In RStudio the easiest way is to open the repository folder as a project.

# Folders for the things this script produces. These warn if they already
# exist, which I ignore, that is not an error.
dir.create("RObjects", showWarnings = FALSE)
dir.create("results/graphs", recursive = TRUE, showWarnings = FALSE)



# ==============================================================
# PART 1: READING THE DATA IN
# ==============================================================

#### The sample information ####

# This is my own sample sheet. It lists the 8 samples, their GEO accession
# (which is what the folders on disk are named), and which developmental
# stage each mouse was in.
sampleinfo <- read_csv("data/raw/sample_info.csv")

# I always print it rather than trusting that it read correctly.
sampleinfo

# The sample names in my sheet use underscores: NP_1, G_1 and so on.
# This matters more than it looks. CreateSeuratObject() decides what to put in
# orig.ident by splitting each cell barcode on "_" and taking the first piece.
# So "NP_1" would become just "NP" and I would silently lose the replicate
# number, which is exactly the information I need later to argue that the
# stage separation is not a batch effect.
# I make a version of the name with a dash instead. This is also the style the
# course samples used, for example "ETV6RUNX1-1".
sampleinfo$SampleName <- str_replace(sampleinfo$sample, "_", "-")

sampleinfo


#### Pointing at the CellRanger matrices ####

# One file path per sample. The folders on disk are named by GEO accession,
# which is the sra_run column of my sample sheet.
list_of_files <- file.path("data/raw/cellranger",
                           sampleinfo$sra_run,
                           "outs/filtered_feature_bc_matrix")

# Naming this vector is what makes Read10X put a sample prefix on every cell
# barcode. Without the names I could not tell the samples apart afterwards.
names(list_of_files) <- sampleinfo$SampleName

# I check the paths exist before trying to read them, because Read10X's error
# message when a path is wrong is not very helpful.
list_of_files
file.exists(list_of_files)


#### Making the gene file name match what Read10X expects ####

# This is the first thing that went wrong for me and it took a while to work
# out. These matrices came from CellRanger version 2, which names the gene
# annotation file "genes.tsv.gz". Version 3 onwards calls it
# "features.tsv.gz", and that is the name Read10X() looks for.
#
# Read10X() does have a fallback for older data, but it only triggers when it
# finds an UNCOMPRESSED "genes.tsv". Mine is gzipped, so the fallback never
# fires and Read10X() stops with:
#   "Gene name or features file missing. Expecting features.tsv.gz"
#
# I fix it by putting a COPY of the gene file under the name Read10X() wants.
# I copy rather than rename so the original files are untouched, in case I
# need to go back to them. file.copy() is vectorised, so one call does all 8.

genes_files <- file.path(list_of_files, "genes.tsv.gz")
features_files <- file.path(list_of_files, "features.tsv.gz")

genes_files

# Make the copies.
# If I run this a second time it returns FALSE for every sample, because the
# copy already exists and file.copy() will not overwrite by default.
# That is fine and does not mean anything failed.
file.copy(from = genes_files, to = features_files)

# I confirm all 8 copies now exist rather than assuming the copy worked
file.exists(features_files)


#### Reading the matrices ####

# One call reads all 8 samples and stitches them into one matrix, because
# list_of_files is a named vector.
expression_matrix <- Read10X(data.dir = list_of_files)

# Genes by cells. I expect roughly 27,000 genes and 25,806 cells.
dim(expression_matrix)


#### Building the Seurat object ####

seurat_object <- CreateSeuratObject(counts = expression_matrix)
seurat_object

# The metadata is reached with double square brackets. At this point it has
# orig.ident, nCount_RNA and nFeature_RNA and nothing else.
head(seurat_object[[]])


#### Adding my sample information to the metadata ####

# I pull the metadata out as a normal data frame, add the two columns I want,
# and put it back. I do it this way rather than in one clever line because it
# is easier to check at each step.
temp_metadata <- seurat_object[[]] %>%
  # I turn the cell barcodes into a column so dplyr does not lose them,
  # because dplyr verbs drop rownames.
  rownames_to_column("Cell") %>%
  # The developmental stage is the part of orig.ident before the dash.
  # The pattern "-.*" matches the first dash and everything after it, so
  # "NP-1" becomes "NP".
  mutate(SampleGroup = str_remove(orig.ident, "-.*")) %>%
  # And I keep the full sample name too, so I can tell the two mice apart.
  mutate(SampleName = orig.ident) %>%
  # Put the barcodes back where Seurat expects them.
  column_to_rownames("Cell")

seurat_object[[]] <- temp_metadata

head(seurat_object[[]])

# I make both new columns factors with the levels in developmental order.
# If I leave them as plain text, R sorts alphabetically and every plot from
# here on reads G, L, NP, PI, which is biologically meaningless. Doing it once
# here means I never have to think about it again.
seurat_object$SampleGroup <- factor(seurat_object$SampleGroup,
                                    levels = c("NP", "G", "L", "PI"))

seurat_object$SampleName <- factor(seurat_object$SampleName,
                                   levels = c("NP-1", "NP-2",
                                              "G-1", "G-2",
                                              "L-1", "L-2",
                                              "PI-1", "PI-2"))


#### Checking the import actually worked ####

# Cells per sample, now in developmental order
table(seurat_object$SampleName)

# Cells per developmental stage.
# I got NP 4376, G 6021, L 9603, PI 5806, which is 25,806 in total and matches
# the number of cells Bach et al. report before their own QC. That agreement
# is the main evidence that I imported the right thing.
table(seurat_object$SampleGroup)

# Are the row names gene symbols rather than Ensembl IDs? They need to be
# symbols for everything downstream, because the marker lists and the
# transcription factor list are both in symbols.
head(rownames(seurat_object))

# Mouse mitochondrial gene symbols start with a LOWERCASE "mt-", not "MT-"
# like the human data used in the course. I check they are actually present
# here, because in part 2 a wrong pattern would give me zero percent
# mitochondrial content for every cell, with no error at all.
str_subset(rownames(seurat_object), "^mt-")


#### CHECKPOINT 1 ####

# I save at the end of every part so that a crash never costs me much.
# This one is about 140 MB. Safe stopping point.
saveRDS(seurat_object, "RObjects/01a_imported.rds")



# ==============================================================
# PART 2: QC METRICS AND PER-SAMPLE DIAGNOSTICS
#
# NOTHING IS FILTERED IN THIS PART.
#
# Why per sample and not one global cutoff: my 8 samples are mice in four very
# different physiological states. A lactating mammary gland is a milk protein
# factory, so those cells genuinely have different library sizes and different
# numbers of detected genes from a nulliparous mouse. If I applied one cutoff
# across everything, I would preferentially delete lactation cells for
# biological reasons and then compare stages that I had myself biased.
# ==============================================================

# If you are resuming here rather than running the whole script:
# seurat_object <- readRDS("RObjects/01a_imported.rds")


#### Removing genes that are never detected ####

# A gene with zero counts in all 25,806 cells carries no information and only
# slows everything down. I take the raw count matrix out first so the next two
# lines are readable.
raw_counts <- seurat_object[["RNA"]]$counts

# How many genes are detected at all?
table(rowSums(raw_counts) > 0)

# Keep only the genes with at least one count somewhere.
# Note this filters GENES, not cells. No cell is removed here.
filtered_seurat_object <- subset(
  seurat_object,
  features = rownames(seurat_object)[rowSums(raw_counts) > 0]
)

filtered_seurat_object

# I verify rather than assume: every remaining gene should have a nonzero sum.
filt_raw_counts <- filtered_seurat_object[["RNA"]]$counts
table(rowSums(filt_raw_counts) > 0)


#### Mitochondrial content ####

# A cell with a very high fraction of mitochondrial reads has usually lost its
# cytoplasm through a hole in the membrane, so the nuclear mRNA leaked out and
# the mitochondria stayed. It is the standard marker of a dying cell.
#
# THE PATTERN IS LOWERCASE. Mouse mitochondrial symbols are "mt-Nd1" and so
# on. The course used "^MT-" because its data were human. If I had used the
# uppercase pattern here every cell would come out at exactly 0 percent, with
# no error and no warning, and my QC would have done nothing at all.
filtered_seurat_object[["percent.mt"]] <-
  PercentageFeatureSet(filtered_seurat_object, pattern = "^mt-")

head(filtered_seurat_object[[]])

# The sanity check that catches the uppercase mistake: this must NOT be all
# zeros. Mine came out with a median of about 0.75 percent and a maximum of
# about 8.2 percent, which is a low-mitochondrial dataset.
summary(filtered_seurat_object$percent.mt)


#### Looking at the metrics, one violin per sample ####

# I look at all three metrics per sample rather than pooled, because the whole
# point is to see whether the samples differ from each other.
# layer = "counts" because I want the raw numbers, not normalised values.
# pt.size = 0 hides the individual points, unreadable at 25,000 cells.
VlnPlot(filtered_seurat_object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "SampleName",
        ncol = 3,
        layer = "counts",
        pt.size = 0)

# The same three metrics grouped by stage instead, which is easier to read
VlnPlot(filtered_seurat_object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "SampleGroup",
        ncol = 3,
        layer = "counts",
        pt.size = 0)

# The metadata as a plain data frame, for the ggplot and the tables below
metadata <- filtered_seurat_object[[]]

# Genes detected against total UMIs, coloured by mitochondrial content, one
# panel per sample. Empty droplets and dying cells sit in the bottom left with
# high mitochondrial content, so this is the plot where they become visible.
ggplot(metadata,
       aes(x = nFeature_RNA,
           y = nCount_RNA,
           colour = percent.mt)) +
  geom_point(size = 0.3) +
  scale_colour_viridis_c() +
  facet_wrap(~ SampleName, nrow = 2) +
  labs(title = "QC metrics per sample",
       x = "Number of genes detected",
       y = "Number of UMIs",
       colour = "Percent\nmitochondrial")


#### The thresholds the course would apply, calculated but NOT applied ####

# The course rule, within each sample:
#   genes:  median - 2 * MAD    (drop cells below)
#   UMIs:   median - 2 * MAD    (drop cells below)
#   mito:   median + 2 * MAD    (drop cells above)
#
# MAD is the median absolute deviation, a spread measure that is not dragged
# around by outliers the way the standard deviation is.
#
# The important part is group_by(SampleName). It makes median() and mad()
# calculate inside each sample, so every cell is judged against its own
# sample's distribution, not against a pooled distribution dominated by
# whichever stage contributed the most cells.

qc_thresholds <- metadata %>%
  group_by(SampleName) %>%
  summarise(n_cells = n(),
            median_genes = median(nFeature_RNA),
            min_genes = median(nFeature_RNA) - 2 * mad(nFeature_RNA),
            median_umis = median(nCount_RNA),
            min_umis = median(nCount_RNA) - 2 * mad(nCount_RNA),
            median_mt = median(percent.mt),
            max_mt = median(percent.mt) + 2 * mad(percent.mt))

# width = Inf stops the tibble printing from cutting columns off
print(qc_thresholds, width = Inf)


#### CHECKPOINT 2 ####

# QC metrics attached, NO cells removed. I save this separately from the
# filtered version so that if I change my mind about the thresholds I do not
# have to recompute anything. Safe stopping point.
saveRDS(filtered_seurat_object, "RObjects/01b_qc_metrics.rds")



# ==============================================================
# PART 3: COMPARING THREE THRESHOLD SCHEMES
#
# STILL NOTHING IS FILTERED.
#
# I wrote this part because I did not want to accept the course's default
# thresholds without knowing what they would cost me on this dataset. A
# MAD-based threshold is RELATIVE: it always cuts a slice off the bottom of
# whatever distribution you hand it, whether or not there is anything wrong
# with the cells. So "2 MADs" is not automatically the right answer, it just
# always removes something.
#
# For reference, the paper's own QC kept 25,010 of 25,806 cells:
#     NP 4223/4376 = 96.5 percent
#     G  5826/6021 = 96.8 percent
#     L  9319/9603 = 97.0 percent
#     PI 5642/5806 = 97.2 percent
# Their total also removes non-unique barcodes to guard against index
# swapping, which I am not doing, so their QC-only retention is a little
# higher than those figures suggest.
# ==============================================================

# If you are resuming here:
# filtered_seurat_object <- readRDS("RObjects/01b_qc_metrics.rds")
# metadata <- filtered_seurat_object[[]]


#### Scheme A: the course version, 2 MADs on all three metrics ####

scheme_a <- metadata %>%
  group_by(SampleName) %>%
  summarise(n_cells = n(),
            min_genes = median(nFeature_RNA) - 2 * mad(nFeature_RNA),
            min_umis  = median(nCount_RNA)   - 2 * mad(nCount_RNA),
            max_mt    = median(percent.mt)   + 2 * mad(percent.mt),
            n_keep = sum(nFeature_RNA > (median(nFeature_RNA) - 2 * mad(nFeature_RNA)) &
                         nCount_RNA   > (median(nCount_RNA)   - 2 * mad(nCount_RNA)) &
                         percent.mt   < (median(percent.mt)   + 2 * mad(percent.mt))))

scheme_a$pct_kept <- round(100 * scheme_a$n_keep / scheme_a$n_cells, 1)

# The column to look at here is max_mt. Mine came out at roughly 1.2 to 1.5
# percent per sample, when the whole dataset only reaches 8.2 percent at its
# very worst. That threshold is not identifying dying cells, it is just
# shaving the top of a narrow distribution because it was told to.
print(scheme_a, width = Inf)


#### Scheme B: 3 MADs on genes and UMIs, fixed 5 percent mitochondrial ####

scheme_b <- metadata %>%
  group_by(SampleName) %>%
  summarise(n_cells = n(),
            min_genes = median(nFeature_RNA) - 3 * mad(nFeature_RNA),
            min_umis  = median(nCount_RNA)   - 3 * mad(nCount_RNA),
            max_mt    = 5,
            n_keep = sum(nFeature_RNA > (median(nFeature_RNA) - 3 * mad(nFeature_RNA)) &
                         nCount_RNA   > (median(nCount_RNA)   - 3 * mad(nCount_RNA)) &
                         percent.mt   < 5))

scheme_b$pct_kept <- round(100 * scheme_b$n_keep / scheme_b$n_cells, 1)

print(scheme_b, width = Inf)


#### Scheme C: the paper's rule, 3 MADs with absolute floors ####

# pmax() compares two values and returns the larger.
# pmax(median - 3 * MAD, 500) therefore means "use three MADs below the
# median, unless that lands below 500 genes, in which case use 500".
# So the floor is a safety net: no matter how tight a sample's distribution
# is, I never keep a cell with fewer than 500 genes or 1000 UMIs.
# This is exactly the rule described in the paper's Methods.

scheme_c <- metadata %>%
  group_by(SampleName) %>%
  summarise(n_cells = n(),
            min_genes = pmax(median(nFeature_RNA) - 3 * mad(nFeature_RNA), 500),
            min_umis  = pmax(median(nCount_RNA)   - 3 * mad(nCount_RNA), 1000),
            max_mt    = 5,
            n_keep = sum(nFeature_RNA > pmax(median(nFeature_RNA) - 3 * mad(nFeature_RNA), 500) &
                         nCount_RNA   > pmax(median(nCount_RNA)   - 3 * mad(nCount_RNA), 1000) &
                         percent.mt   < 5))

scheme_c$pct_kept <- round(100 * scheme_c$n_keep / scheme_c$n_cells, 1)

print(scheme_c, width = Inf)


#### The three side by side ####

comparison <- data.frame(
  SampleName    = scheme_a$SampleName,
  n_cells       = scheme_a$n_cells,
  A_course_2mad = scheme_a$pct_kept,
  B_3mad        = scheme_b$pct_kept,
  C_paper       = scheme_c$pct_kept
)

comparison

# Totals out of 25,806.
# I got A = 23,696 (so 2,110 removed) and C = 25,802 (so 4 removed).
# That is the number that decided it. Scheme A was going to delete 8.2 percent
# of my cells, with retention varying from 88.8 to 93.5 percent BETWEEN
# samples, which is precisely the stage-dependent bias I was trying to avoid.
# And it was doing that almost entirely on a mitochondrial threshold of about
# 1.3 percent, which is not a biologically meaningful cutoff for dying cells.
sum(scheme_a$n_keep)
sum(scheme_b$n_keep)
sum(scheme_c$n_keep)


#### Where do the thresholds actually sit on the distributions? ####

# I did not want to pick a scheme from a table alone, so I drew the paper's
# thresholds onto the real per-sample violins to see where they land.

VlnPlot(filtered_seurat_object,
        features = "nFeature_RNA",
        group.by = "SampleName",
        layer = "counts",
        pt.size = 0) +
  geom_hline(yintercept = 500, colour = "red") +
  ggtitle("Genes detected per cell, with the paper's floor of 500")

VlnPlot(filtered_seurat_object,
        features = "percent.mt",
        group.by = "SampleName",
        layer = "counts",
        pt.size = 0) +
  geom_hline(yintercept = 5, colour = "red") +
  ggtitle("Mitochondrial percentage per cell, with the paper's cutoff of 5 percent")

# Seeing these two plots is what convinced me. The 500-gene line sits below
# the bottom of every sample's violin, and the 5 percent line sits above
# almost every cell. In other words this dataset was already clean, and the
# honest thing to report is that QC removed almost nothing, not to invent a
# threshold that removes a respectable-looking number of cells.



# ==============================================================
# PART 4: APPLYING THE FILTER I CHOSE
#
# Decision, after part 3: I use the thresholds published in Bach et al.,
# computed separately within each sample:
#     genes: max(median - 3 MAD, 500)
#     UMIs:  max(median - 3 MAD, 1000)
#     mito:  fixed 5 percent
#
# I am not claiming the paper's numbers are better in general. I am saying
# that on THIS data an absolute floor plus a fixed mitochondrial cutoff
# reflects what the distributions actually look like, and a purely relative
# cutoff does not.
# ==============================================================

#### For the record: what drove the course scheme's removals? ####

# This is NOT part of my filter. I calculate it so that I can say in the
# presentation exactly which criterion was responsible for the 2,110 cells
# the 2-MAD scheme would have removed, instead of just asserting it was the
# mitochondrial one.
# The three counts will not sum to the total removed, because one cell can
# fail more than one criterion at once.

scheme_a_breakdown <- metadata %>%
  group_by(SampleName) %>%
  summarise(n_cells = n(),
            fail_genes = sum(nFeature_RNA <= (median(nFeature_RNA) - 2 * mad(nFeature_RNA))),
            fail_umis  = sum(nCount_RNA   <= (median(nCount_RNA)   - 2 * mad(nCount_RNA))),
            fail_mt    = sum(percent.mt   >= (median(percent.mt)   + 2 * mad(percent.mt))))

print(scheme_a_breakdown, width = Inf)

sum(scheme_a_breakdown$fail_genes)
sum(scheme_a_breakdown$fail_umis)
sum(scheme_a_breakdown$fail_mt)


#### The thresholds I am actually using, written down ####

# I print this table and keep it, because "we used the paper's thresholds" is
# not a defensible answer on its own. If someone asks me in the presentation
# what the cutoff was for the L-1 sample, I need the number.
qc_thresholds_used <- metadata %>%
  group_by(SampleName) %>%
  summarise(n_cells = n(),
            median_genes = median(nFeature_RNA),
            min_genes = pmax(median(nFeature_RNA) - 3 * mad(nFeature_RNA), 500),
            median_umis = median(nCount_RNA),
            min_umis = pmax(median(nCount_RNA) - 3 * mad(nCount_RNA), 1000),
            median_mt = median(percent.mt),
            max_mt = 5)

print(qc_thresholds_used, width = Inf)


#### Applying it ####

# This follows the course's "Apply filtering across all samples" block.
# group_by(SampleName) is doing the real work: the median() and mad() inside
# filter() are computed within each sample, so a lactation cell is compared
# with other lactation cells and not with the dataset as a whole.
# pmax(..., 500) adds the paper's absolute floor on top.

all_keep_cells <- metadata %>%
  # dplyr drops rownames, so I move the barcodes into a column first
  rownames_to_column("cell") %>%
  # thresholds are computed inside each sample
  group_by(SampleName) %>%
  # the three criteria
  filter(nFeature_RNA > pmax(median(nFeature_RNA) - 3 * mad(nFeature_RNA), 500),
         nCount_RNA   > pmax(median(nCount_RNA)   - 3 * mad(nCount_RNA), 1000),
         percent.mt   < 5) %>%
  # and I pull out just the barcodes of the survivors
  pull(cell)

# How many cells passed? I got 25,802.
length(all_keep_cells)

qc_seurat_object <- subset(
  filtered_seurat_object,
  cells = all_keep_cells
)

qc_seurat_object


#### Reporting retention PER SAMPLE ####

# The whole justification for per-sample QC falls apart if I do not check
# that no single sample lost far more than the others. This is that check.
qc_summary <- data.frame(
  SampleName = names(table(metadata$SampleName)),
  before = as.numeric(table(metadata$SampleName)),
  after  = as.numeric(table(qc_seurat_object$SampleName))
)

qc_summary$removed  <- qc_summary$before - qc_summary$after
qc_summary$pct_kept <- round(100 * qc_summary$after / qc_summary$before, 1)

qc_summary

# Totals: 25,806 before, 25,802 after, so 4 cells removed.
sum(qc_summary$before)
sum(qc_summary$after)
sum(qc_summary$removed)

# Cells per stage after filtering, to line up against the paper's post-QC
# counts (they report NP 4223, G 5826, L 9319, PI 5642). Mine are higher,
# because their pipeline also removes non-unique barcodes as an index-swapping
# precaution and mine does not.
table(qc_seurat_object$SampleGroup)


#### The metrics after filtering ####

# These violins should look almost identical to the ones in part 2, because I
# only removed 4 cells. If they looked different I would have made a mistake.
VlnPlot(qc_seurat_object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "SampleName",
        ncol = 3,
        layer = "counts",
        pt.size = 0)


#### An aside: how big is the index-swapping problem here? ####

# The paper removed barcodes that appear in more than one sample, because on
# a patterned flow cell reads can be misassigned between multiplexed samples
# ("index swapping"), and a barcode showing up twice is the symptom.
# I did not implement their removal, so I measured how many cells it would
# have affected instead of guessing. This is the honest version of "our count
# is slightly higher than theirs".

# Strip the sample prefix I added in part 1, leaving the bare 10x barcode
bare_barcodes <- str_remove(colnames(qc_seurat_object), ".*_")

# How many barcode sequences appear more than once across the 8 samples?
sum(duplicated(bare_barcodes))

# And how many cells in total are involved, counting every copy?
# This is the number the paper's rule would have removed.
sum(bare_barcodes %in% bare_barcodes[duplicated(bare_barcodes)])


#### CHECKPOINT 3: END OF SCRIPT 1 ####

# This is the object script 2 reads.
saveRDS(qc_seurat_object, "RObjects/01c_filtered.rds")
