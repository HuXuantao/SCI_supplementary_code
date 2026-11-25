# Build and integrate the five-condition SCI single-cell object from gene-level RNA counts.
# Inputs are already quantified. Mouse gene symbols in the input are retained.
#
# Set SCI_SCRNA_INPUT to either a directory with one 10x GEX folder per condition
# or the merged SCI_srt_clean.rds. For the merged RDS, set SCI_SCRNA_METADATA to
# its cell metadata CSV, with columns 'barcode' and 'group'. The merge suffix
# on a cell barcode is not used as a time-point label.
# SCI_SCRNA_OUTPUT selects the directory for the three Seurat checkpoints.

suppressPackageStartupMessages(library(Seurat))

sample_ids <- c("Ctrl", "D1", "D3", "D6", "D10")
input_path <- Sys.getenv("SCI_SCRNA_INPUT", file.path("data", "scRNA"))
metadata_file <- Sys.getenv("SCI_SCRNA_METADATA", file.path(dirname(input_path), "metadata.csv"))
output_dir <- Sys.getenv("SCI_SCRNA_OUTPUT", file.path("results", "scRNA"))

check_gene_ids <- function(ids) {
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids))) {
    stop("Gene identifiers are missing.")
  }
  ensembl_id <- grepl("^ENS[A-Z]*[0-9]+", ids)
  mouse_gene <- grepl("^ENSMUSG[0-9]+(\\.[0-9]+)?$", ids)
  transcript <- grepl("^(NM|NR|XM|XR)_[0-9]+", ids)
  if (any((ensembl_id & !mouse_gene) | transcript)) {
    stop("Supply mouse gene-level features. Other Ensembl species IDs and transcript IDs cannot be relabelled as mouse genes.")
  }
}

read_mouse_10x <- function(path) {
  # Read10X uses column 2 for the supplied symbols. Check column 1 as well:
  # symbols shared across species cannot identify the underlying reference.
  feature_file <- file.path(path, "genes.tsv")
  if (!file.exists(feature_file)) feature_file <- file.path(path, "features.tsv.gz")
  if (!file.exists(feature_file)) stop("No 10x features file in ", path)
  features <- read.delim(feature_file, header = FALSE, colClasses = "character",
                         quote = "", comment.char = "")
  if (ncol(features) >= 3L) features <- features[features[[3L]] == "Gene Expression", , drop = FALSE]
  if (ncol(features) < 2L || !nrow(features)) stop("No gene expression features in ", path)
  check_gene_ids(features[[1L]])
  counts <- Read10X(data.dir = path, gene.column = 2L)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  counts
}

make_sample <- function(counts, id) {
  check_gene_ids(rownames(counts))
  obj <- CreateSeuratObject(counts = counts, project = id, assay = "RNA")
  obj$orig.ident <- id
  Idents(obj) <- id
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-", assay = "RNA")
  obj
}

read_samples <- function(input_path, metadata_file, sample_ids) {
  if (dir.exists(input_path)) {
    samples <- lapply(sample_ids, function(id) {
      counts <- read_mouse_10x(file.path(input_path, id))
      RenameCells(make_sample(counts, id), add.cell.id = id)
    })
  } else {
    if (!file.exists(input_path)) stop("Missing RNA counts input: ", input_path)
    if (!file.exists(metadata_file)) stop("Missing barcode-to-group CSV: ", metadata_file)
    merged <- readRDS(input_path)
    if (!inherits(merged, "Seurat") || !"RNA" %in% Assays(merged) ||
        !"counts" %in% Layers(merged[["RNA"]], search = NA)) {
      stop("Expected a Seurat object with a joined RNA counts layer.")
    }
    counts <- LayerData(merged, assay = "RNA", layer = "counts")
    meta <- read.csv(metadata_file, colClasses = "character", check.names = FALSE)
    if (!all(c("barcode", "group") %in% names(meta)) || anyNA(meta$barcode) ||
        any(!nzchar(meta$barcode)) || anyDuplicated(meta$barcode)) {
      stop("The CSV must have group and barcode columns, with one unique barcode per row.")
    }
    group <- meta$group[match(colnames(counts), meta$barcode)]
    if (anyNA(group) || any(!nzchar(group))) stop("Each RNA cell must have a named group in the CSV.")
    samples <- lapply(sample_ids, function(id) {
      columns <- which(group == id)
      if (!length(columns)) stop("Missing cells for condition ", id)
      make_sample(counts[, columns, drop = FALSE], id)
    })
  }
  setNames(samples, sample_ids)
}

save_checkpoint <- function(object, filename) {
  target <- file.path(output_dir, filename)
  tmp <- tempfile(pattern = paste0(filename, "-"), tmpdir = output_dir)
  on.exit(unlink(tmp))
  saveRDS(object, tmp)
  if (!file.rename(tmp, target)) stop("Could not save checkpoint: ", target)
}

sce_list <- read_samples(input_path, metadata_file, sample_ids)
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
  stop("Could not create output directory: ", output_dir)
}
save_checkpoint(sce_list, "seurat_sample_list_qc.rds")

sce_list <- lapply(sce_list, SCTransform, vst.flavor = "v2",
                   vars.to.regress = "percent.mt", verbose = FALSE)
save_checkpoint(sce_list, "seurat_sample_list_sct.rds")

features <- SelectIntegrationFeatures(object.list = sce_list)
sce_list <- PrepSCTIntegration(object.list = sce_list, anchor.features = features, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = sce_list, normalization.method = "SCT",
                                  anchor.features = features, verbose = FALSE)
seurat_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
seurat_obj$group <- factor(seurat_obj$orig.ident, levels = sample_ids, ordered = TRUE)
save_checkpoint(seurat_obj, "seurat_integrated.rds")
sessionInfo()
