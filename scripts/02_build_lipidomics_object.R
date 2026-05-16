# Build lipidomics analysis object
  sample_names <- colnames(mat)
  groups <- gsub("_\d+$", "", sample_names)
  is_qc <- grepl(qc_regex, groups, ignore.case = TRUE)
  keep_idx <- !is_qc

  list(
    keep_idx = keep_idx,
    group_info = droplevels(factor(groups[keep_idx])),
    sample_names = sample_names[keep_idx]
  )
}

remove_qc_samples <- function(obj) {
  group_info <- make_group_info(obj$scaled)

  obj$scaled_no_qc <- obj$scaled[, group_info$keep_idx, drop = FALSE]
  obj$group_info <- group_info$group_info
  obj$sample_names <- group_info$sample_names

  obj
}

# -----------------------------------------------------------------------------
# Read and process lipidomics data
# -----------------------------------------------------------------------------

if (!file.exists(data_file_all)) {
  stop("Input file does not exist: ", data_file_all)
}

if (!file.exists(data_file_ms2)) {
  stop("Input file does not exist: ", data_file_ms2)
}

norm_intensity <- readxl::read_excel(data_file_all)
norm_ms2_intensity <- readxl::read_excel(data_file_ms2)

lipid_all <- process_lipidomics_file(
  data = norm_intensity,
  tag = "all_features",
  id_col = "MS1.name"
)

lipid_ms2 <- process_lipidomics_file(
  data = norm_ms2_intensity,
  tag = "ms2_features",
  id_col = "MS2.name"
)

lipid_all <- remove_qc_samples(lipid_all)
lipid_ms2 <- remove_qc_samples(lipid_ms2)

lipidomics_object <- list(
  all_features = lipid_all,
  ms2_features = lipid_ms2,
  preprocessing = list(
    sample_missing_threshold = sample_missing_threshold,
    metabolite_missing_threshold = metabolite_missing_threshold,
    transformation = "log2(x + 1)",
    normalization = "sample-wise sum normalization",
    scaling = "metabolite-wise autoscaling",
    qc_sample_rule = "sample names beginning with QC are excluded from scaled_no_qc"
  )
)

saveRDS(
  lipidomics_object,
  file.path(output_dir, "lipidomics_object.rds")
)

write.csv(
  lipid_all$metadata,
  file.path(output_dir, "metadata_all_features.csv"),
  row.names = FALSE
)

write.csv(
  lipid_ms2$metadata,
  file.path(output_dir, "metadata_ms2_features.csv"),
  row.names = FALSE
)

sessionInfo()
