```r
  biomart = "ensembl",
  dataset = "rnorvegicus_gene_ensembl"
)

sce_list <- setNames(
  lapply(sample_ids, read_one_sample, input_dir = input_dir, rat_mart = rat_mart),
  sample_ids
)

saveRDS(sce_list, file.path(output_dir, "seurat_sample_list_qc.rds"))

# -----------------------------------------------------------------------------
# SCTransform and integration
# -----------------------------------------------------------------------------

future::plan(future::multisession)

sce_list <- future.apply::future_lapply(
  sce_list,
  function(x) {
    SCTransform(
      x,
      vst.flavor = "v2",
      vars.to.regress = "percent.mt",
      verbose = FALSE
    )
  }
)

future::plan(future::sequential)

saveRDS(sce_list, file.path(output_dir, "seurat_sample_list_sct.rds"))

features <- SelectIntegrationFeatures(
  object.list = sce_list,
  nfeatures = integration_features
)

sce_list <- PrepSCTIntegration(
  object.list = sce_list,
  anchor.features = features,
  verbose = FALSE
)

anchors <- FindIntegrationAnchors(
  object.list = sce_list,
  normalization.method = "SCT",
  anchor.features = features,
  verbose = FALSE
)

seurat_obj <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  verbose = FALSE
)

seurat_obj$group <- factor(
  seurat_obj$orig.ident,
  levels = sample_ids,
  ordered = TRUE
)

saveRDS(seurat_obj, file.path(output_dir, "seurat_integrated.rds"))

sessionInfo()
