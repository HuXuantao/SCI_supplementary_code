# Build spatial transcriptomics Seurat object

    slice_cells <- colnames(seurat_obj)[seurat_obj[[spot_column]][, 1] == slice_id]
    needs_prefix <- any(startsWith(slice_cells, paste0(slice_id, "-")))

    if (needs_prefix && !is.null(image_obj@boundaries$centroids)) {
      image_obj@boundaries$centroids@cells <- paste0(
        slice_id,
        "-",
        image_obj@boundaries$centroids@cells
      )
    }

    image_obj <- image_obj[Cells(seurat_obj)]
    seurat_obj[[slice_id]] <- image_obj
  }

  seurat_obj
}

# -----------------------------------------------------------------------------
# Build spatial Seurat object
# -----------------------------------------------------------------------------

if (!dir.exists(counts_dir)) {
  stop("Count matrix directory does not exist: ", counts_dir)
}

counts <- Read10X(counts_dir, gene.column = 1)

spatial_obj <- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)

spatial_obj <- RenameAssays(
  spatial_obj,
  RNA = "Spatial"
)

metadata <- read_spatial_metadata(metadata_file)

if (nrow(metadata) != ncol(spatial_obj)) {
  warning(
    "The number of metadata rows does not match the number of spots. ",
    "Please confirm that metadata rows are ordered consistently with the count matrix."
  )
}

spatial_obj@meta.data <- data.frame(row.names = colnames(spatial_obj))
spatial_obj <- AddMetaData(spatial_obj, metadata = metadata)

spatial_obj <- attach_spatial_images(
  seurat_obj = spatial_obj,
  image_dir = image_dir,
  spot_column = spot_column,
  resolution = image_resolution
)

if (run_sctransform) {
  spatial_obj <- SCTransform(
    spatial_obj,
    assay = "Spatial",
    verbose = FALSE
  )
}

saveRDS(
  spatial_obj,
  file.path(output_dir, "spatial_seurat_object.rds")
)

sessionInfo()
