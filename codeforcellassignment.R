#1.1. Run function to assign cell types
assign_celltypes_by_score <- function(seurat_obj, marker_lists, score_name_prefix = "CellTypeScore", min_score = 0) {
  # 1. Run AddModuleScore on the gene sets
  seurat_obj <- AddModuleScore(seurat_obj, features = marker_lists, name = score_name_prefix)
  
  # 2. Extract score columns
  score_cols <- paste0(score_name_prefix, 1:length(marker_lists))
  scores <- seurat_obj@meta.data[, score_cols]
  
  # 3. Cell type names (in the same order as marker_lists)
  cell_types <- names(marker_lists)
  
  # 4. Assign each cell to cell type with max score
  max_idx <- max.col(scores, ties.method = "first")
  max_scores <- apply(scores, 1, max)
  
  predicted <- ifelse(max_scores >= min_score, cell_types[max_idx], "Unknown")
  
  # 5. Add predicted cell type to metadata
  seurat_obj$predicted_celltype <- predicted
  
  return(seurat_obj)
}

#1.2. Define different markers for the respective cell types
marker_lists <- list(
  Pluripotent = c("POU5F1","NANOG","SOX2"),
  Mesendoderm = c("MIXL1","EOMES"),
  CardiacMesoderm = c("MESP1","MESP2","KDR","PDGFRA"),
  CardiacProg = c("ISL1","NKX2-5","TBX5"),
  Cardiomyocytes = c("TNNT2","TNNI3","ACTN2","MYH6","MYH7"),
  Epicardial = c("WT1","TBX18"),
  Fibroblast = c("COL1A1","COL3A1","DCN","LUM","VIM"),
  Endothelial = c("PECAM1","CDH5","VWF","KDR"),
  SmoothMusc = c("ACTA2","TAGLN","MYH11","CNN1"),
  Endoderm = c("SOX17","FOXA2","EPCAM"),
  NeuralCrest = c("SOX10","PAX3")
)


#1.3. Run the function on the seurat object
merged2<- assign_celltypes_by_score(merged2, marker_lists, min_score = 0.3)
