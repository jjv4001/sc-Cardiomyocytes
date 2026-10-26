#Load and create all Seurat objects
Pool1.data<- Read10X(data.dir = "Pool1.merged.filtered.matrix")
Pool1 <- CreateSeuratObject(counts = Pool1.data, min.cells = 3, min.features = 200)
Pool1metadata<-read.csv("Pool1.merged.allCells.csv")
Pool1$pool<-rep("Pool1")
Pool1$scaleplex<-Pool1metadata$assigned_scaleplex
Pool1sub <- subset(
  Pool1,
  !(scaleplex %in% c("Indeterminate", "Max_Fail","Unexpected"))
)
Pool1sub$scaleplex[Pool1sub$scaleplex == "XXX"] <- "XXX"

Pool2.data<- Read10X(data.dir = "Pool2.merged.filtered.matrix")
Pool2 <- CreateSeuratObject(counts = Pool2.data, min.cells = 3, min.features = 200)
Pool2metadata<-read.csv("Pool2.merged.allCells.csv")
Pool2$pool<-rep("Pool2")
Pool2$scaleplex<-Pool2metadata$assigned_scaleplex
Pool2sub <- subset(
  Pool2,
  !(scaleplex %in% c("Indeterminate", "Max_Fail","Unexpected"))
)
Pool2sub$scaleplex[Pool2sub$scaleplex == "XXX"] <- "XXX"


Pool2sub$chemical[Pool2sub$chemical == "XXX"] <- "XXX"
Pool3.data<- Read10X(data.dir = "Pool3.merged.filtered.matrix")
Pool3 <- CreateSeuratObject(counts = Pool3.data, min.cells = 3, min.features = 200)
Pool3metadata<-read.csv("Pool3.merged.allCells.csv")
Pool3$pool<-rep("Pool3")
Pool3$scaleplex<-Pool3metadata$assigned_scaleplex
Pool3sub <- subset(
  Pool3,
  !(scaleplex %in% c("Indeterminate", "Max_Fail","Unexpected"))
)
Pool3sub$scaleplex[Pool3sub$scaleplex == "XXX"] <- "XXX"


Pool4.data<- Read10X(data.dir = "Pool4.merged.filtered.matrix")
Pool4 <- CreateSeuratObject(counts = Pool4.data, min.cells = 3, min.features = 200)
Pool4metadata<-read.csv("Pool4.merged.allCells.csv")
Pool4$pool<-rep("Pool4")
Pool4$scaleplex<-Pool4metadata$assigned_scaleplex
Pool4sub <- subset(
  Pool4,
  !(scaleplex %in% c("Indeterminate", "Max_Fail","Unexpected"))
)
Pool4sub$scaleplex[Pool4sub$scaleplex == "XXX"] <- "XXX"


Pool5.data<- Read10X(data.dir = "Pool5.merged.filtered.matrix")
Pool5 <- CreateSeuratObject(counts = Pool5.data, min.cells = 3, min.features = 200)
Pool5metadata<-read.csv("Pool5.merged.allCells.csv")
Pool5$pool<-rep("Pool5")
Pool5$scaleplex<-Pool5metadata$assigned_scaleplex
Pool5sub <- subset(
  Pool5,
  !(scaleplex %in% c("Indeterminate", "Max_Fail","Unexpected"))
)
Pool5sub$scaleplex[Pool5sub$scaleplex == "XXX"] <- "XXX"

Pool6.data<- Read10X(data.dir = "Pool6.merged.filtered.matrix")
Pool6 <- CreateSeuratObject(counts = Pool6.data, min.cells = 3, min.features = 200)
Pool6metadata<-read.csv("Pool6.merged.allCells.csv")
Pool6$pool<-rep("Pool6")
Pool6$scaleplex<-Pool6metadata$assigned_scaleplex
Pool6sub <- subset(
  Pool6,
  !(scaleplex %in% c("Indeterminate", "Max_Fail","Unexpected"))
)
Pool6sub$scaleplex[Pool6sub$scaleplex == "XXX"] <- "XXX"

#merge all objects
merged2<-merge(Pool1sub, Pool2sub, Pool3sub, Pool4sub, Pool5sub, Pool6sub)
merged2 <- NormalizeData(merged2)
merged2 <- FindVariableFeatures(merged2, selection.method = "vst", nfeatures = 2000)
merged2 <- ScaleData(merged2)
merged2 <- RunPCA(merged2, features = VariableFeatures(object = merged2))
RunUMAP(merged2, dim=1:10)
merged2 <- FindNeighbors(merged2, dims = 1:10)
merged2 <- FindClusters(merged2, resolution = 0.5)

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
Idents(merged2)<-merged2$predicted.celltype

#saveallcelldata and cardiomyocyte only data
saveRDS(merged2, "cardiomyocyteassigned.rds")

cardiomyocytes<-subset(merged2, idents="Cardiomyocytes")

saveRDS(merged2, "cardiomyocytesonly.rds")

#Load cardiomyocyteonly data
cardiomyocytes<-readRDS("/athena/chenlab/scratch/jjv4001/Cardiomyocytescreening/cardiomyocytesonly.rds")
cardiomyocytes <- NormalizeData(cardiomyocytes)
cardiomyocytes <- FindVariableFeatures(cardiomyocytes, selection.method = "vst", nfeatures = 2000)
cardiomyocytes <- ScaleData(cardiomyocytes)
cardiomyocytes <- RunPCA(cardiomyocytes, features = VariableFeatures(object = cardiomyocytes))
RunUMAP(cardiomyocytes, dim=1:10)

#Load table with number of cells for each chemical treatment and plot
chemical<-table(cardiomyocytes$chemical)
chemical <- chemical %>% arrange(desc(Freq))
chemical$Index <- 1:359
ggplot(chemical, aes(x = Index, y = Freq)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = 100, color = "red", linetype = "dotted", size = 1) +
  labs(y = "Number of cells", x = "Index") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

#subset only ventricular cardiomyocytes data
VM<-subset(cardiomyocytes, idents="Ventricular cardiomyocytes")
VM <- NormalizeData(VM)
VM <- FindVariableFeatures(VM, selection.method = "vst", nfeatures = 2000)
VM <- ScaleData(VM)
VM <- RunPCA(VM, features = VariableFeatures(object = VM))
RunUMAP(VM, dim=1:10)
saveRDS(VM, "/athena/chenlab/scratch/jjv4001/Cardiomyocytescreening/cardiomyocytesonlyvm.rds")

#subset only atrial cardiomyocytes data
AM<-subset(cardiomyocytes, idents="Atrial cardiomyocytes")
AM <- NormalizeData(AM)
AM <- FindVariableFeatures(AM, selection.method = "vst", nfeatures = 2000)
AM <- ScaleData(AM)
AM <- RunPCA(AM, features = VariableFeatures(object = AM))
RunUMAP(AM, dim=1:10)
saveRDS(AM, "/athena/chenlab/scratch/jjv4001/Cardiomyocytescreening/cardiomyocytesonlyam.rds")

#Run differential gene expression analysis for all chemicals

dir.create("cardiomyocytes_deseq", showWarnings = FALSE)

# Loop through chemicals
for (chem in chemical$chemical) {
  
  # Run FindMarkers
  markers <- FindMarkers(cardiomyocytes, ident.1 = chem, ident.2 = "negative_control")
  
  
  # Add gene_symbol column
  markers$gene_symbol <- rownames(markers)
  
  # Move gene_symbol to the first column
  markers <- markers[, c("gene_symbol", setdiff(colnames(markers), "gene_symbol"))]
  
  # Clean chemical name for file name
  clean_name <- gsub("[^[:alnum:]_]", "_", chem)
  
  # Save to CSV
  write.csv(markers, file = paste0("cardiomyocytes_deseq/", clean_name, "_markers.csv"), row.names = FALSE)
}


dir.create("cardiomyocytes_deseq_VM", showWarnings = FALSE)

# Loop through chemicals
for (chem in chemical$chemical) {
  
  # Run FindMarkers
  markers <- FindMarkers(VM, ident.1 = chem, ident.2 = "negative_control")
  
  
  # Add gene_symbol column
  markers$gene_symbol <- rownames(markers)
  
  # Move gene_symbol to the first column
  markers <- markers[, c("gene_symbol", setdiff(colnames(markers), "gene_symbol"))]
  
  # Clean chemical name for file name
  clean_name <- gsub("[^[:alnum:]_]", "_", chem)
  
  # Save to CSV
  write.csv(markers, file = paste0("cardiomyocytes_deseq_VM/", clean_name, "_markers.csv"), row.names = FALSE)
}

dir.create("cardiomyocytes_deseq_AM", showWarnings = FALSE)

# Loop through chemicals
for (chem in chemical$chemical) {
  
  # Run FindMarkers
  markers <- FindMarkers(AM, ident.1 = chem, ident.2 = "negative_control")
  
  
  # Add gene_symbol column
  markers$gene_symbol <- rownames(markers)
  
  # Move gene_symbol to the first column
  markers <- markers[, c("gene_symbol", setdiff(colnames(markers), "gene_symbol"))]
  
  # Clean chemical name for file name
  clean_name <- gsub("[^[:alnum:]_]", "_", chem)
  
  # Save to CSV
  write.csv(markers, file = paste0("cardiomyocytes_deseq_AM/", clean_name, "_markers.csv"), row.names = FALSE)
}