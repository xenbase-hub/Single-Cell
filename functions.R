# Take a single file location. Must be in .rds format
# Returns a Seurat object
load_seurat <- function(dir=NULL) {
	require(Seurat)
	obj <- readRDS(dir)
	print(paste0(get_seurat_basename(dir), " | ", ncol(obj), " cells, ", nrow(obj), " genes, ~ ", format(object.size(obj),"Gb")))
	return(obj)
}

# Accepts either a directory or a vector of file locations. Must be in .rds format
# Returns a list of named Seurat objects
load_multiple_seurat <- function(dir=NULL,file_list=NULL) {
	require(Seurat)
	if (is.null(file_list)) {
		seurat_files <- list.files(path=dir,pattern=".rds$",recursive=FALSE,full.names=TRUE)
	} else if (is.null(dir)) {
		seurat_files <- file_list
	}
	seurat.objects <- c()
	seurat.names <- c()
	num_files <- length(seurat_files)
	print(paste0("Loading ", num_files, " Seurat objects"))
	for (i in 1:num_files) {
		seurat.object <- load_seurat(seurat_files[i])
		seurat.objects <- c(seurat.objects,seurat.object)
		seurat.names <- c(seurat.names,get_seurat_basename(seurat_files[i]))
	}
	names(seurat.objects) <- seurat.names
	return(seurat.objects)
}

# Takes a single file location. Must be in .rds format
# Returns the file name
get_seurat_basename <- function(dir = NULL) {
	filename <- basename(dir)
	basename <- substr(filename,1,nchar(filename)-4)
	return(basename)
}

# Takes a directory of outs/filtered_feature_bc_matrix from Cellranger
# Returns a unprocessed seurat object
create_seurat <- function(dir=NULL) {
	require(Seurat)
    tryCatch(
        expr = {
			data <- Read10X(data.dir = dir)
        },
        error = function(e){
            data <- Read10X(data.dir = dir, gene.column=1)
        }
    )
	seurat.object <- CreateSeuratObject(counts = data, project = "seurat", min.cells = 3, min.features = 200)
	return(seurat.object)
}

# Takes a single Seurat object, and mode of "SCT" or "default"
# Returns the object with standard processing applied
process_seurat <- function(seurat=NULL,mode="SCT") {
	require(Seurat)
	require(sctransform)
	require(glmGamPoi)
	if (mode == "default") {
		seurat.object <- NormalizeData(seurat,verbose = FALSE)
		seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
		seurat.object <- ScaleData(seurat.object)
		seurat.object <- RunPCA(seurat.object)
		seurat.object <- FindNeighbors(seurat.object, dims = 1:20)
		seurat.object <- FindClusters(seurat.object, resolution = 0.8)
		seurat.object <- RunUMAP(seurat.object, dims = 1:20)
		return(seurat.object)
	} else if (mode == "SCT") {
		seurat.object <- SCTransform(seurat.object, vst.flavor = "v2")
		seurat.object <- RunPCA(seurat.object)
		seurat.object <- FindNeighbors(seurat.object, dims = 1:20)
		seurat.object <- FindClusters(seurat.object, resolution = 0.8)
		seurat.object <- RunUMAP(seurat.object, dims = 1:20)
		return(seurat.object)
	}
}

# Takes a list of Seurat objects
# Returns the CCA integrated Seurat object
integrate_seurat <- function(object_list=NULL, mode = "SCT") {
	require(Seurat)
	if (mode == "default") {
		Integration.anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:30, reduction="cca")
		Integration.combined <- IntegrateData(anchorset = Integration.anchors, dims = 1:30)
		DefaultAssay(Integration.combined) <- "integrated"
		Integration.combined <- ScaleData(Integration.combined)
		Integration.combined <- RunPCA(Integration.combined, verbose = TRUE)
		Integration.combined <- FindNeighbors(Integration.combined, reduction = "pca", dims = 1:20)
		Integration.combined <- FindClusters(Integration.combined, resolution = 0.5)
		Integration.combined <- RunUMAP(Integration.combined, reduction = "pca", dims = 1:20, return.model=TRUE)
		return(Integration.combined)
	} else if (mode == "SCT") {
		features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
		object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = features)
		anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:30, normalization.method = "SCT", anchor.features = features)
		Integration.combined <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "SCT")
		DefaultAssay(Integration.combined) <- "integrated"
		Integration.combined <- RunPCA(Integration.combined, verbose = TRUE)
		Integration.combined <- FindNeighbors(Integration.combined, reduction = "pca", dims = 1:20)
		Integration.combined <- FindClusters(Integration.combined, resolution = 0.5)
		Integration.combined <- RunUMAP(Integration.combined, reduction = "pca", dims = 1:20, return.model=TRUE)
		return(Integration.combined)
	}
}

# Takes a reference and query Seurat object and a column from the reference metadata
# Returns the query Seurat object with labels transferred 
transfer_labels <- function(reference=NULL,query=NULL,meta=NULL) {
	require(Seurat)
	transfer.anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30)
	predictions <- TransferData(anchorset = transfer.anchors, reference = reference, refdata = meta, dims = 1:30)
	seurat.object <- AddMetaData(query, metadata = predictions)
	return(seurat.object)
}

# Takes a reference and query Seurat object and a column from the reference metadata
# Returns the query Seurat object mapped to the reference 
map_query <- function(reference=NULL,query=NULL,meta=NULL) {
	require(Seurat)
	anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30, reference.reduction = "pca")
	mapped <- MapQuery(anchorset = anchors, reference = reference, query = query, refdata = meta, reference.reduction = "pca", reduction.model = "umap")
	return(mapped)
}

# Takes the output of map_query, a Seurat object with predicitons and reference UMAP applied
# Optionally, assign a name to the file
# Saves a UMAP plot of the query projected onto the reference
plot_projection <- function(seurat=NULL,name="map_query") {
	require(Seurat)
	require(ggplot2)
	DimPlot(seurat, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,label.size = 6, repel = T) + NoLegend()
	ggsave(file=paste0(name,".png"), width = 10, height = 10)
}

# Takes a Seurat object, a column in the metadata, option for labels or legend, and an optional name
# Saves a UMAP plot of the object
save_umap <- function(seurat=NULL,group=NULL,label=FALSE,name="umap") {
	require(Seurat)
	require(ggplot2)
	if (label) {
		DimPlot(seurat, group.by = group, label = T, repel = TRUE) + NoLegend()
	} else {
		DimPlot(seurat, group.by = group, label = F)
	}
	ggsave(file=paste0(name,".png"), width = 10, height = 10)
}

# Takes a Seurat object and a vector of gene names, and an optional name
# Saves a feature plot of the object
save_featureplot <- function(seurat=NULL,features=NULL,name="features") {
	require(Seurat)
	require(ggplot2)
	FeaturePlot(seurat,features=features,pt.size = 1, min.cutoff ="q10", order=TRUE)
	ggsave(file=paste0(name,".png"), width = 10, height = 10)
}

# Takes a Seurat object, a column in the metadata, and an optional name
# Wries to file the number of cells for each category
save_cell_counts <- function(seurat=NULL,ident=NULL,name="cell_counts") {
	require(Seurat)
	cell_counts <- as.data.frame(table(seurat@meta.data[,ident]))
	colnames(cell_counts) <- c(ident,"Frequency")
	print(cell_counts)
	save_table(cell_counts,name,rows=F)
}

# Takes a Seurat object and an optional name
# Writes a file containing the entirety of its metadata
save_metadata <- function(seurat=NULL, name="metadata") {
	require(Seurat)
	metadata <- seurat@meta.data
	save_table(metadata,name)
}

# Takes a Seurat object and an optional name
# Writes a file containing the raw counts data
save_raw_counts <- function(seurat=NULL,name="counts") {
	require(Seurat)
	counts <- seurat@assays$RNA@counts
	save_table(counts,name)
}

# Takes a Seurat object, a column in the metadata, and optional downsampling to X cells per ident
# Returns the markers for the given idents
get_markers <- function(seurat=NULL,idents=NULL,downsample=Inf,mode="SCT") {
	require(Seurat)
	require(MAST)
	Idents(seurat) <- idents
	if (mode == "default") {
		markers <- FindAllMarkers(seurat,logfc.threshold = 0.5,max.cells.per.ident=downsample,test.use = "MAST")
		return(markers)
	} else if (mode == "SCT") {
		seurat <- PrepSCTFindMarkers(seurat)
		markers <- FindAllMarkers(seurat,assay = "SCT",logfc.threshold = 0.5,max.cells.per.ident=downsample,test.use = "MAST")
		return(markers)
	}
}

# Takes the markers output by get_markers, and the amount of markers that should be kept for each cluster
# Returns the subset markers
subset_markers <- function(markers = NULL,amount=100) {
	require(tidyverse)
	markers_subset <- markers %>% arrange(desc(avg_log2FC)) %>% group_by(cluster) %>% slice(1:amount)
	return(markers_subset)
}

# Takes a Cellranger output directory
# Returns the corrected counts matrix 
run_soupx <- function(dir=NULL) {
	require(Seurat)
	require(SoupX)
	toc <- Seurat::Read10X(file.path(dir, "filtered_feature_bc_matrix"))
	tod <- Seurat::Read10X(file.path(dir, "raw_feature_bc_matrix"))
	clusters <- read.csv(file.path(dir,'/analysis/clustering/graphclust/clusters.csv'), header=T, row.names=1)
	projections <- read.csv(file.path(dir,'/analysis/tsne/2_components/projection.csv'), header=T, row.names=1)
	metadata <- cbind(clusters,projections)
	sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
	sc <- setClusters(sc, setNames(clusters$Cluster, rownames(clusters)))
	sc <- setDR(sc, projections, c("RD1", "RD2"))
	sc <- estimateSoup(sc)
	sc <- autoEstCont(sc, forceAccept = TRUE)
	out <- adjustCounts(sc)
	return(out)
}
