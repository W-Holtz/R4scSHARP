#' Runs Scina
#'
#' @param data Gene expression matrix
#' @param markers Marker genes for cell classification
#' @param ref Reference dataset for predictions.
#'
#' @return Returns Scina predictions as an R data frame
#' @importFrom Seurat GetAssayData
#' @export
run_scina <- function(data, markers = NULL, ref = NULL) {
  # can return "unknown"

  if (is.null(markers)) {
    return(FALSE)
  }
  print("Running SCINA")
  norm_counts <- GetAssayData(object = data, slot = "data")
  results <- SCINA::SCINA(norm_counts, markers, rm_overlap = FALSE)
  scina_preds <- results$cell_labels

  scina_preds <- replace(scina_preds, scina_preds == "unknown", NA)
  return(scina_preds)

}

#' Runs scSorter
#'
#' @param data Gene expression matrix
#' @param markers Marker genes for cell classification
#' @param ref Reference dataset for predictions.
#'
#' @return Returns scSorter predictions as an R data frame
#' @importFrom Seurat GetAssayData VariableFeatures
#' @export
run_scsorter <- function(data, markers = NULL, ref = NULL) {
  #need top variable genes
  #can return "Unknown"

  if (is.null(markers)) {
    return(FALSE)
  }

  print("Running scSorter")
  types <- list()
  for (i in 1:length(names(markers))) {
    types <- append(types, rep(names(markers)[i], length(unlist(markers[i]))))
  }

  types <- unlist(types)
  markers <- unname(unlist(markers))
  anno <- data.frame(Type = types, Marker = markers)

  topgenes <- head(VariableFeatures(data), 2000)
  picked_genes <- unique(c(topgenes, anno$Marker))
  #expr <- as.matrix(data@assays$RNA@data)
  expr <- as.matrix(GetAssayData(object = data, slot = "data"))
  expr <- expr[rownames(expr) %in% picked_genes, ]

  rts <- scSorter::scSorter(expr, anno)
  scsort_preds <- rts$Pred_Type
  scsort_preds <- replace(scsort_preds, scsort_preds == "Unknown", NA)

  return(scsort_preds)
}

# GNU General Public License v3.0
# (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
#' gene_sets_prepare: prepare gene sets and calculate marker sensitivity
#' from input Cell Type excel file
#'
#' @param path_to_db_file - DB file with cell types
#' @param cell_type - cell type (e.g. Immune system, Liver, Pancreas, Kidney,
#' Eye, Brain)
#' @importFrom stats na.omit
#' @importFrom HGNChelper checkGeneSymbols
#
gene_sets_prepare <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(HGNChelper::checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(HGNChelper::checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
   
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

# GNU General Public License v3.0
# (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
#' sctype_score: calculate ScType scores and assign cell types
#'
#' @param scRNAseqData - input scRNA-seq matrix
#' (rownames - genes, column names - cells)
#' @param scaled - indicates whether the matrix is scaled (TRUE by default)
#' @param gs - list of gene sets positively expressed in the cell type
#' @param gs2 - list of gene sets that should not be expressed in the cell type
#' @param ... - additional parameter input
#' (NULL if not applicable)
#' @param gene_names_to_uppercase - if true, returns gene names in uppercase
#
sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                      gene_ = names(marker_stat), stringsAsFactors = !1)

  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
 
  es.max
}


#' Runs scType
#'
#' @param data Gene expression matrix
#' @param markers Marker genes for cell classification
#' @param ref Reference dataset for predictions.
#'
#' @return Returns scType predictions as an R data frame
#' @importFrom Seurat GetAssayData VariableFeatures
#' @export
run_sctype <- function(data, markers = NULL, ref = NULL) {
  if (is.null(markers)) {
    return(FALSE)
  }
  print("Running scType")
  norm_counts <- as.matrix(GetAssayData(object = data, slot = "data"))
  es.max = sctype_score(scRNAseqData = norm_counts, scaled = F,
                        gs = markers, gs2 = NULL, gene_names_to_uppercase = F)
  cL_resutls = do.call("rbind", lapply(unique(data@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$seurat_clusters==cl)), 10)
  }))

  cluster <- NULL
  scores <- NULL
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) <
    sctype_scores$ncells/4] = NA

  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster==j, ]
    data@meta.data$customclassif[data@meta.data$seurat_clusters == j] <-
      as.character(cl_type$type[1])
  }

  sctype_preds <- data$customclassif

  return(sctype_preds)
}

#' Runs Singler
#'
#' @param data Gene expression matrix
#' @param ref  Reference dataset for predictions
#' @param ref_labels Labels of reference dataset
#'
#' @return Returns SingleR predictions as an R data frame
#' @importFrom Seurat GetAssayData
#' @export
run_singler <- function(data, ref, ref_labels) {
  norm_counts <- as.matrix(GetAssayData(object = data, slot = "data"))
  results <- SingleR::SingleR(norm_counts, as.matrix(ref@assays$RNA@data),
    ref_labels[, 1])
  return(results$pruned.labels)
}

#' Runs scPred
#'
#' @param data Gene expression matrix
#' @param ref Reference dataset for predictions
#' @param ref_labels Labels of reference dataset
#'
#' @return Returns scPred predictions as an R data frame
#' @importFrom Seurat GetAssayData FindVariableFeatures ScaleData
#' RunPCA RunUMAP AddMetaData
#' @export
run_scpred <- function(data, ref, ref_labels) {
  # can return "unassigned"
  ref <- FindVariableFeatures(ref)
  ref <- ScaleData(ref)
  ref <- RunPCA(ref)
  ref <- RunUMAP(ref, dims = 1:30)
  ref <- AddMetaData(ref, ref_labels, col.name = "celltype")
  ref <- scPred::getFeatureSpace(ref, "celltype")
  ref <- scPred::trainModel(ref)
  data <- scPred::scPredict(data, ref)

  scpreds <- data$scpred_prediction
  scpreds <- replace(scpreds, scpreds == "unassigned", NA)

  return(scpreds)
}

#' Run Component Tools (Scina, scSorter, scType, scPred, Singler)
#' @description This function runs up to five cell classification
#' tools (SCINA, scSorter, scType, SingleR, and scPred) on an gene
#' expression matrix and returns the compiled predictions.
#' @param data_path Path to gene expression matrix stored as csv
#' @param tools Tool you would like to run
#' @param min_cells Filters data so cells with a sample size less then the
#' given amount are ignored (value is 0 by default)
#' @param min_feats Filters data so features with a sample size less then the
#' given amount are ignored (value is 0 by default)
#' @param markers Marker genes for cell classification
#' @param marker_names Names of the marker genes
#' @param ref_path Path to reference dataset for predictions
#' @param ref_labels_path Path to reference dataset labels
#'
#' @return Returns the tool predictions in a single R data frame object
#' where the rows are cells, and the columns are tool prediction output values.
#' @importFrom utils head read.csv
#' @import Seurat
#' @import dplyr
#' @export
run_tools <- function(data_path, tools, min_cells, min_feats, markers = NULL,
  marker_names = NULL, ref_path = NULL, ref_labels_path = NULL) {
  # add ability to take in seurat counts object

  counts <- read.csv(data_path, header = TRUE, row.names = 1)
  print(dim(counts))
  data <- CreateSeuratObject(t(counts), min.cells = min_cells,
    min.features = min_feats)
  print(dim(data@assays$RNA@counts))
  data <- NormalizeData(data)
  data <- ScaleData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data <- RunPCA(data)
  data <- FindNeighbors(data)
  data <- FindClusters(data)

  markers <- unname(markers)
  for (i in 1:length(markers)){
    markers[i] <- list(unlist(markers[i]))
  }

  names(markers) <- marker_names

  tools <- unlist(tools)
  results_df <- data.frame(start = rep(0, ncol(x = data)))
  if ("scina" %in% tools) {
    if (!requireNamespace("SCINA", quietly = TRUE)) {
      message("ERROR: The package SCINA has not been installed! 
        Please install that tool first, then try to run the program again.")
      message("For more help with installation, please see 
        https://github.com/W-Holtz/R4scSHARP#installation")
    } else {
    results_df$scina <- run_scina(data, markers)
    }
  }
  if ("scsorter" %in% tools) {
    if (!requireNamespace("scSorter", quietly = TRUE)) {
      message("ERROR: The package scSorter has not been installed! 
        Please install that tool first, then try to run the program again.")
      message("For more help with installation, please see 
        https://github.com/W-Holtz/R4scSHARP#installation")
    } else {
    results_df$scsorter <- run_scsorter(data, markers)
    }
  }
  if ("sctype" %in% tools) {
    results_df$sctype <- run_sctype(data, markers)
  }

  if ((!is.null(ref_path)) && (!is.na(ref_path))) {
    ref_counts <- read.csv(ref_path, header = TRUE, row.names = 1)
    ref_labels <- read.csv(ref_labels_path, header = TRUE, row.names = 1)
    ref <- CreateSeuratObject(t(ref_counts))
    ref <- NormalizeData(ref)

    if ("singler" %in% tools) {
      if (!requireNamespace("SCINA", quietly = TRUE)) {
      message("ERROR: The package SCINA has not been installed! 
        Please install that tool first, then try to run the program again.")
      message("For more help with installation, please see 
        https://github.com/W-Holtz/R4scSHARP#installation")
      } else {
      results_df$singler <- run_singler(data, ref, ref_labels)
      }
    }
    if ("scpred" %in% tools) {
      if (!requireNamespace("SCINA", quietly = TRUE)) {
      message("ERROR: The package SCINA has not been installed! 
        Please install that tool first, then try to run the program again.")
      message("For more help with installation, please see 
        https://github.com/W-Holtz/R4scSHARP#installation")
      } else {
      results_df$scpred <- run_scpred(data, ref, ref_labels)
      }
    }

  }

  results_df <- results_df[-1]
  row.names(results_df) <- row.names(data@meta.data)
  return(results_df)
}


#' R4scSHARP
#' @description This function is the primary function in the R4scSHARP
#' package. It calls another tool to run up to five cell classification
#' tools (SCINA, scSorter, scType, SingleR, and scPred) on an gene
#' expression matrix and saves the output file so scSHARP can run
#' predictions using that information.
#' @param data_path Path represented as a string to gene expression
#' matrix (or DGE matrix) stored as csv.
#' @param marker_path Path represented as a string to marker genes
#' for cell classification
#' @param ref_path Path represented as a string to reference dataset
#' for predictions
#' @param ref_label_path Path represented as a string to reference
#' dataset labels
#' @param out_path Path represented as a string to desired location for
#' the results to be saved. If no out_path parameter is given, no output
#' file will be generated.
#' @param tools Tools you would like to run (runs all tools by default).
#' Example Inputs:  "scina,scpred,singler", "sctype"
#' @param min_cells Filters data so cells with a sample size less then the
#' given amount are ignored (value is 0 by default)
#' @param min_feats Filters data so features with a sample size less then the
#' given amount are ignored (value is 0 by default)
#'
#' @return Returns the tool predictions in a single R data frame object
#' where the rows are cells, and the columns are tool prediction output values.
#' Save the returned data frame as a csv to a given output location.
#' @importFrom utils head read.csv write.csv
#' @import Seurat
#' @import dplyr
#' @export
run_r4scsharp <- function(data_path, marker_path, ref_path, ref_label_path,
  out_path = NULL, tools = "scina,scsorter,sctype,singler,scpred",
  min_cells = 0, min_feats = 0) {

  tools <- unlist(strsplit(tools, ","))

  print(tools)
  marker_df <- read.csv(marker_path, header = FALSE, row.names = 1)
  markers <- as.list(as.data.frame(t(marker_df)))
  markers <- lapply(markers, function(z) {
    z[!is.na(z) & z != ""]
    })
  marker_names <- row.names(marker_df)
  print(markers)
  print(marker_names)
  print(ref_path)
  results <- run_tools(data_path, tools, min_cells, min_feats,
    markers, marker_names, ref_path, ref_label_path)
  if (!is.null(out_path)) {
    write.csv(results, out_path)
  }
  return(results)
}