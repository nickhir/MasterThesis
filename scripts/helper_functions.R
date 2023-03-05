# simple function, which runs Seurats FindConservedMarkers for all clusters
FindAllConservedMarkers <- function(object, cores, clusters=NULL,...){
  library(pbmcapply)
  
  # if no clusters specified, all are used. One after another.
  if (is.null(clusters)){
    clusters <- stringr::str_sort(levels(Idents(object)), numeric = T)
  }
  message(stringr::str_interp("Looking for conserved markers in ${length(clusters)} clusters"))
  df <- pbmclapply(clusters, function(cluster_id){
    tmp_marker <- FindConservedMarkers(object, 
                                       ident.1 = cluster_id, 
                                       verbose = F,
                                       ...)
    tmp_marker <- tmp_marker %>%
      dplyr::mutate(cluster=cluster_id) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::relocate("gene")
  
    return(tmp_marker)
  }, mc.cores = cores) %>% 
    bind_rows()
  return(df)
}

# simple function which plots a violin plot and feature plot site by site.
# take one gene as input, can be used in lapply and then the plots can be stiched together with pathwork
plotFeatureExpression <- function(object, gene, pt.size_=0, all_plots=F, reduction_="umap", order_=T){
  library(patchwork)
  vln_plot <- VlnPlot(object, features = gene, ncol=1, pt.size=pt.size_) + 
    xlab("Cluster") + NoLegend()
  vln_plot <- vln_plot %>% ggrastr::rasterise(dpi=150)
  
  feature_plot <- FeaturePlot(object, features = gene, ncol=1,reduction = reduction_, order = order_) %>% 
    ggrastr::rasterise(dpi=150)
  
  if (all_plots){
    cluster_plot <- DimPlot(object,reduction = reduction_,label = T,label.box = T,repel = T)+ NoLegend() 
    cluster_plot <- cluster_plot %>% ggrastr::rasterise(dpi=150)
    return(wrap_plots(cluster_plot,
                      feature_plot,
                      vln_plot,ncol = 3) +
             plot_layout(widths = c(1,1, 1.6)))
  } else{
    return(wrap_plots(feature_plot,vln_plot,ncol = 2)+
             plot_layout(widths = c(1, 1.6)))
  }
}

plot_integrated_clusters = function (srat, normalize=T, show_ideal_proportion=F,
                                     levels_=NULL) {
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  
  
  
  count_table <- table(Idents(srat), srat@meta.data$orig.ident)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)

  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  sorted_labels <- cluster_size %>% arrange(value) %>% pull(cluster)
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)

  colnames(melt_mtx)[2] <- "Experiment"
  
  # if the number of cells differs strongly between experiments, it makes sense to "normalize" them
  if (normalize==T){
    #get the library sizes of experimnt
    lib_size <- table(srat@meta.data$orig.ident) %>% enframe(name="Experiment",value="num_cells")
    melt_mtx <- melt_mtx %>%
      data.frame() %>%
      left_join(.,lib_size, by="Experiment") %>%
      mutate(value = value/num_cells)
  }
  
  if (!is.null(levels_)){
    melt_mtx$Experiment <- factor(melt_mtx$Experiment, levels=levels_)
  }
  
  p1 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=Experiment)) + 
    geom_bar(position="fill", stat="identity") + 
    theme_bw() %+%
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
          legend.text = element_text(size=16.5))+ 
    scale_y_continuous(expand = c(0.01, 0.01))+
    coord_flip() + 
    scale_fill_brewer(palette = ifelse(n_distinct(srat@meta.data$orig.ident)>8,"Set3","Set2")) +
    ylab("Fraction of cells in each experiment") + xlab("Cluster") + theme(legend.position="top")
  
  # add the expected proportion, given that all clusters are perfectly mixed
  if (show_ideal_proportion==T){
    interval <- seq(0,1,1/n_distinct(srat@meta.data$orig.ident))
    interval <- interval[2:(length(interval)-1)]
    
    p1 <- p1 + geom_hline(yintercept = interval, alpha=0.3)
    
  }
  p2 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_point(color= "grey60",shape=18,size=4) + 
    theme(plot.margin = unit(c(0,0,0,0), "cm")) %+% 
    theme_bw() + scale_x_log10() + xlab("Cells per cluster") + ylab("")
  
  
  p1 + p2 + plot_layout(widths = c(2.8,1.2))
}


# simply plot the expression of two genes against one another. 
# remove all cells that show 0 expression.
plotCorr <- function(object,gene1,gene2){
  data <- object@assays$RNA@data
  data_oi <- rbind(data[gene1,], data[gene2,])
  # only keep values that are not 0
  data_oi <- data_oi[,colSums(data_oi)!=0]
  plot_data <- data_oi %>% 
    t() %>% 
    data.frame()
  colnames(plot_data) <- c(gene1, gene2)
  p<- ggplot(plot_data,aes(x=get(gene1), y=get(gene2)))+
    geom_point(size=0.2)
  return(p)
}


# for 2 given genes. returns a tibble which specifies the percentage of cells that expressed gene 1, gene 2 or both
qunatify_expression <- function(object,gene1,gene2){
  data <- object@assays$RNA@data
  data_oi <- rbind(data[gene1,], data[gene2,])
  
  both_not_expressed=sum(colSums(data_oi)==0)/ncol(data)
  gene1_expressed = sum(data_oi[1,] != 0 & data_oi[2,]==0)/ncol(data)
  gene2_expressed = sum(data_oi[1,] == 0 & data_oi[2,]!=0)/ncol(data)
  both_expressed = sum(data_oi[1,] != 0 & data_oi[2,]!=0)/ncol(data)
  res <- tibble(gene1=gene1,
                gene2=gene2,
                both_not_expressed=both_not_expressed,
                gene1_expressed=gene1_expressed,
                gene2_expressed=gene2_expressed,
                both_expressed=both_expressed)
  return(res)
}



# easy function to plot heatmap for set of genes. each column is one cell
marker_heatmap <- function(seurat, markers, celltype, group.by, cap_value=NULL){
  library(ComplexHeatmap)
  library(circlize)
  # plot a heatmap, where each column is one cell! The cells are grouped and split according to the cluster they are from 
  
  ordered_index <- order(seurat@meta.data[[group.by]])
  c_split <- sort(seurat@meta.data[[group.by]])
  
  if (!all(markers %in% rownames(seurat@assays$RNA))){
    mising_genes <- paste(markers[!markers %in% rownames(seurat@assays$RNA)], sep=", ")
    celltype <- celltype[markers %in% rownames(seurat@assays$RNA)]
    markers <- markers[markers %in% rownames(seurat@assays$RNA)]
    print(str_interp("Not all markers occur in the data matrix of the seurat object. Specifically, ${mising_genes} are missing. Removing it and continuing"))
  }
  
  # get counts and perform gene wise scaling
  cnts_scaled <- as.matrix(seurat@assays$RNA[markers, ordered_index]) %>% 
    t() %>% 
    scale() %>% 
    t()
  
  
  
  if(!is.null(cap_value)){
    message(str_interp("Zscores > |${cap_value}| are set to ${cap_value} (or -${cap_value})"))
    changed_values <- sum(cnts_scaled>cap_value) + sum(cnts_scaled < -cap_value)
    cnts_scaled[cnts_scaled>cap_value] = cap_value 
    cnts_scaled[cnts_scaled < -cap_value] = -cap_value 
    message(str_interp("This was the case for ${changed_values} values"))
    col_fun = colorRamp2(c(-cap_value, 0, cap_value), c("blue", "white", "red"))
    
  } else {col_fun = colorRamp2(breaks=c(min(cnts_scaled, na.rm=T), mean(cnts_scaled, na.rm=T), max(cnts_scaled, na.rm=T)),
                     colors=c("blue", "white", "red"))}

  p <- Heatmap(cnts_scaled,
               cluster_rows = F, 
               cluster_columns = F, 
               show_row_names = T, 
               show_column_names = F,
               column_split = c_split,
               row_split = celltype,
               name="Z-score",
               col=col_fun,
               row_gap = unit(4, "mm"),
               border = "black",
               column_title_gp = grid::gpar(fontsize = 9, fontface="bold"),
               column_title_rot = 90,
               row_names_gp = grid::gpar(fontsize = 10, fontface="bold"),
               use_raster = F,
               heatmap_legend_param = list(at = c(-cap_value, 0, cap_value)))
  p <- ggplotify::as.ggplot(p)
  
  return(p)
  
}

# runs the "normal" seurat pipeline
run_seurat_steps <- function(seurat_object, include_leiden=F, include_norm=T, include_tsne=F, include_louvain=F){
  
  if (include_norm){
    print("Running normalization and variable feature extraction")
    seurat_object <- seurat_object %>% 
      NormalizeData(.,  normalization.method = "LogNormalize", scale.factor = 10000) %>% 
      FindVariableFeatures(., selection.method = "vst", nfeatures = 3000) 
  }
  
  print("Running PCA and UMAP")
  seurat_object <- seurat_object %>% 
    ScaleData(.) %>% 
    RunPCA(., npcs = 50, verbose = T) %>% 
    RunUMAP(., reduction = "pca", dims = 1:20, verbose = T) 
  
  if(include_tsne){
    print("Running Tsne")
    seurat_object <- seurat_object %>% 
      RunTSNE(., reduction = "pca", dims = 1:20, verbose = T) 
  }
    
  
  if (include_leiden){
    print("Running neighborhood detection")
    seurat_object <- seurat_object %>% 
      FindNeighbors(., dims = 1:20, k.param = 10, verbose = F) %>% 
      FindClusters(., 
                   algorithm = 4, 
                   resolution = c(0.3, 0.5,
                                  0.6, 0.7), 
                   method="igraph") 
    
  }
  if (include_louvain){
    print("Running neighborhood detection")
    seurat_object <- seurat_object %>% 
      FindNeighbors(., dims = 1:20, k.param = 10, verbose = F) %>% 
      FindClusters(., 
                   resolution = c(0.3, 0.5,
                                  0.6, 0.7)) 
    
  }
  
  return(seurat_object)
}











