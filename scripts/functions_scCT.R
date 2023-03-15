# convert list of objects to dataframe of metadata
listToMeta <- function(obj) {
  for (i in 1:length(obj)) {
    mdata <- obj[[i]][[]]
    if (i==1) {
      final <- mdata
    } else {
      final <- rbind(final,mdata)
    }
  }
  return(final)
}


# Plot UMI counts per each modality and sample
plotUMIcounts <- function(obj,quantiles,feature) {
  dataf <- listToMeta(obj)
  # plot
  pp=ggplot(dataf, aes(x="x",y=dataf[,feature],fill=sample)) +
    theme_bw() +
    geom_half_violin(draw_quantiles = .5) +
    geom_half_point(shape=1,aes(color=sample)) +
    facet_grid(sample~modality,scales = "free_x") +
    theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0, face= 'bold')) +
    theme(strip.text.y = element_text(size = 11, colour = "black", face= 'bold')) +
    # x axis
    xlab("Samples") +
    theme(axis.text.x=element_text(size = 0,angle = 0, hjust = .5)) +
    theme(axis.title.x = element_text(size=0)) +
    # y axis
    ylab(feature) +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size=14)) +
    theme(legend.position = "none") +
    # add quantiles
    stat_summary(fun = "quantile", fun.args = list(probs = quantiles), 
                 geom = "hline", aes(yintercept = ..y..), linetype = "dashed")
    return(pp)
}

# modified version of the Signac DepthCor - it just adds to the plot the modality name
DepthCorMulMod <- function(obj) {
  nm <- unique(obj[[]][,"modality"])
  p1 <- DepthCor(obj) +
    ggtitle(nm) +
    theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  return(p1)
}


# Umap with connected modalities 
plotConnectModal <- function(seurat,group) {
  nmodalities <- length(seurat)
  
  # First, get coords of UMAP, modality name, cluster name and cell barcode for each modality
  umap_embeddings <- list()
  i <- 0
  for(mod in names(seurat)) {
    i <- i + 1
    # get modality name without barcode
    mod2 <- strsplit(mod,"_")[[1]][1]
    # get UMAP 1 and 2 coords
    umap_embeddings[[mod2]]               <- as.data.frame(seurat[[mod]]@reductions[['umap']]@cell.embeddings)
    # adjust UMAP1 for plotting pursoses
    umap_embeddings[[mod2]]$UMAP_1        <- umap_embeddings[[mod2]]$UMAP_1 + (i-1)*40 
    # add modality name
    umap_embeddings[[mod2]]$modality      <- unique(seurat[[mod]]$modality) 
    # add cluster name
    umap_embeddings[[mod2]]$cluster       <- seurat[[mod]]@meta.data[,group] 
    # add cell barcode
    umap_embeddings[[mod2]]$cell_barcode  <- rownames(umap_embeddings[[mod2]]) 
  }
  
  # convert the list to dataf
  umap.embeddings.merge <- purrr::reduce(umap_embeddings,rbind)
  # get the names of the common cells among the analysed modalities
  common.cells                        <- table(umap.embeddings.merge$cell_barcode)
  common.cells                        <- names(common.cells[common.cells==nmodalities])
  # subset the umap_embeddings, selecting only the info from the common set of cells
  umap.embeddings.merge               <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]
  
  # label modality - get coords
  coords <- aggregate(umap.embeddings.merge$UMAP_1,by=list(umap.embeddings.merge$modality),max)
  names(coords) <- c("modality","max")
  coords$min <- aggregate(umap.embeddings.merge$UMAP_1,by=list(umap.embeddings.merge$modality),min)[,"x"]
  coords$mid_point <- (coords$max + coords$min) / 2
  
  plot <- ggplot(data=umap.embeddings.merge,aes(x=UMAP_1,y=UMAP_2,col=cluster)) + 
    geom_point(size=0.2) + 
    geom_line(data=umap.embeddings.merge, aes(group=cell_barcode,col=cluster),alpha=0.2,size=0.02) + 
    theme_classic() + NoAxes() +
    guides(color=guide_legend(title="")) +
    geom_text(data=coords,aes(label=modality,x=mid_point,y=max(umap.embeddings.merge$UMAP_2+.1*umap.embeddings.merge$UMAP_2)),
              colour='black', fontface = "bold",size=4)
  return(plot)
  
}





