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
plotCounts <- function(obj,quantiles,feature,ylabel=feature) {
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
    ylab(ylabel) +
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
    guides(color = guide_legend(override.aes = list(size=5),title="")) +
    geom_text(data=coords,aes(label=modality,x=mid_point,y=max(umap.embeddings.merge$UMAP_2+.1*umap.embeddings.merge$UMAP_2)),
              colour='black', fontface = "bold",size=4)
  return(plot)
  
}


# Function to plot upset on peaks of each modality
getUpsetPeaks <- function(modalities, samples, combined_peaks_ls, input_ls) {
  
  list_up <- list()
  for (mod in modalities) {
    i <- 0
    for (smp in samples) {
      i <- i + 1
      
      # overlap
      combined_mod <- as.data.frame(combined_peaks_ls[[mod]])
      overlap <- GenomicRanges::findOverlaps( toGRanges(combined_mod), input_ls[[paste0(mod,"_",smp)]], select = "first")
      # convert NA to 0s and make it binary
      overlap[is.na(overlap)] = 0
      overlap <- ifelse(overlap>0,1,0)
      combined_mod[,smp] <- overlap
      
      # add new values
      if (i==1) {
        final <- combined_mod
      } else {
        final[,smp] <- overlap
      }
      
      # if last sample, append to list
      if (smp==samples[length(samples)]) {
        list_up[[mod]] <- final
      }
      
      
    }
  }
  
  if (length(list_up)==1) {
    pfinal=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
                 set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
                 base_annotations=list('Size'=(intersection_size(counts=FALSE))))  + xlab(names(list_up)[1]) + theme(axis.title.x = element_text(size=12))
  } else if (length(list_up)==2) {
    p1=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[1])+ theme(axis.title.x = element_text(size=12)) 
    p2=upset(list_up[[2]][,6:ncol(list_up[[1]])],colnames(list_up[[2]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[2]) + theme(axis.title.x = element_text(size=12))
    pfinal=ggarrange(p1,p2,ncol=2)
  } else if (length(list_up)==3) {
    p1=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE))))  + xlab(names(list_up)[1]) + theme(axis.title.x = element_text(size=12))
    p2=upset(list_up[[2]][,6:ncol(list_up[[1]])],colnames(list_up[[2]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[2]) + theme(axis.title.x = element_text(size=12))
    p3=upset(list_up[[3]][,6:ncol(list_up[[1]])],colnames(list_up[[3]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[3]) + theme(axis.title.x = element_text(size=12))
    pfinal=ggarrange(p1,p2,p3,ncol=3)
  } else {
    warning("Upset plot function implemented in this vignette is not implemented to work on more than 3 modalities")
  }
  return(pfinal)
}


plotPassed <- function(mdata.list,xaxis_text=NULL,angle_x=NULL) {
  if (is.null(xaxis_text)) { xaxis_text=12 }
  if (is.null(angle_x)) { angle_x=30 }
  i <- 0
  for (exp in names(mdata.list)) {
    i <- i + 1
    df <- mdata.list[[exp]]
    counts <- as.data.frame(table(df$passedMB))
    counts$Perc <- counts$Freq / sum(counts$Freq) * 100
    counts$sample <- exp
    # bind
    if (i==1) {
      metadata.df <- counts
    } else {
      metadata.df <- rbind(metadata.df,counts)
    }
  }
  names(metadata.df)[1] <- "passedMB"
  # now plot
  pl=ggplot(metadata.df,aes(x=sample,y=Freq,fill=passedMB)) +
    theme_bw() +
    geom_bar(stat="identity",color='black',alpha=.7) +
    # x-axis
    theme(axis.text.x=element_text(angle=angle_x, hjust=1, size = xaxis_text)) + 
    theme(axis.title.x = element_text(size= 0)) +
    xlab("") +
    theme(axis.title.x = element_text(size= 14)) +
    # y-axis
    ylab("Number of cells") +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size= 14)) 
  return(pl)
}


plotPassedCells <- function(obj,sample_list,mod_list) {
  # convert obj to df
  i <- 0
  for (nm in names(obj)) {
    i <- i + 1
    # get samples and barcode name
    for (s in sample_list) {
      if (grepl(s,nm)) { 
        sample <- s
      }
    }
    for (m in mod_list) {
      if (grepl(m,nm)) { 
        mod <- m
      }
    }
    obj[[nm]]$sample <- sample
    obj[[nm]]$modality <- mod
    if (i==1) {
      df <- obj[[nm]]
    } else {
      df <- rbind(df,obj[[nm]])
    }
  }
  plot <- ggplot(df,aes(x=all_unique_MB,y=peak_ratio_MB,fill=passedMB)) + 
    theme_bw() +
    geom_point(shape=21,size=.4) +
    scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
    #coord_cartesian(ylim = c(0,1),xlim = c(10,1000000)) +
    facet_grid(modality~sample) +
    theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0, face= 'bold')) +
    theme(strip.text.y = element_text(size = 11, colour = "black", face= 'bold')) + 
    scale_fill_manual(values = c("#F8766D","#00BA38")) +
    # x-axis
    ylab("Fractio reads in peaks") +
    theme(axis.text.x=element_text(angle=0, hjust=.5, size = 12)) + 
    theme(axis.title.x = element_text(size= 14)) +
    # y-axis
    xlab("UMI") +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size= 14))+
    guides(fill = guide_legend(override.aes = list(size=5)))
  plot <- ggarrange(plot,legend='bottom')
  return(plot)
}




commonCellHistonMarks <- function(mod1,name_mod1,mod2,name_mod2,mod3=NULL,name_mod3=NULL,sample) {
  
  # for 2 modalities
  x <- list(name_mod1=mod1$barcode,
            name_mod2=mod2$barcode)
  names(x) <- c(name_mod1,name_mod2)
  pp=ggVennDiagram(x) + 
    scale_fill_gradient(low = "white", high = "white") +
    theme(legend.position = "none") +
    ggtitle(sample) +
    theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  
  # 3 modalities
  if (!is.null(mod3)) {
    x <- list(name_mod1=mod1$barcode,
              name_mod2=mod2$barcode,
              name_mod3=mod3$barcode)
    names(x) <- c(name_mod1,name_mod2,name_mod3)
    pp=ggVennDiagram(x) + 
      scale_fill_gradient(low = "white", high = "white") +
      theme(legend.position = "none") +
      ggtitle(sample) +
      theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  }
  return(pp)
  
}


