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



load_ensembl_annot <- function(version = 'mm10') {
  library(ensembldb)
  
  if (version == 'mm10'){
    cat("*** Loading mm10 annotation of genes \n")
    library(EnsDb.Mmusculus.v79)
    ensdb = EnsDb.Mmusculus.v79
  } else if (version == 'hg38'){
    cat("*** Loading hg38 annotation of genes \n")
    library(EnsDb.Hsapiens.v86)
    ensdb = EnsDb.Hsapiens.v86
  } else {
    cat("*** ERROR: Only mm10 and hg38 supported\n")
    return(NULL)
  }
  
  seqlevelsStyle(ensdb) <- 'UCSC'
  
  gene.coords <- ensembldb::genes(ensdb, filter = ~ gene_biotype == "protein_coding")
  lncRNA.coords <- ensembldb::genes(ensdb, filter = ~ gene_biotype == "lincRNA")
  gene.coords <- c(gene.coords,lncRNA.coords)
  
  genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
  
  # Flatten the overlapping genes and extend by 2kb upstream of promoters
  genebody.coords.flat <- GenomicRanges::reduce(x = genebody.coords)
  genebodyandpromoter.coords.flat <- Signac::Extend(genebody.coords.flat,upstream = 2000)
  
  # Retrieve gene names from the original annotation (lost because of flatenning)
  genebodyandpromoter.coords.flat$name<- gene.coords[nearest(genebodyandpromoter.coords.flat,genebody.coords)]$gene_name
  return(genebodyandpromoter.coords.flat)
}
