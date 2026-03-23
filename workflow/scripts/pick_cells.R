cat("*** Loading libraries*** \n")
library(argparse)
library(ggplot2)
library(funr)
library(patchwork)
library(mclust)
library(rtracklayer)
library(ggrepel)
library(scales)




# Source aux functions
source(paste0(dirname(funr::sys.script()),"/func.R"))
set.seed(1234)

########### Arguments parser

parser <- ArgumentParser()

parser$add_argument("-s", "--sample", type="character", default='foo', 
                    help="sample name [as in config file key]")

parser$add_argument("-m", "--modality   ", type="character", default='foo',
                    help="antibody name [as in config file key]")

parser$add_argument("-o", "--out_prefix", type="character", default="/",
                    help="folder for the output in clustering_snakemake folder")

parser$add_argument("--bcd_all", type="character",
                    help="Path to barcodes summary statistics per cell")

parser$add_argument("--bcd_peak", type="character",
                    help="Path to barcodes in peaks summary statistics per cell")

parser$add_argument("--metadata", type="character",
                    help="Path to the cellranger metadata singlecell.csv file")

parser$add_argument("--fragments", type="character",
                    help="Path to the cellranger fragments.tsv.gz file")

parser$add_argument("--min_reads", type="double", default='3.0',
                    help="Minimum number of reads per cell")

parser$add_argument("--max_reads", type="double", default='5.5',
                    help="Maximum number of reads per cell")

parser$add_argument("--peak_fraction_min", type="double", default='0.2',
                    help="Minimum fraction of reads within peak")

parser$add_argument("--peak_fraction_max", type="double", default='1',
                    help="Minimum fraction of reads within peak")



args      <- parser$parse_args()
# saveRDS(object=args,file='results/arguments.Rds')

cutoff_reads_min            = args$min_reads
cutoff_reads_max            = args$max_reads
cutoff_peak_percentage_low  = args$peak_fraction_min
cutoff_peak_percentage_high = args$peak_fraction_max

########## Filter the barcodes
cat("*** Reading barcode statistics files \n")

all_barcodes_file  <- args$bcd_all
peak_barcodes_file <- args$bcd_peak
metadata_file      <- args$metadata

metadata = read.csv(metadata_file, header = 1)
metadata = metadata[2:nrow(metadata),]
metadata$logUMI = log10(metadata$passed_filters + 1)
metadata$promoter_ratio = (metadata$promoter_region_fragments+1) / (metadata$passed_filters + 1)
metadata$peak_region_ratio = (metadata$peak_region_fragments+1) / (metadata$passed_filters + 1)

# Fix metadata cell barcodes
# metadata$barcode <- paste0(args$sample,"_",metadata$barcode)

# Read barcode statistics files
all_barcodes <- read.table(file=all_barcodes_file)
peak_barcodes <- read.table(file=peak_barcodes_file)
bcd <- merge(all_barcodes,peak_barcodes,by="V2")  
colnames(bcd) <- c("barcode","all_unique_MB","peak_MB")
bcd$peak_ratio_MB <- bcd$peak_MB/bcd$all_unique_MB
bcd$sample <- args$sample


# Merge 10x metadata with barcode statistics
metadata <- merge(metadata,bcd,by='barcode')


metadata$is__cell_barcode <- as.factor(metadata$is__cell_barcode)

################ MB filtering
cat("*** Filtering cells \n")

# First hard-coded filtration for droplets with super little signal and <0.1 fraction in peaks
metadata <- metadata[metadata$all_unique_MB > 50 & metadata$peak_ratio_MB > 0.1,]

# Legacy filtering
metadata[,"passedMB"] <- FALSE
metadata[metadata$all_unique_MB > 10^cutoff_reads_min &
         metadata$all_unique_MB < 10^cutoff_reads_max &
         metadata$peak_ratio_MB > cutoff_peak_percentage_low &
         metadata$peak_ratio_MB < cutoff_peak_percentage_high,"passedMB"] <- TRUE


metadata$passedMB_legacy <- metadata$passedMB
################ Metadata picking by fiting GMM

metadata.pick <- pick_cells_MB(log_unique_reads = log10(metadata$all_unique_MB), frip = metadata$peak_ratio_MB)

metadata$class    <- metadata.pick$class
metadata$passedMB <- metadata.pick$pass_model


################ Cell picking scatterplot nreads ~ percent in peaks
dir.create(args$out_prefix,recursive=TRUE)

safe_density2d <- function(p, color = "black") {
  tryCatch(
    {
      p2 <- p + geom_density2d(col = color)
      ggplot_build(p2)
      p2
    },
    error = function(e) {
      message("Skipping geom_density2d(): ", conditionMessage(e))
      p
    }
  )
}

label_df <- metadata %>%
  filter(!is.na(class)) %>%
  group_by(class) %>%
  summarise(
    n_cells = n(),
    mean_fragments = mean(all_unique_MB, na.rm = TRUE),
    mean_frip = mean(peak_ratio_MB, na.rm = TRUE),
    x = mean(log10(all_unique_MB), na.rm = TRUE),
    y = mean(peak_ratio_MB, na.rm = TRUE),
    pass = unique(passedMB)[1],
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "cluster ", class,
      "\nstatus: ", ifelse(pass, "passed cell", "empty droplets"),
      "\nN = ", n_cells,
      "\nmean fragments = ", round(mean_fragments, 0),
      "\nmean FRiP = ", round(mean_frip, 3)
    )
  )

# same palette for points and passed labels
class_levels <- sort(unique(na.omit(metadata$class)))
class_cols <- setNames(hue_pal()(length(class_levels)), class_levels)

# passed labels keep class color, failed labels become grey
label_df$label_color <- ifelse(
  label_df$pass,
  class_cols[as.character(label_df$class)],
  "grey55"
)

p0 <- ggplot(data = metadata) +
  geom_point(
    aes(x = log10(all_unique_MB), y = peak_ratio_MB, color = factor(class)),
    size = 0.5, alpha = 0.2
  ) +
  geom_label_repel(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    color = label_df$label_color,
    fill = alpha("white", 0.85),
    size = 3,
    label.size = 0.25,
    box.padding = 0.5,
    point.padding = 0.5,
    force = 2,
    max.overlaps = Inf,
    segment.color = label_df$label_color
  ) +
  scale_color_manual(values = class_cols, name = "class") +
  labs(
    x = "log10(all_unique_MB)",
    y = "FRiP"
  ) +
  theme_bw()

ggsave(
  plot = p0,
  filename = paste0(args$out_prefix, "cells_clustering_for_picking.png"),
  width = 14, height = 14, units = "in"
)

n_picked <- sum(as.numeric(as.character(metadata$is__cell_barcode)), na.rm = TRUE)

cellranger_label <- paste0("Picked: ", n_picked, " cells")

p1 <- ggplot(data = metadata, aes(x = log10(all_unique_MB), y = peak_ratio_MB)) +
  geom_point(aes(col = is__cell_barcode), size = 0.1) +
  annotate(
    "label",
    x = Inf, y = Inf,
    label = cellranger_label,
    hjust = 1.05, vjust = 1.1,
    size = 6
  ) +
  labs(title = "Default Cell Ranger cell picking") +
  theme(
    legend.position = "bottom",
    text = element_text(size = 26),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p2 <- ggplot(data = metadata, aes(x = log10(passed_filters), y = peak_region_fragments / passed_filters)) +
  geom_point(aes(col = is__cell_barcode), size = 0.1) +
  annotate(
    "label",
    x = Inf, y = Inf,
    label = cellranger_label,
    hjust = 1.05, vjust = 1.1,
    size = 6
  ) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 26),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p1 <- safe_density2d(p1)
p2 <- safe_density2d(p2)

ggsave(
  plot = p1 + p2,
  filename = paste0(args$out_prefix, "cells_10x.png"),
  width = 20, height = 10, units = "in"
)
n_picked <- sum(metadata$passedMB, na.rm = TRUE)
label_text <- paste0("Picked: ", n_picked, " cells")

p1 <- ggplot(data = metadata, aes(x = log10(all_unique_MB), y = peak_ratio_MB)) +
  geom_point(aes(col = passedMB), size = 0.1) +
  annotate(
    "label",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1.05, vjust = 1.1,
    size = 6
  ) +
  labs(title = "Nanoscope cell picking") +
  theme(legend.position = "bottom", text = element_text(size = 26))

p2 <- ggplot(data = metadata, aes(x = log10(passed_filters), y = peak_region_fragments / passed_filters)) +
  geom_point(aes(col = passedMB), size = 0.1) +
  annotate(
    "label",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1.05, vjust = 1.1,
    size = 6
  ) +
  theme(legend.position = "bottom", text = element_text(size = 26))

p1 <- safe_density2d(p1)
p2 <- safe_density2d(p2)

ggsave(plot = p1+p2,
       filename=paste0(args$out_prefix,'cells_picked.png'),width = 20,height = 10,units = 'in')


# ################# Export bw selected / unselected
# cat("*** Reading fragments file \n")
#
# fragments         <- args$fragments
# fragments_gr      <- rtracklayer::import(fragments,format = "bed")
#
#
# cat("*** Exporting merged bw files \n")
# barcode_pass   <- metadata$barcode[metadata$passedMB]
# barcode_nopass <- metadata$barcode[!metadata$passedMB]
#
# fragments.pass   <- fragments_gr[fragments_gr$name %in% barcode_pass]
# fragments.nopass <- fragments_gr[fragments_gr$name %in% barcode_nopass]
#
# cat("*** Calculating coverage \n")
# coverage.pass   <- GenomicRanges::coverage(fragments.pass)
# coverage.nopass <- GenomicRanges::coverage(fragments.nopass)
#
# cat("*** Normalizing \n")
# coverage.pass <- coverage.pass/length(fragments.pass)
# coverage.nopass <- coverage.nopass/length(fragments.nopass)
#
# cat("*** Exporting \n")
# dir.create(args$out_prefix)
# rtracklayer::export(object=coverage.pass,  con = paste0(args$out_prefix,'/cells_picked.bw'))
# rtracklayer::export(object=coverage.nopass,con = paste0(args$out_prefix,'/cells_not_picked.bw'))

cat("*** Writing metadata \n")
write.csv(x = metadata,file = paste0(args$out_prefix,'metadata.csv'))
