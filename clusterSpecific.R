##annotate for tissue specific enhancers and promoters
## need to give arguments

# workdir="/media/Data/zhangz/usb1/chip/analysis/Rhinopithecus_roxellana/anno/peaks/Liverspecific"
# species="Rhinopithecus_roxellana"
# tissue="Liver"
# enhancer_peak="/media/Data/zhangz/usb1/chip/analysis/Rhinopithecus_roxellana/anno/peaks/Liverspecific/Rhinopithecus_roxellana_Liver__overlap_optimal_enhancer_specific.narrowPeak"
# promoter_peak="/media/Data/zhangz/usb1/chip/analysis/Rhinopithecus_roxellana/anno/peaks/Liverspecific/Rhinopithecus_roxellana_Liver__overlap_optimal_promoter_specific.narrowPeak"
# actProm_peak="/media/Data/zhangz/usb1/chip/analysis/Rhinopithecus_roxellana/anno/peaks/Liverspecific/Rhinopithecus_roxellana_Liver__overlap_optimal_actProm_specific.narrowPeak"
# gff="/media/Data/zhangz/chip/genomes/Rhinopithecus_roxellana/Rhinopithecus_roxellana.gtf"

library(optparse)

option_list <- list(
  make_option(c("-s", "--species"), type = "character", default = NULL, action = "store", help = "species name"),
  make_option(c("-t", "--tissue"), type = "character", default = NULL, action = "store", help = "tissue name"),
  make_option(c("-e", "--enhancer"), type = "character", default = NULL, action = "store", help = "enhancer peak file path"),
  make_option(c("-p", "--promoter"), type = "character", default = NULL, action = "store", help = "promoter peak file path"),
  make_option(c("-g", "--gff"), type = "character", default = NULL, action = "store", help = "genome annotation file, can be gff3 or gtf"),
  make_option(c("-w", "--workdir"), type = "character", default = NULL, action = "store", help = "working directory")
)
opt_parser <- parse_args(OptionParser(option_list = option_list, usage = "%prog [options]"))
print(opt_parser)
species <- opt_parser$species
tissue <- opt_parser$tissue
enhancer_peak <- opt_parser$enhancer
promoter_peak <- opt_parser$promoter
# actProm_peak <- opt_parser$actProm
gff <- opt_parser$gff
workdir <- opt_parser$workdir

library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(GenomeInfoDb)
library(ggimage)
library("AnnotationHub")
library(stringr)

setwd(workdir)

sp_info <- read.csv("/media/Data/zhangz/chip/scripts/info/species.csv", header = T, row.names = 11)
if (species == "Neophocaena_asiaeorientalis") {
  return("Neophocaena_asiaeorientalis have no enough tissue, exit")
  quit()
}
if (species == "Hipposideros_larvatus") {
  txdb <- makeTxDbFromGFF(gff, format = "gff3", organism = gsub("_", " ", species))
} else if (species == 'Mus_musculus') {
  txdb <- makeTxDbFromGFF('/media/Data/zhangz/chip/genomes/Mus_musculus/ncbi_dataset/data/GCF_000001635.27/genomic.gff',organism='Mus musculus')
} else {
  txdb <- makeTxDbFromGFF(gff, organism = gsub("_", " ", species))
}
enhancer <- readPeakFile(enhancer_peak)
promoter <- readPeakFile(promoter_peak)
# actProm <- readPeakFile(actProm_peak)
peaks <- GRangesList(enhancer = enhancer, promoter = promoter)
# coverage <- covplot(peaks, weightCol = "V5") + facet_grid(chr ~ .id)
# main = paste("Coverage plot of enhancer, promoter and active promoter of", species, tissue, sep = " ")) +
annotate <- lapply(peaks, annotatePeak,
  TxDb = txdb, addFlankGeneInfo = TRUE,
  tssRegion = c(-1000, 1000), flankDistance = 5000
)
# enhancerHeat = peakHeatmap(enhancer, TxDb = txdb, upstream = 3000, downstream = 3000, color = rainbow(length(peaks)))
# promoterHeat = peakHeatmap(promoter, TxDb = txdb, upstream = 3000, downstream = 3000, color = rainbow(length(peaks)))
# actPromHeat = peakHeatmap(actProm, TxDb = txdb, upstream = 3000, downstream = 3000, color = rainbow(length(peaks)))
enhancer_df <- as.data.frame(annotate$enhancer)
promoter_df <- as.data.frame(annotate$promoter)
# actProm_df <- as.data.frame(annotate$actProm)
png(file = paste(species, tissue, "heat.png", sep = "_"), width = 6000, height = 12000, res = 600)
heat <- lapply(peaks, peakHeatmap, TxDb = txdb, upstream = 3000, downstream = 3000) # , color = rainbow(length(peaks))
dev.off()

genes <- lapply(peaks, function(i) seq2gene(i, c(-1000, 1000), 3000, txdb))
write.csv(enhancer_df, file = paste(species, tissue, "enhancer.csv", sep = "_"), row.names = FALSE)
write.csv(promoter_df, file = paste(species, tissue, "promoter.csv", sep = "_"), row.names = FALSE)
# write.csv(actProm_df, file = paste(species, tissue, "actProm.csv", sep = "_"), row.names = FALSE)
write.csv(genes$enhancer, file = paste(species, tissue, "enhancer_genes.csv", sep = "_"), row.names = FALSE)
write.csv(genes$promoter, file = paste(species, tissue, "promoter_genes.csv", sep = "_"), row.names = FALSE)
ognsm <- sp_info[species, "organism"]
ognsmDB <- sp_info[species, "OrgDb"]

# write.csv(genes$actProm, file = paste(species, tissue, "actProm_genes.csv", sep = "_"), row.names = FALSE)
# png(file = paste(species, tissue, "specific_chrCov.png", sep = "_"), width = 6000, height = 12000, res = 600)
# plot(coverage)
# dev.off()
if (ognsm == "") {
  print("No organism information")
} else {
  cc <- compareCluster(genes, fun = "enrichKEGG", pvalueCutoff = 0.05, pAdjustMethod = "BH", organism = ognsm)
  ccGO <- compareCluster(genes, fun = "enrichGO", pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = loadDb(ognsmDB))
  png(file = paste(species, tissue, "specific_kegg.png", sep = "_"), width = 6000, height = 12000, res = 600)
  keggdot <- plotAvgProf(cc, showCategory = 10)
  dev.off()
}
# cc = compareCluster(genes, fun = "enrichGO", pvalueCutoff = 0.05, pAdjustMethod = "BH")



# png(file=paste(species, tissue, "heat.png",sep = '_'),width=6000,height=12000,res=600)
# plot(heat)
# dev.off()
# library(getopt)
# spec = matrix(data = c('species', 's', 2, 'character', 'species',
#                         'tissue', 't', 2, 'character', 'tissue',
#                         "enhancer", 'e', 2, 'character', 'enhancer peak file path',
#                         'promoter', 'p', 2, 'character', 'promoter peak file path',
#                         'actProm', 'ap', 2, 'character', 'active promoter peak file path',
#                         'gff', 'g', 2, 'character', 'genome annotation file, can be gff3 or gtf',
#                         'workdir', 'w', 2, 'character', 'working directory'),
#                         ncol = 5, byrow = TRUE)
# opt = getopt(spec)
# species = opt$species
# tissue = opt$tissue
# enhancer_peak = opt$enhancer
# promoter_peak = opt$promoter
# actProm_peak = opt$actProm
# gff = opt$gff
# workdir = opt$workdir
