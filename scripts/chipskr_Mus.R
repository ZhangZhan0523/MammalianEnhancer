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
library(optparse)

option_list <- list(
    # make_option(c("-s", "--species"), type = "character", default = NULL, action = "store", help = "species name"),
    make_option(c("-t", "--tissue"), type = "character", default = NULL, action = "store", help = "tissue name")
    # make_option(c("-e", "--enhancer"), type = "character", default = NULL, action = "store", help = "enhancer peak file path"),
    # make_option(c("-p", "--promoter"), type = "character", default = NULL, action = "store", help = "promoter peak file path"),
    # make_option(c("-a", "--actProm"), type = "character", default = NULL, action = "store", help = "active promoter peak file path"),
    # make_option(c("-g", "--gff"), type = "character", default = NULL, action = "store", help = "genome annotation file, can be gff3 or gtf"),
    # make_option(c("-w", "--workdir"), type = "character", default = NULL, action = "store", help = "working directory")
)
opt_parser <- parse_args(OptionParser(option_list = option_list, usage = "%prog [options]"))
print(opt_parser)
# species = opt_parser$species
tissue <- opt_parser$tissue
# enhancer_peak = opt_parser$enhancer
# promoter_peak = opt_parser$promoter
# actProm_peak = opt_parser$actProm
# gff = opt_parser$gff
# workdir = opt_parser$workdir
species <- "Mus_musculus"
enhancer_peak <- paste0("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_", tissue, "__overlap_optimal_enhancer.narrowPeak")
promoter_peak <- paste0("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_", tissue, "__overlap_optimal_promoter.narrowPeak")
actProm_peak <- paste0("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_", tissue, "__overlap_optimal_actProm.narrowPeak")
workdir <- "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/"

library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(GenomeInfoDb)
library(ggimage)
library(AnnotationHub)

setwd(workdir)
common_enhancer <- readPeakFile("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/common/Mus_musculus_Liver_Kidney_Brain_enhancer_common_liver.bed")
common_promoter <- readPeakFile("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/common/Mus_musculus_Liver_Kidney_Brain_promoter_common_liver.bed")
# if (species=='Hipposideros_larvatus'){
#     txdb = makeTxDbFromGFF(gff, format='gff3', organism = gsub('_', ' ', species))
# }else{
#     txdb = makeTxDbFromGFF(gff, organism = gsub('_', ' ', species))
# }
txdb <- loadDb("/media//Data/zhangz//chip/genomes/Mus_musculus/mm39_txdb_AH84139.sqlite")
chr2acc <- read.table("/media/Data/zhangz/chip/genomes/Mus_musculus/GRCm39chr2acc.txt", header = T)
new_seqnames <- chr2acc$Chromosome
names(new_seqnames) <- chr2acc$Accession.version
enhancer <- readPeakFile(enhancer_peak)
promoter <- readPeakFile(promoter_peak)
actProm <- readPeakFile(actProm_peak)
Peaks <- GRangesList(enhancer = enhancer, promoter = promoter, actProm = actProm)
Peaks <- lapply(Peaks, renameSeqlevels, new_seqnames)
coverage <- covplot(Peaks, weightCol = "V5") + facet_grid(chr ~ .id)
# main = paste("Coverage plot of enhancer, promoter and active promoter of", species, tissue, sep = " ")) +
annotate = lapply(Peaks, annotatePeak,
    TxDb = txdb, addFlankGeneInfo = TRUE,
    tssRegion = c(-1000, 1000), flankDistance = 5000
)
commons <- GRangesList(enhancer = common_enhancer, promoter = common_promoter)
commons <- lapply(commons, renameSeqlevels, new_seqnames)
annotate_common <- lapply(commons, annotatePeak,
    TxDb = txdb, addFlankGeneInfo = TRUE,
    tssRegion = c(-1000, 1000), flankDistance = 5000
)
common_enhancer_df <- as.data.frame(annotate_common$enhancer)
common_promoter_df <- as.data.frame(annotate_common$promoter)
common_enhancer_df$geneId <- mapIds(org.Mm.eg.db, keys = as.character(common_enhancer_df$geneId), column = "SYMBOL", keytype = "ENTREZID")
common_promoter_df$geneId <- mapIds(org.Mm.eg.db, keys = as.character(common_promoter_df$geneId), column = "SYMBOL", keytype = "ENTREZID")
setwd("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/common")
common_genes <- lapply(commons, seq2gene, c(-1000, 1000), 3000, txdb)
common_genes <- lapply(common_genes, function(i) {
    mapIds(org.Mm.eg.db, keys = as.character(i), column = "SYMBOL", keytype = "ENTREZID")
})
write.csv(common_enhancer_df, file = paste("Mus_musculus", "common_enhancer.csv", sep = "_"), row.names = FALSE)
write.csv(common_promoter_df, file = paste("Mus_musculus", "common_promoter.csv", sep = "_"), row.names = FALSE)
write.csv(common_genes$enhancer, file = paste("Mus_musculus", "common_enhancer_genes.csv", sep = "_"), row.names = FALSE)
write.csv(common_genes$promoter, file = paste("Mus_musculus", "common_promoter_genes.csv", sep = "_"), row.names = FALSE)
cc <- compareCluster(common_genes, fun = "enrichKEGG", pvalueCutoff = 0.05, pAdjustMethod = "BH", organism = "mmu")
png(file = paste("Mus_musculus", "tissue_common_kegg.png", sep = "_"), width = 6000, height = 12000, res = 600)
keggdot <- dotplot(cc, showCategory = 10)
dev.off()

enhancer_GO <- enrichGO(common_genes$enhancer, OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05, pAdjustMethod = "BH")
png(file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks/common/Mus_musculus_common_enhancer_GO.png", width = 6000, height = 12000, res = 600)
dotplot(enhancer_GO, showCategory = 10)
dev.off()
promoter_GO <- enrichGO(common_genes$promoter, OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05, pAdjustMethod = "BH")
png(file = "/media/Data/zhangz/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_common_promoter_GO.png", width = 6000, height = 12000, res = 600)
# plotAvgProf(promoter_GO, showCategory = 10)
dotplot(promoter_GO)
dev.off()
ccGO <- compareCluster(common_genes, fun = "enrichGO", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = org.Mm.eg.db)
png(file = "Mus_musculus_tissue_common_GO.png", width = 6000, height = 12000, res = 600)
dotplot(ccGO, showCategory = 10) + facet_grid(ONTOLOGY ~ .)
dev.off()
# enhancerHeat = peakHeatmap(enhancer, TxDb = txdb, upstream = 3000, downstream = 3000, color = rainbow(length(Peaks)))
# promoterHeat = peakHeatmap(promoter, TxDb = txdb, upstream = 3000, downstream = 3000, color = rainbow(length(Peaks)))
# actPromHeat = peakHeatmap(actProm, TxDb = txdb, upstream = 3000, downstream = 3000, color = rainbow(length(Peaks)))
enhancer_df <- as.data.frame(annotate$enhancer)
promoter_df <- as.data.frame(annotate$promoter)
actProm_df <- as.data.frame(annotate$actProm)
png(file = paste(species, tissue, "heat.png", sep = "_"), width = 6000, height = 12000, res = 600)
heat <- lapply(Peaks, peakHeatmap, TxDb = txdb, upstream = 3000, downstream = 3000, color = rainbow(length(Peaks)))
dev.off()
# spInfo = read.csv('/media/Data/zhangz/chip/scripts/info/species.csv',header=T,row.names=5)
genes <- lapply(Peaks, function(i) seq2gene(i, c(-1000, 1000), 3000, txdb))
# ognsm = spInfo[species, 'organism']
# ognsm = 'mmu'
# if (ognsm == ''){
#     print('No organism information')
# }else{
#     cc = compareCluster(genes, fun = "enrichKEGG", pvalueCutoff = 0.05, pAdjustMethod = "BH", organism = ognsm)
#     png(file=paste(species, tissue, "specific_kegg.png",sep = '_'),width=6000,height=12000,res=600)
#     keggdot = plotAvgProf(cc, showCategory = 10)
#     dev.off()
# }
# cc = compareCluster(genes, fun = "enrichGO", pvalueCutoff = 0.05, pAdjustMethod = "BH")
write.csv(enhancer_df, file = paste(species, tissue, "enhancer.csv", sep = "_"), row.names = FALSE)
write.csv(promoter_df, file = paste(species, tissue, "promoter.csv", sep = "_"), row.names = FALSE)
write.csv(actProm_df, file = paste(species, tissue, "actProm.csv", sep = "_"), row.names = FALSE)
write.csv(genes$enhancer, file = paste(species, tissue, "enhancer_genes.csv", sep = "_"), row.names = FALSE)
write.csv(genes$promoter, file = paste(species, tissue, "promoter_genes.csv", sep = "_"), row.names = FALSE)
write.csv(genes$actProm, file = paste(species, tissue, "actProm_genes.csv", sep = "_"), row.names = FALSE)
# png(file=paste(species, tissue, "chrCov.png",sep = '_'),width=6000,height=12000,res=600)
# plot(coverage)
# dev.off()
# png(file=paste(species, tissue, "heat.png",sep = '_'),width=6000,height=12000,res=600)
# plot(heat)
# dev.off()

## chipseeker annotation with 100,000bp flanking region
## read in files
Brain_enhancer <- readPeakFile("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_Brain_enhancer_optimal.narrowPeak.bed4", header = FALSE)
Brain_promoter <- readPeakFile("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_Brain_promoter_optimal.narrowPeak.bed4", header = FALSE)
Kidney_enhancer <- readPeakFile("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_Kidney_enhancer_optimal.narrowPeak.bed4", header = FALSE)
Kidney_promoter <- readPeakFile("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_Kidney_promoter_optimal.narrowPeak.bed4", header = FALSE)
Liver_enhancer <- readPeakFile("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_Liver_enhancer_optimal.narrowPeak.bed4", header = FALSE)
Liver_promoter <- readPeakFile("/media/usb1/chip/analysis/Mus_musculus/anno/peaks/Mus_musculus_Liver_promoter_optimal.narrowPeak.bed4", header = FALSE)

enhancers <- GRangesList(Brain_enhancer = Brain_enhancer, Kidney_enhancer = Kidney_enhancer, Liver_enhancer = Liver_enhancer)
promoters <- GRangesList(Brain_promoter = Brain_promoter, Kidney_promoter = Kidney_promoter, Liver_promoter = Liver_promoter)
enhancers <- lapply(enhancers, renameSeqlevels, new_seqnames)
promoters <- lapply(promoters, renameSeqlevels, new_seqnames)
enhancer_anno <- lapply(enhancers, annotatePeak, TxDb = txdb, addFlankGeneInfo = TRUE, tssRegion = c(-10, 10), flankDistance = 100000)
promoter_annos <- lapply(promoters, annotatePeak, TxDb = txdb, addFlankGeneInfo = TRUE, tssRegion = c(-10, 10), flankDistance = 100000)
enhancer_Brain_anno_df <- as.data.frame(enhancer_anno$Brain_enhancer)
enhancer_Kidney_anno_df <- as.data.frame(enhancer_anno$Kidney_enhancer)
enhancer_Liver_anno_df <- as.data.frame(enhancer_anno$Liver_enhancer)
promoter_Brain_anno_df <- as.data.frame(promoter_annos$Brain_promoter)
promoter_Kidney_anno_df <- as.data.frame(promoter_annos$Kidney_promoter)
promoter_Liver_anno_df <- as.data.frame(promoter_annos$Liver_promoter)
write.csv(enhancer_Brain_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Brain_enhancer_anno.csv", row.names = FALSE)
write.csv(enhancer_Kidney_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Kidney_enhancer_anno.csv", row.names = FALSE)
write.csv(enhancer_Liver_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Liver_enhancer_anno.csv", row.names = FALSE)
write.csv(promoter_Brain_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Brain_promoter_anno.csv", row.names = FALSE)
write.csv(promoter_Kidney_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Kidney_promoter_anno.csv", row.names = FALSE)
write.csv(promoter_Liver_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Liver_promoter_anno.csv", row.names = FALSE)
enhancer_Brain_anno_df$geneId <- mapIds(org.Mm.eg.db, keys = as.character(enhancer_Brain_anno_df$geneId), column = "SYMBOL", keytype = "ENTREZID")
enhancer_Kidney_anno_df$geneId <- mapIds(org.Mm.eg.db, keys = as.character(enhancer_Kidney_anno_df$geneId), column = "SYMBOL", keytype = "ENTREZID")
enhancer_Liver_anno_df$geneId <- mapIds(org.Mm.eg.db, keys = as.character(enhancer_Liver_anno_df$geneId), column = "SYMBOL", keytype = "ENTREZID")
promoter_Brain_anno_df$geneId <- mapIds(org.Mm.eg.db, keys = as.character(promoter_Brain_anno_df$geneId), column = "SYMBOL", keytype = "ENTREZID")
promoter_Kidney_anno_df$geneId <- mapIds(org.Mm.eg.db, keys = as.character(promoter_Kidney_anno_df$geneId), column = "SYMBOL", keytype = "ENTREZID")
promoter_Liver_anno_df$geneId <- mapIds(org.Mm.eg.db, keys = as.character(promoter_Liver_anno_df$geneId), column = "SYMBOL", keytype = "ENTREZID")
write.csv(enhancer_Brain_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Brain_enhancer_anno_symbol.csv", row.names = FALSE)
write.csv(enhancer_Kidney_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Kidney_enhancer_anno_symbol.csv", row.names = FALSE)
write.csv(enhancer_Liver_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Liver_enhancer_anno_symbol.csv", row.names = FALSE)
write.csv(promoter_Brain_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Brain_promoter_anno_symbol.csv", row.names = FALSE)
write.csv(promoter_Kidney_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Kidney_promoter_anno_symbol.csv", row.names = FALSE)
write.csv(promoter_Liver_anno_df, file = "/media/usb1/chip/analysis/Mus_musculus/anno/peaks2/test2/Mus_musculus_Liver_promoter_anno_symbol.csv", row.names = FALSE)
