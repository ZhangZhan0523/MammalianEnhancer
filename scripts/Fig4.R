## draw plots for Fig4, the theme is the regulatory entropies of different species

library(ggplot2)
library(reshape2)
library(UpSetR)

pgls_res <- read.csv("/media/Data/zhangz/chip/analysis/summary2/pgls/pgls_res.csv", header = TRUE)
pgls_res[["P-value"]] <- as.numeric(pgls_res[["P-value"]])
pgls_res[["Max_P-value"]] <- as.numeric(pgls_res[["Max_P-value"]])
pgls_res[["R2"]] <- as.numeric(pgls_res[["R2"]])
head(pgls_res)
summary(pgls_res)
head(pgls_res[order(pgls_res[["Max_P-value"]]), ])
nrow(pgls_res[pgls_res[["Max_P-value"]] < 0.05, ])
sig_pgls_res <- pgls_res[pgls_res[["Max_P-value"]] < 0.05, ]
summary(sig_pgls_res)
tapply(sig_pgls_res$R2, sig_pgls_res$trait, summary)
tapply(sig_pgls_res$R2, sig_pgls_res$trait, length)
sig_pgls_res <- sig_pgls_res[!is.na(sig_pgls_res$OGID), ]
pgls_res[!is.na(pgls_res$error), ]

## there are four ogid fail to fit PGLS model in Temperature.K.
## OG0009116 OG0010155 OG0011506 OG0003532
## only OG0010155 has a significant p-value in Temperature.K., with a R2 = 0.4054453
## the other three have a p-value > 0.05
## so we just ignore them

## plot the significant PGLS results
sig_pgls_res <- sig_pgls_res %>% arrange(trait, `P-value`)
tapply(sig_pgls_res$trait, sig_pgls_res$Model, length)
#  BM  OU
#   8 329
sig_pgls_res <- sig_pgls_res[, 3:11]
sig_pgls_res$gene <- sapply(sig_pgls_res$OGID, find_mus_ortho)
sig_pgls_res <- sig_pgls_res[!is.na(sig_pgls_res$Trait), ]
sig_pgls_res[sig_pgls_res$Trait == "MLres.y", ]
sig_pgls_res[sig_pgls_res$Trait %in% c("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y"), ]
sig_pgls_res[which(sig_pgls_res$Trait %in% c("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y") & sig_pgls_res$Slope > 0), ]
intersect(sig_pgls_res[sig_pgls_res$Trait == "ML.yrs.y", "OGID"], sig_pgls_res[sig_pgls_res$Trait == "MLres.y", "OGID"])
upset_matrix <- table(sig_pgls_res$OGID, sig_pgls_res$Trait)
upset_df <- as.data.frame(upset_matrix)
colnames(upset_df) <- c("OGID", "Trait", "Frequency")

# long df upset_df to wide df
wide_df <- spread(upset_df, Trait, Frequency)
p1 <- upset(wide_df,
    sets = c("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y", "AW.g"),
    nset = 23, nintersects = 60,
    queries = list(
        list(
            query = intersects, params = list("ML.yrs.y", "MLres.y"),
            color = "red", active = TRUE
        ),
        list(
            query = intersects, params = list("MLres.y", "FTMres.y"),
            color = "blue", active = TRUE
        ),
        list(
            query = intersects, params = list("ML.yrs.y", "FTM.d.y"),
            active = TRUE
        ),
        list(
            query = intersects, params = list("FTMres.y", "FTM.d.y"),
            active = TRUE
        )
    ),
)
pdf("/media/Data/zhangz/chip/analysis/summary2/pgls/upset_longevity.pdf", width = 10, height = 10)
print(p1)
dev.off()
p2 <- upset(wide_df,
    nset = 23, nintersects = 60,
    queries = list(
        # list(
        #     query = intersects, params = list("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y"),
        #     color = "red", active = TRUE
        # ),
        list(
            query = intersects, params = list("ML.yrs.y", "MLres.y"),
            active = TRUE
        ),
        list(
            query = intersects, params = list("MLres.y", "FTMres.y"),
            active = TRUE
        )
    ),
)
pdf("/media/Data/zhangz/chip/analysis/summary2/pgls/upset_all.pdf", width = 16, height = 9)
print(p2)
dev.off()

## plot the distribution of mouse entropy
sp <- "Mus_musculus"

ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
plot(density(ent$ent1))
ent1_dist <- ggplot(ent[ent$ent1 > 0, ], aes(x = ent1)) +
    geom_histogram(bins = 100, fill = "#324777", alpha = 0.5) +
    theme_minimal() +
    labs(x = "Gene regulatory entropy of mouse", y = "Count") +
    scale_x_continuous(breaks = c(0.25, 0.75, 1.25)) +
    scale_y_continuous(breaks = c(0, 200, 400)) +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8))
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/mouse_ent1_dist.pdf", ent1_dist, width = 8, height = 8, units = "cm")

## ecdf plot
ent_ecdf <- ecdf(ent[ent$ent1 > 0, ]$ent1)
ent_ecdf(0.95)
ent95 <- quantile(ent[ent$ent1 > 0, ]$ent1, 0.95)
ent1_ecdf <- ggplot(ent[ent$ent1 > 0, ], aes(x = ent1)) +
    stat_ecdf(geom = "step", color = "#324777") +
    # geom_segment(aes(x = 0, y = 0.95, xend = ent95, yend = 0.95), linetype = "dashed", color = "red") +
    # geom_segment(aes(x = ent95, y = 0, xend = ent95, yend = 0.95), linetype = "dashed", color = "red") +
    annotate("text", x = 0.25, y = 0.99, label = "95%", color = "#c01d2e", size = 6) +
    annotate("text", x = ent95 + 0.07, y = 0.05, label = round(ent95, 2), color = "#c01d2e", size = 6) +
    annotate("segment", x = 0, xend = ent95, y = 0.95, yend = 0.95, color = "#c01d2e", linetype = "dashed") +
    annotate("segment", x = ent95, xend = ent95, y = 0, yend = 0.95, color = "#c01d2e", linetype = "dashed") +
    theme_classic() +
    labs(x = "Gene regulatory entropy", y = "Cumulative probability", title = "ECDF plot of murine gene regulatory entropy") +
    scale_x_continuous(breaks = c(0.25, 0.75, 1.25), limits = c(0, 1.3), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.05), expand = c(0, 0)) +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/mouse_ent1_ecdf.pdf", ent1_ecdf, width = 6, height = 6)

# ent2_dist <- ggplot(ent[ent$ent2 > 0, ], aes(x = ent2)) +
#     geom_histogram(fill = "blue", alpha = 0.5) +
#     theme_minimal() +
#     labs(title = "Mouse Entropy Distribution", x = "Entropy", y = "Density")
# most ent2 = 1

## gene ontology analysis of high entropy genes in Mus_musculus
library(clusterProfiler)
library(org.Mm.eg.db)

library(KEGG.db)
setwd("/media/Data/zhangz/chip/analysis/summary2/ent_gene")

species <- c(
    "Ovis_aries",
    "Bos_taurus", "Neophocaena_asiaeorientalis",
    "Sus_scrofa", "Lama_glama", "Mustela_putorius",
    "Canis_lupus", "Felis_catus", "Equus_asinus",
    "Equus_caballus", "Rhinolophus_pusillus", "Myotis_chinensis",
    "Atelerix_albiventris", "Mus_musculus", "Rattus_norvegicus",
    "Cavia_porcellus", "Oryctolagus_cuniculus", "Macaca_mulatta",
    "Rhinopithecus_roxellana", "Tupaia_belangeri",
    "Procavia_capensis", "Petaurus_breviceps"
)

## load data
sp <- "Mus_musculus"
mus_ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
mus_ent <- mus_ent[mus_ent$ent1 > 0, ]
summary(mus_ent)
plot(density(mus_ent$ent1))
mus_ent <- mus_ent[order(mus_ent$ent1, decreasing = TRUE), ]

## get the top 100 high entropy genes
top100 <- mus_ent$gene[1:100]

## do GO and KEGG analysis
top100 <- unique(top100)
top100_GO <- enrichGO(top100, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
top100_id <- bitr(top100, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
top100_KEGG <- enrichKEGG(top100_id$ENTREZID,
    organism = "mmu",
    pvalueCutoff = 0.05, qvalueCutoff = 0.05,
    use_internal_data = TRUE
)
# plot top 10 GO terms
go_dot <- dotplot(top100_GO, showCategory = 15)
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/go_top100_dot.pdf", go_dot, width = 8, height = 9)
go_cnet_p <- cnetplot(top100_GO,
    categorySize = "pvalue",
    color.params = list(foldChange = "qvalue"), showCategory = 15
)
go_cnet_p
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/go_top100_cnet_p.pdf", go_cnet_p, width = 12, height = 9)

library(plyr)
library(magrittr)
doMC::registerDoMC(cores = 4)
## try to find the overlap of orthologs with high regulatory entropy in different species
high_ent_ogs <- adply(species, 1, function(sp) {
    ent_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
    ent_df <- ent_df[ent_df$ent1 > 0, ]
    ent_df <- ent_df[order(ent_df$ent1, decreasing = TRUE), ]

    ## get the top 100 high entropy genes
    top100_og <- ent_df[1:100, c("OGID", "gene")]

    ## do GO and KEGG analysis
    top100_og <- unique(top100_og)
    top100_og$species <- sp
    return(top100_og)
}, .parallel = TRUE)
# count occurrence of orthologs in top 100 high entropy genes
high_ent_ogs_count <- high_ent_ogs %>%
    count("OGID")
high_ent_ogs_count <- high_ent_ogs_count[order(high_ent_ogs_count$freq, decreasing = TRUE), ]
plot(density(high_ent_ogs_count$freq))
nrow(high_ent_ogs_count[high_ent_ogs_count$freq > 3, ])

# plot count bar plot of high_ent_ogs_count
## x axis is the number that high entropy orthologs occur in different species
## y axis is the frequency of orthologs
colnames(high_ent_ogs_count) <- c("OGID", "count")
high_ent_ogs_count$count <- as.numeric(high_ent_ogs_count$count)
species_count <- unique(high_ent_ogs_count) %>%
    count("count")
high_ent_fre_plot <- ggplot(species_count, aes(x = count, y = freq)) +
    geom_bar(stat = "identity", fill = "#324777", alpha = 0.5) +
    geom_text(aes(label = freq), vjust = -0.5, size = 5) +
    theme_classic() +
    labs(title = "The number of species in which the regulatory entropy of an orthogroup ranks among the top 100", y = "Count") +
    scale_x_continuous(breaks = c(1, 6, 11)) +
    scale_y_continuous(n.breaks = 3, limits = c(0, 1400), expand = c(0, 0)) +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), axis.title.x = element_blank())
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/orthologs_count.pdf", high_ent_fre_plot, width = 8, height = 8)

ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)
high_ent_ogs_count$mus_gene <- ortho_all$Mus_musculus[match(high_ent_ogs_count$OGID, ortho_all$Orthogroup)] %>%
    sapply(., function(x) {
        stringr::str_split(x, "[|]")[[1]][1]
    })
# high_ent_ogs_count$mus_id <- bitr(high_ent_ogs_count$mus_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
high_ent_id <- bitr(high_ent_ogs_count$mus_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
high_ent_ogs_count$mus_id <- high_ent_id$ENTREZID[match(high_ent_ogs_count$mus_gene, high_ent_id$SYMBOL)]
top48_mus_go <- enrichGO(high_ent_ogs_count$mus_gene[1:48],
    OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
    ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05
)
top48_mus_go_dot <- dotplot(top48_mus_go, showCategory = 15)
top48_mus_go_cnet_p <- cnetplot(top48_mus_go,
    categorySize = "pvalue",
    color.params = list(foldChange = "qvalue", edge = TRUE), showCategory = 13,
    layout = "graphopt", cex.params = list(gene_node = 0.5, category_label = 1.25, gene_label = 1)
)
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/overlap_top48_cnet.pdf", top48_mus_go_cnet_p, width = 13, height = 9)

top48_mus_kegg <- enrichKEGG(high_ent_ogs_count$mus_id[1:48],
    organism = "mmu",
    pvalueCutoff = 0.05, qvalueCutoff = 0.05,
    use_internal_data = TRUE
)

MLtraits <- c("ML.yrs.y", "MLres.y")
ML_ent_ogs <- c("Lactb2" = "OG0009639", "Abcf3" = "OG0007610", "Pafah2" = "OG0006213")

ent_all <- adply(species, 1, function(sp) {
    # read ent1 from file and concate together, no matter gene copy number
    ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
    res <- data.frame(t(ent[, c("OGID", "ent1")]))
    colnames(res) <- res[1, ]
    res <- res[-1, ]
    rownames(res) <- sp
    new_col <- data.frame(species = sp)
    res <- cbind(new_col, res)
    return(res)
}, .parallel = TRUE)
ent_all[is.na(ent_all)] <- 0
ent_all <- ent_all[, -1]
ent_all[, 2:ncol(ent_all)] <- apply(ent_all[, 2:ncol(ent_all)], 2, function(x) as.numeric(x))
head(ent_all[, 1:5])
ent_non0 <- data.frame(OGID = colnames(ent_all)[2:ncol(ent_all)], non0 = apply(ent_all[, 2:ncol(ent_all)], 2, function(x) sum(x != 0)))
nrow(ent_non0[ent_non0$non0 >= 15, ])
# 1425
ortho_15 <- ent_non0$OGID[ent_non0$non0 >= 15]
ent_15 <- ent_all[, c("species", ent_non0$OGID[ent_non0$non0 >= 15])]

ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)
og_mark <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_mark.csv", header = TRUE)


# read tree
tre <- read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
tre <- keep.tip(tre, species)

# read life history traits
traits <- read.csv("/media/Data/zhangz/chip/analysis/summary2/tissue_ent/traits.csv", header = TRUE)
traits <- traits[traits$species %in% species, c(1, 2, 12:42)]
traits$svg_path <- paste0("/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/", traits$species, ".svg")
traits$svg_path[which(traits$species %in% c("Myotis_chinensis", "Myotis_ricketti", "Rhinolophus_ferrumequinum", "Rhinolophus_pusillus"))] <- "/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/Hipposideros_larvatus.svg"
for (svg in traits$svg_path) {
    if (!file.exists(svg)) {
        print(svg)
    }
}

ent_trait <- merge(ent_15[, c("species", ML_ent_ogs)], traits[, c("species", MLtraits)], by = "species")
