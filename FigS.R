library(plyr)
library(doMC)
registerDoMC(4)
library(ggplot2)
library(ggtree)
library(treeio)
library(aplot)
library(ggpubr)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(reshape2)
library(ape)
library(RRphylo) ## tree operating
library(ggunchained) ## half split
library(ggrepel)

# read in data
setwd("/media/Data/zhangz/chip/analysis/summary2/rep_ele")
species <- c(
    "Ovis_aries",
    "Bos_taurus", "Neophocaena_asiaeorientalis",
    "Sus_scrofa", "Lama_glama", "Mustela_putorius",
    "Canis_lupus", "Felis_catus", "Equus_asinus",
    "Equus_caballus", "Rhinolophus_pusillus", "Rhinolophus_ferrumequinum",
    "Hipposideros_larvatus", "Myotis_ricketti", "Myotis_chinensis",
    "Atelerix_albiventris", "Mus_musculus", "Rattus_norvegicus",
    "Cavia_porcellus", "Oryctolagus_cuniculus", "Macaca_mulatta",
    "Rhinopithecus_roxellana", "Tupaia_belangeri",
    "Procavia_capensis", "Petaurus_breviceps"
)

species <- species[!species %in% c(
    "Hipposideros_larvatus", "Myotis_ricketti", "Rhinolophus_ferrumequinum"
)]
rep_ele_class_enrichment_res <- adply(species, 1, function(sp) {
    if (sp == "Neophocaena_asiaeorientalis") {
        tiss <- c("Brain")
    } else if (sp == "Rhinopithecus_roxellana") {
        tiss <- c("Liver", "Kidney")
    } else {
        tiss <- c("Brain", "Kidney", "Liver")
    }
    read_rep <- function(tis, ele) read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/rep_ele/", sp, "_", tis, "_", ele, "_rep_ele_class_count.csv"), header = TRUE)
    print_res <- function(tis, ele) print(paste(tis, ele))
    res <- mdply(expand.grid(tis = tiss, ele = c("enhancer", "promoter")), read_rep)
    res$species <- sp
    return(res)
})
rep_class_box <- ggplot(rep_ele_class_enrichment_res, aes(x = repeat_class, y = ele_ratio, fill = tis)) +
    geom_boxplot(outliers = FALSE) +
    scale_fill_manual(values = c("Brain" = "#FFD2A8", "Kidney" = "#84C990", "Liver" = "#7E8BB4")) +
    # geom_point(data = rep_ele_class_enrichment_res[rep_ele_class_enrichment_res$species %in% c("Rhinolophus_pusillus", "Myotis_chinensis"), ], aes(x = repeat_class, y = ele_ratio, color = tis), position = position_jitter(width = 0.2), size = 2, alpha = 1, shape = 17) +
    facet_grid(ele ~ ., scales = "free_y") +
    scale_y_continuous(n.breaks = 3, labels = scales::percent) +
    labs(x = "Repeat class", y = "Proportion", fill = "Tissue") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "top")
ggsave("/media/Data/zhangz/chip/analysis/summary2/rep_ele/fig/rep_class_box.pdf", rep_class_box, width = 10, height = 10)


ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)
og_mark <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_mark.csv", header = TRUE)

# read tree
tree <- read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
tree <- keep.tip(tree, species)
# group_info <- split(species_info$species, species_info$order)
# tree <- groupOTU(tree, group_info)
# max_length <- max(tree$edge.length)
base_tre_horizontal <- ggtree(tree, layout = "rectangular") + coord_flip()

rep_ele_enrichment_res <- adply(species, 1, function(sp) {
    if (sp == "Neophocaena_asiaeorientalis") {
        tiss <- c("Brain")
    } else if (sp == "Rhinopithecus_roxellana") {
        tiss <- c("Liver", "Kidney")
    } else {
        tiss <- c("Brain", "Kidney", "Liver")
    }
    read_rep <- function(tis, ele) read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/rep_ele/", sp, "_", tis, "_", ele, "_rep_ele_count.csv"), header = TRUE)
    # print_res <- function(tis, ele) print(paste(tis, ele))
    res <- mdply(expand.grid(tis = tiss, ele = c("enhancer", "promoter")), read_rep)
    res$species <- sp
    return(res)
})

rep_ele_enrichment_res$repeat_class <- sapply(strsplit(rep_ele_enrichment_res$repeat_classfamily, "/"), "[", 1)
rep_ele_enrichment_sig <- rep_ele_enrichment_res[rep_ele_enrichment_res$p_value < 0.05, ]
rep_ele_enrichment_sig$repeat_class <- factor(rep_ele_enrichment_sig$repeat_class, levels = c("SINE", "LINE", "LTR", "DNA", "others"))
rep_ele_sele <- rep_ele_enrichment_sig[rep_ele_enrichment_sig$repeat_class %in% c("SINE", "LINE", "LTR", "DNA"), ]

class_family_order <- data.frame(
    repeat_classfamily = unique(rep_ele_sele$repeat_classfamily),
    repeat_class = sapply(strsplit(unique(rep_ele_sele$repeat_classfamily), "/"), "[", 1),
    mean_enrichment = sapply(unique(rep_ele_sele$repeat_classfamily), function(x) sum(rep_ele_sele[rep_ele_sele$repeat_classfamily == x, ]$enrichment / 25)),
    # mean_enh_ratio = sapply(unique(rep_ele_sele$repeat_classfamily), function(x) sum(rep_ele_sele[rep_ele_sele$repeat_classfamily == x, ]$enh_ratio / 75)),
    median_enrichment = sapply(unique(rep_ele_sele$repeat_classfamily), function(x) median(rep_ele_sele[rep_ele_sele$repeat_classfamily == x, ]$enrichment))
)
class_family_order$repeat_classfamily <- factor(class_family_order$repeat_classfamily, levels = class_family_order[order(class_family_order$repeat_class, -class_family_order$mean_enrichment), ]$repeat_classfamily)

class_family_order <- class_family_order[order(class_family_order$repeat_class, -class_family_order$mean_enrichment), ]


rep_enrich_point_enh <- ggplot(rep_ele_sele[rep_ele_sele$ele == "enhancer", ], aes(x = species, y = repeat_classfamily, size = -log10(q_value + 0.000001), color = log(enrichment))) +
    geom_point() +
    facet_grid(tis ~ .) +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(class_family_order$repeat_classfamily))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "#c01d2e", midpoint = 0) +
    labs(title = "Repeat element enrichment in enhancers", x = "Species", y = "Repeat element", size = "-log10(qvalue)", color = "log(Enrichment)") +
    theme(legend.position = "top", axis.title = element_blank())
rep_enh_p <- rep_enrich_point_enh %>% insert_bottom(base_tre_horizontal, height = 0.1)

rep_enrich_point_pro <- ggplot(rep_ele_sele[rep_ele_sele$ele == "promoter", ], aes(x = species, y = repeat_classfamily, size = -log10(q_value + 0.000001), color = log(enrichment))) +
    geom_point() +
    facet_grid(tis ~ ele) +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(class_family_order$repeat_classfamily))) +
    # scale_y_discrete(limits = rev(levels(rep_ele_enrichment_sig$repeat_class))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "#03324f", midpoint = 0) +
    labs(title = "Repeat element enrichment in promoters", x = "Species", y = "Repeat element", size = "-log10(qvalue)", color = "log(Enrichment)") +
    theme(legend.position = "top", axis.title = element_blank())
rep_pro_p <- rep_enrich_point_pro %>% insert_bottom(base_tre_horizontal, height = 0.1)

ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/rep_ele_enrichment_point_enh.pdf", rep_enh_p, width = 20, height = 30)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/rep_ele_enrichment_point_pro.pdf", rep_pro_p, width = 20, height = 30)

library(ggplotify)
species_func <- function(sp) {
    sp <- gsub("Macaca_mulatta", "Macaca_fascicularis", sp)
    sp <- gsub("_", " ", sp)
    return(sp)
}
rep_enrich_point_enh_Brain <- ggplot(rep_ele_sele[rep_ele_sele$ele == "enhancer" & rep_ele_sele$tis == "Brain", ], aes(x = species, y = repeat_classfamily, size = -log10(q_value + 0.000001), color = log(enrichment))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(class_family_order$repeat_classfamily))) +
    scale_x_discrete(labels = species_func) +
    theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "#c01d2e", midpoint = 0) +
    labs(title = "Repeat element enrichment in brain enhancers", size = "-log10(qvalue)", color = "log(Enrichment)") +
    theme(legend.position = "top", axis.title = element_blank()) +
    theme(title = element_text(size = 15))
base_tre_horizontal_B <- ggtree(drop.tip(tree, "Rhinopithecus_roxellana"), layout = "rectangular") + coord_flip()
rep_enh_p_Brain <- rep_enrich_point_enh_Brain %>% insert_bottom(base_tre_horizontal_B, height = 0.1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/rep_ele_enrichment_point_enh_Brain.pdf", rep_enh_p_Brain, width = 15, height = 10)

rep_enrich_point_enh_Kidney <- ggplot(rep_ele_sele[rep_ele_sele$ele == "enhancer" & rep_ele_sele$tis == "Kidney", ], aes(x = species, y = repeat_classfamily, size = -log10(q_value + 0.000001), color = log(enrichment))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(class_family_order$repeat_classfamily))) +
    scale_x_discrete(labels = species_func) +
    theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "#c01d2e", midpoint = 0) +
    labs(title = "Repeat element enrichment in kidney enhancers", size = "-log10(qvalue)", color = "log(Enrichment)") +
    theme(legend.position = "top", axis.title = element_blank()) +
    theme(title = element_text(size = 15))
base_tre_horizontal_KL <- ggtree(drop.tip(tree, "Neophocaena_asiaeorientalis"), layout = "rectangular") + coord_flip()
rep_enh_p_Kidney <- rep_enrich_point_enh_Kidney %>% insert_bottom(base_tre_horizontal_KL, height = 0.1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/rep_ele_enrichment_point_enh_Kidney.pdf", rep_enh_p_Kidney, width = 15, height = 10)

rep_enrich_point_enh_Liver <- ggplot(rep_ele_sele[rep_ele_sele$ele == "enhancer" & rep_ele_sele$tis == "Liver", ], aes(x = species, y = repeat_classfamily, size = -log10(q_value + 0.000001), color = log(enrichment))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(class_family_order$repeat_classfamily))) +
    scale_x_discrete(labels = species_func) +
    theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "#c01d2e", midpoint = 0) +
    labs(title = "Repeat element enrichment in liver enhancers", size = "-log10(qvalue)", color = "log(Enrichment)") +
    theme(legend.position = "top", axis.title = element_blank()) +
    theme(title = element_text(size = 15))
rep_enh_p_Liver <- rep_enrich_point_enh_Liver %>% insert_bottom(base_tre_horizontal_KL, height = 0.1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/rep_ele_enrichment_point_enh_Liver.pdf", rep_enh_p_Liver, width = 15, height = 10)

rep_enrich_point_pro_Brain <- ggplot(rep_ele_sele[rep_ele_sele$ele == "promoter" & rep_ele_sele$tis == "Brain", ], aes(x = species, y = repeat_classfamily, size = -log10(q_value + 0.000001), color = log(enrichment))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(class_family_order$repeat_classfamily))) +
    scale_x_discrete(labels = species_func) +
    theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "#03324f", midpoint = 0) +
    labs(title = "Repeat element enrichment in brain promoters", size = "-log10(qvalue)", color = "log(Enrichment)") +
    theme(legend.position = "top", axis.title = element_blank()) +
    theme(title = element_text(size = 15))
rep_pro_p_Brain <- rep_enrich_point_pro_Brain %>% insert_bottom(base_tre_horizontal_B, height = 0.1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/rep_ele_enrichment_point_pro_Brain.pdf", rep_pro_p_Brain, width = 15, height = 10)

rep_enrich_point_pro_Kidney <- ggplot(rep_ele_sele[rep_ele_sele$ele == "promoter" & rep_ele_sele$tis == "Kidney", ], aes(x = species, y = repeat_classfamily, size = -log10(q_value + 0.000001), color = log(enrichment))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(class_family_order$repeat_classfamily))) +
    scale_x_discrete(labels = species_func) +
    theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "#03324f", midpoint = 0) +
    labs(title = "Repeat element enrichment in kidney promoters", size = "-log10(qvalue)", color = "log(Enrichment)") +
    theme(legend.position = "top", axis.title = element_blank()) +
    theme(title = element_text(size = 15))
rep_pro_p_Kidney <- rep_enrich_point_pro_Kidney %>% insert_bottom(base_tre_horizontal_KL, height = 0.1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/rep_ele_enrichment_point_pro_Kidney.pdf", rep_pro_p_Kidney, width = 15, height = 10)

rep_enrich_point_pro_Liver <- ggplot(rep_ele_sele[rep_ele_sele$ele == "promoter" & rep_ele_sele$tis == "Liver", ], aes(x = species, y = repeat_classfamily, size = -log10(q_value + 0.000001), color = log(enrichment))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(class_family_order$repeat_classfamily))) +
    scale_x_discrete(labels = species_func) +
    theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "#03324f", midpoint = 0) +
    labs(title = "Repeat element enrichment in liver promoters", size = "-log10(qvalue)", color = "log(Enrichment)") +
    theme(legend.position = "top", axis.title = element_blank()) +
    theme(title = element_text(size = 15))
rep_pro_p_Liver <- rep_enrich_point_pro_Liver %>% insert_bottom(base_tre_horizontal_KL, height = 0.1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/rep_ele_enrichment_point_pro_Liver.pdf", rep_pro_p_Liver, width = 15, height = 10)


# ent_all <- adply(species, 1, function(sp) {
#     # read ent1 from file and concate together, no matter gene copy number
#     ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
#     res <- data.frame(t(ent[, c("OGID", "ent1")]))
#     colnames(res) <- res[1, ]
#     res <- res[-1, ]
#     rownames(res) <- sp
#     new_col <- data.frame(species = sp)
#     res <- cbind(new_col, res)
#     return(res)
# }, .parallel = TRUE)
# ent_all[is.na(ent_all)] <- 0
# ent_all <- ent_all[, -1]

ent_all <- adply(species, 1, function(sp) {
    # read ent1 from file and concate together, no matter gene copy number
    ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
    res <- data.frame(ent[which(ent$ent1 != 0), c("OGID", "ent1")])
    res$species <- sp
    return(res)
}, .parallel = TRUE)

species_func <- function(sp) {
    sp <- gsub("Macaca_mulatta", "Macaca_fascicularis", sp)
    sp <- gsub("_", " ", sp)
    return(sp)
}
ent_dis_plot <- ggplot(ent_all, aes(x = ent1)) +
    geom_density(fill = "#c6d068", linewidth = 0.3) +
    facet_wrap(~species, scales = "free", ncol = 3, labeller = as_labeller(species_func)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Regulatory entropy", y = "Density") +
    scale_y_continuous(n.breaks = 3) +
    scale_x_continuous(n.breaks = 3) +
    theme(strip.text = element_text(face = "bold.italic")) +
    theme(strip.background = element_rect(fill = "grey85", color = "transparent")) +
    theme(axis.line = element_line(linewidth = 0.25))
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/ent_dis.pdf", width = 9, height = 24)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/sup/ent_dis.png", width = 9, height = 24)


calc_ratio <- function(x, df) {
    x / sum(df$count)
}

genome_rep_count <- adply(species, 1, function(sp) {
    if (sp == "Macaca_mulatta") {
        sp_g <- "Macaca_fascicularis"
    } else {
        sp_g <- sp
    }
    gen_bg <- read.table(paste0("/media/Data/zhangz/chip/genomes/", sp_g, "/", sp_g, ".fa.out"), header = FALSE, comment.char = "*", skip = 2)
    colnames(gen_bg) <- c("score", "div", "del", "ins", "chr", "start", "end", "left", "strand", "rep_target", "rep_classfamily", "rep_start", "rep_end", "rep_left", "id")
    # split classfamily
    gen_bg$rep_classfamily <- as.character(gen_bg$rep_classfamily)
    gen_bg$rep_class <- sapply(strsplit(gen_bg$rep_classfamily, "/"), "[", 1)
    # calc count and ratio
    # gen_bg_count <- as.data.frame(table(gen_bg$rep_classfamily))
    # colnames(gen_bg_count) <- c("repeat_classfamily", "count")
    # gen_bg_count$ratio <- sapply(gen_bg_count$count, calc_ratio, gen_bg_count)
    # colnames(gen_bg_count) <- c("repeat_classfamily", "gen_bg_count", "gen_bg_ratio")

    # gen_bg_count_class <- as.data.frame(table(gen_bg$rep_class))
    # colnames(gen_bg_count_class) <- c("repeat_class", "count")
    res <- table(gen_bg$rep_classfamily) %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame()
    colnames(res) <- res[1, ]
    res <- res[-1, ]
    rownames(res) <- sp
    return(res)
}, .parallel = TRUE)

library(dplyr)

genome_rep_count <- data.frame(rep_classfamily = unique(gennome_rep_count$repeat_classfamily))
genome_rep_class_count <- data.frame(rep_class = colnames(gen_bg_count_class[, -1]))
for (sp in species) {
    if (sp == "Macaca_mulatta") {
        sp_g <- "Macaca_fascicularis"
    } else {
        sp_g <- sp
    }
    gen_bg <- read.table(paste0("/media/Data/zhangz/chip/genomes/", sp_g, "/", sp_g, ".fa.out"), header = FALSE, comment.char = "*", skip = 2)
    colnames(gen_bg) <- c("score", "div", "del", "ins", "chr", "start", "end", "left", "strand", "rep_target", "rep_classfamily", "rep_start", "rep_end", "rep_left", "id")
    # split classfamily
    gen_bg$rep_classfamily <- as.character(gen_bg$rep_classfamily)
    gen_bg$rep_class <- sapply(strsplit(gen_bg$rep_classfamily, "/"), "[", 1)
    res <- table(gen_bg$rep_classfamily) %>%
        as.data.frame()
    colnames(res) <- c("rep_classfamily", sp)
    genome_rep_count <- full_join(genome_rep_count, res, by = "rep_classfamily")
    res2 <- table(gen_bg$rep_class) %>%
        as.data.frame()
    colnames(res2) <- c("rep_class", sp)
    genome_rep_class_count <- full_join(genome_rep_class_count, res2, by = "rep_class")
}
genome_rep_count[is.na(genome_rep_count)] <- 0
genome_rep_class_count[is.na(genome_rep_class_count)] <- 0
write.csv(genome_rep_count, "/media/Data/zhangz/chip/analysis/summary2/tab/sup/genome_rep_count.csv", row.names = FALSE)
write.csv(genome_rep_class_count, "/media/Data/zhangz/chip/analysis/summary2/tab/sup/genome_rep_class_count.csv", row.names = FALSE)

gen_bg_count_class <- adply(species, 1, function(sp) {
    if (sp == "Macaca_mulatta") {
        sp_g <- "Macaca_fascicularis"
    } else {
        sp_g <- sp
    }
    gen_bg <- read.table(paste0("/media/Data/zhangz/chip/genomes/", sp_g, "/", sp_g, ".fa.out"), header = FALSE, comment.char = "*", skip = 2)
    colnames(gen_bg) <- c("score", "div", "del", "ins", "chr", "start", "end", "left", "strand", "rep_target", "rep_classfamily", "rep_start", "rep_end", "rep_left", "id")
    # split classfamily
    gen_bg$rep_classfamily <- as.character(gen_bg$rep_classfamily)
    gen_bg$rep_class <- sapply(strsplit(gen_bg$rep_classfamily, "/"), "[", 1)
    res <- table(gen_bg$rep_class) %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame()
    colnames(res) <- res[1, ]
    res <- res[-1, ]
    rownames(res) <- sp
    return(res)
    return(res)
}, .parallel = TRUE)
