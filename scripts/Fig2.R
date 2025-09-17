## draw plot about evolutionary pleiotropy and specificity

library(ggplot2)
library(plyr)
library(doMC)
library(reshape2)
library(scales)
library(stringr)
library(ggtree)
library(ggstance)
library(ggrepel)
library(magrittr)
library(ggh4x)
library(aplot)
library(tibble)
library(tidyverse)
library(ape)
library(ggpubr)
library(cowplot)
library(ggpubr)
library(cowplot)

doMC::registerDoMC(cores = 4)


# basic data
sp_info <- read.csv("/media/Data/zhangz/chip/scripts/info/species.csv", header = TRUE, sep = ",")
species <- sp_info$species
species <- setdiff(species, c("Myotis_ricketti", "Rhinolophus_ferrumequinum", "Hipposideros_larvatus"))

# tree
tree <- read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
tree <- keep.tip(tree, species)

## first, let's draw a plot about the acitivity age and DNA sequence age
## summary both age of each species
act_dna_age <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    anno <- anno[, c("unique_id", "element", "ele_pattern", "species", "tissue", "align_div_time", "align_level", "ec_level", "ec_div_time", "act_sp", "act_time")]
    return(anno)
}, .parallel = TRUE)

nrow(act_dna_age[which(act_dna_age$align_div_time < act_dna_age$act_time), ]) / nrow(act_dna_age)
unique(act_dna_age[which(act_dna_age$align_div_time < act_dna_age$act_time), "act_time"])
head(which(anno$act_time > anno$align_div_time))
nrow(anno[anno$act_time > anno$align_div_time, ]) / nrow(anno)
head(anno[anno$act_time - anno$align_div_time > 10, ])
unique(act_dna_age[which(act_dna_age$act_time - act_dna_age$align_div_time > 10), "species"])

# draw hex heat plot
time_hex_e <- ggplot(act_dna_age[act_dna_age$element == "enhancer", ], aes(x = align_div_time, y = act_time)) +
    geom_hex(bins = 16, color = "grey70", aes(fill = after_stat(count))) +
    scale_fill_gradient2(low = "white", high = "#c01d2e", mid = "#ef4343", midpoint = 200000, labels = label_number(scale = 1 / 1e3)) +
    theme_classic() +
    theme(axis.text = element_text(size = 6)) +
    xlab("DNA sequence age (Ma)") +
    ylab("Activity age (Ma)") +
    labs(fill = "Element counts (1000)") +
    theme(legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = 6), legend.text = element_text(size = 5)) +
    theme(legend.position = "inside", legend.position.inside = c(0.3, 0.8), legend.margin = margin(0, 0, 0, 0)) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2), axis.ticks = element_line(color = "black", size = 0.2))
ggsave("/media/Data/zhangz/chip/analysis/summary2/evol/act_dna_age_hex_e.pdf", time_hex_e, width = 6, height = 5.5, units = "cm")

time_hex_p <- ggplot(act_dna_age[act_dna_age$element == "promoter", ], aes(x = align_div_time, y = act_time)) +
    geom_hex(bins = 16, color = "grey70") +
    scale_fill_gradient2(low = "white", high = "#02263e", mid = "#73b8d5", midpoint = 100000, labels = label_number(scale = 1 / 1e3)) +
    theme_classic() +
    theme(axis.text = element_text(size = 6)) +
    xlab("DNA sequence age (Ma)") +
    ylab("Activity age (Ma)") +
    labs(fill = "Element counts (1000)") +
    theme(legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = 6), legend.text = element_text(size = 5)) +
    theme(legend.position = "inside", legend.position.inside = c(0.3, 0.8), legend.margin = margin(0, 0, 0, 0)) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2), axis.ticks = element_line(color = "black", size = 0.2))
ggsave("/media/Data/zhangz/chip/analysis/summary2/evol/act_dna_age_hex_p.pdf", time_hex_p, width = 6, height = 5.5, units = "cm")

test_tre <- read.tree("/media/Data/zhangz/chip/analysis/summary2/evol/test.newick")
## plot tree
library(ggtree)
mod_tre <- ggtree(test_tre, layout = "rectangular") + geom_tiplab(aes(label = label), size = 2) +
    theme(legend.position = "none")
ggsave("/media/Data/zhangz/chip/analysis/summary2/evol/test_tre.pdf", mod_tre, width = 3, height = 3.5, units = "cm")

rep_enrich <- read.csv("/media/Data/zhangz/chip/analysis/summary2/exap/rep_enrich.csv", header = TRUE, sep = ",")
## draw bubble plot
five_sp <- c("Macaca_mulatta", "Mus_musculus", "Bos_taurus", "Canis_lupus", "Myotis_chinensis")
rec_young_e_rep_p <- ggplot(rep_enrich[rep_enrich$species %in% five_sp, ], aes(x = species, y = rep_class, size = -log10(rec_young_e_p), color = log2(rec_young_e_enrich))) +
    geom_point() +
    theme_light() +
    scale_size_continuous(range = c(0.5, 4), limits = c(0, 300)) +
    theme(text = element_text(size = 8)) +
    scale_color_gradient(low = "white", high = "#c01d2e") +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("young DNA")
rec_anc_e_rep_p <- ggplot(rep_enrich[rep_enrich$species %in% five_sp, ], aes(x = species, y = rep_class, size = -log10(rec_anc_e_p), color = log2(rec_anc_e_enrich))) +
    geom_point() +
    theme_light() +
    scale_size_continuous(range = c(0.5, 4), limits = c(0, 300)) +
    theme(text = element_text(size = 8)) +
    scale_color_gradient(low = "white", high = "#c01d2e") +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_blank()) +
    labs(title = "ancestral DNA", y = "")
rec_e_rep_p <- rec_young_e_rep_p %>% insert_right(rec_anc_e_rep_p, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/rec_e_rep_p.pdf", rec_e_rep_p, height = 12, width = 14, unit = "cm")

rec_young_p_rep_p <- ggplot(rep_enrich[rep_enrich$species %in% five_sp, ], aes(x = species, y = rep_class, size = -log10(rec_young_p_p), color = log2(rec_young_e_enrich))) +
    geom_point() +
    theme_light() +
    scale_size_continuous(range = c(0.5, 2), limits = c(0, 3)) +
    theme(text = element_text(size = 8)) +
    scale_color_gradient(low = "white", high = "#033250") +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("young DNA")
rec_anc_p_rep_p <- ggplot(rep_enrich[rep_enrich$species %in% five_sp, ], aes(x = species, y = rep_class, size = -log10(rec_anc_p_p), color = log2(rec_anc_e_enrich))) +
    geom_point() +
    theme_light() +
    scale_size_continuous(range = c(0.5, 2), limits = c(0, 3)) +
    theme(text = element_text(size = 8)) +
    scale_color_gradient(low = "white", high = "#033250") +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_blank()) +
    labs(title = "ancestral DNA", y = "")
rec_p_rep_p <- rec_young_p_rep_p %>% insert_right(rec_anc_p_rep_p, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/rec_p_rep_p.pdf", rec_p_rep_p, height = 12, width = 14, unit = "cm")

cons_young_e_rep_p <- ggplot(rep_enrich[rep_enrich$species %in% five_sp, ], aes(x = species, y = rep_class, size = -log10(cons_young_e_p), color = log2(cons_young_e_enrich))) +
    geom_point() +
    theme_light() +
    scale_size_continuous(range = c(0.5, 2), limits = c(0, 100)) +
    theme(text = element_text(size = 8)) +
    scale_color_gradient(low = "white", high = "#c01d2e") +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("young DNA")
cons_anc_e_rep_p <- ggplot(rep_enrich[rep_enrich$species %in% five_sp, ], aes(x = species, y = rep_class, size = -log10(cons_anc_e_p), color = log2(cons_anc_e_enrich))) +
    geom_point() +
    theme_light() +
    scale_size_continuous(range = c(0.5, 2), limits = c(0, 100)) +
    theme(text = element_text(size = 8)) +
    scale_color_gradient(low = "white", high = "#c01d2e") +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_blank()) +
    labs(title = "ancestral DNA", y = "")
cons_e_rep_p <- cons_young_e_rep_p %>% insert_right(cons_anc_e_rep_p, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/cons_e_rep_p.pdf", cons_e_rep_p, height = 12, width = 14, unit = "cm")

cons_young_p_rep_p <- ggplot(rep_enrich[rep_enrich$species %in% five_sp, ], aes(x = species, y = rep_class, size = -log10(cons_young_p_p), color = log2(cons_young_p_enrich))) +
    geom_point() +
    theme_light() +
    scale_size_continuous(range = c(0.5, 2), limits = c(0, 2)) +
    theme(text = element_text(size = 8)) +
    scale_color_gradient(low = "white", high = "#033250") +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("young DNA")
cons_anc_p_rep_p <- ggplot(rep_enrich[rep_enrich$species %in% five_sp, ], aes(x = species, y = rep_class, size = -log10(cons_anc_p_p), color = log2(cons_anc_p_enrich))) +
    geom_point() +
    theme_light() +
    scale_size_continuous(range = c(0.5, 2), limits = c(0, 2)) +
    theme(text = element_text(size = 8)) +
    scale_color_gradient(low = "white", high = "#033250") +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_blank()) +
    labs(title = "ancestral DNA", y = "")
cons_p_rep_p <- cons_young_p_rep_p %>% insert_right(cons_anc_p_rep_p, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/cons_p_rep_p.pdf", cons_p_rep_p, height = 12, width = 14, unit = "cm")

pseudo_ts <- adply(species, 1, function(sp) {
    ## the ratio of pseudo ts in enhancer and promoter, which means the ratio of ts in enhancer and promoter in single species
    ## turn out to be evolutionary pleiotropy in other species
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    pseudo_ts_e <- nrow(anno[which(anno$ele_pattern %in% c("B", "K", "L") & (anno$s0_1 != 0 | anno$s1_1 != 0) & anno$element == "enhancer"), ])
    ts_e <- nrow(anno[which(anno$ele_pattern %in% c("B", "K", "L") & anno$element == "enhancer"), ])
    pseudo_ts_p <- nrow(anno[which(anno$ele_pattern %in% c("B", "K", "L") & (anno$s0_1 != 0 | anno$s1_1 != 0) & anno$element == "promoter"), ])
    ts_p <- nrow(anno[which(anno$ele_pattern %in% c("B", "K", "L") & anno$element == "promoter"), ])
    pseudo_ratio_e <- pseudo_ts_e / ts_e
    pseudo_ratio_p <- pseudo_ts_p / ts_p
    res <- data.frame(
        species = sp, pseudo_ts_e = pseudo_ts_e, ts_e = ts_e, pseudo_ratio_e = pseudo_ratio_e,
        pseudo_ts_p = pseudo_ts_p, ts_p = ts_p, pseudo_ratio_p = pseudo_ratio_p
    )
}, .parallel = TRUE)
pseudo_ts <- pseudo_ts[, -1]
summary(pseudo_ts)
write.csv(pseudo_ts, "/media/Data/zhangz/chip/analysis/summary2/exap/pseudo_ts.csv", row.names = FALSE)
pseudo_ts[which(pseudo_ts$pseudo_ratio_e == 0), ]
## boxplot of pseudo ratio in enhancer and promoter
pseudo_ts_melt <- reshape2::melt(pseudo_ts[, c("species", "pseudo_ratio_e", "pseudo_ratio_p")], id.vars = "species", variable.name = "element", value.name = "pseudo_ratio")
pseudo_ts_melt$element <- pseudo_ts_melt$element %>%
    gsub("pseudo_ratio_e", "enhancer", .) %>%
    gsub("pseudo_ratio_p", "promoter", .)
pseudo_ts_box <- ggplot(pseudo_ts_melt, aes(x = element, y = pseudo_ratio, fill = element)) +
    geom_boxplot(color = "black", linewidth = 0.4, outlier.size = 0.2) +
    # geom_line(aes(group = species), color = "grey", linewidth = 0.2) +
    scale_fill_manual(values = c("enhancer" = "#ed4343", "promoter" = "#73b8d5")) +
    scale_x_discrete(labels = c("enhancer" = "Enhancer", "promoter" = "Promoter")) +
    theme_classic() +
    ylab("Ratio") +
    scale_y_continuous(labels = scales::percent) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2), axis.ticks = element_line(color = "black", size = 0.2))
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/pseudo_ts_box.pdf", pseudo_ts_box, height = 6, width = 6, unit = "cm")

false_exap1 <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    exap_e <- nrow(anno[which(anno$element == "enhancer" & (anno$align_div_time - anno$ec_div_time > 10)), ])
    exap_p <- nrow(anno[which(anno$element == "promoter" & (anno$align_div_time - anno$ec_div_time > 10)), ])
    false_exap_e <- nrow(anno[which(anno$element == "enhancer" & (anno$align_div_time - anno$ec_div_time > 10) & (anno$align_div_time - anno$act_time <= 10)), ])
    false_exap_p <- nrow(anno[which(anno$element == "promoter" & (anno$align_div_time - anno$ec_div_time > 10) & (anno$align_div_time - anno$act_time <= 10)), ])
    false_ratio_e <- false_exap_e / exap_e
    false_ratio_p <- false_exap_p / exap_p
    res <- data.frame(
        species = sp, exap_e = exap_e, false_exap_e = false_exap_e, false_ratio_e = false_ratio_e,
        exap_p = exap_p, false_exap_p = false_exap_p, false_ratio_p = false_ratio_p
    )
    return(res)
}, .parallel = TRUE)
false_exap1 <- false_exap1[, -1]
summary(false_exap1)
# wirte.csv(false_exap1, "/media/Data/zhangz/chip/analysis/summary2/exap/false_exap1.csv", row.names = FALSE)
false_exap2 <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    exap_e <- nrow(anno[which(anno$element == "enhancer" & anno$ec_div_time == 0 & anno$align_div_time >= 80), ])
    exap_p <- nrow(anno[which(anno$element == "promoter" & anno$ec_div_time == 0 & anno$align_div_time >= 80), ])
    false_exap_e <- nrow(anno[which(anno$element == "enhancer" & anno$ec_div_time == 0 & anno$align_div_time >= 80 & anno$act_time != 0), ])
    false_exap_p <- nrow(anno[which(anno$element == "promoter" & anno$ec_div_time == 0 & anno$align_div_time >= 80 & anno$act_time != 0), ])
    false_ratio_e <- false_exap_e / exap_e
    false_ratio_p <- false_exap_p / exap_p
    res <- data.frame(
        species = sp, exap_e = exap_e, false_exap_e = false_exap_e, false_ratio_e = false_ratio_e,
        exap_p = exap_p, false_exap_p = false_exap_p, false_ratio_p = false_ratio_p
    )
    return(res)
}, .parallel = TRUE)
false_exap2 <- false_exap2[, -1]
summary(false_exap2)
write.csv(false_exap2, "/media/Data/zhangz/chip/analysis/summary2/exap/false_exap2.csv", row.names = FALSE)
false_exap2_melt <- reshape2::melt(false_exap2[, c("species", "false_ratio_e", "false_ratio_p")], id.vars = "species", variable.name = "element", value.name = "false_ratio")
false_exap2_melt$element <- false_exap2_melt$element %>%
    gsub("false_ratio_e", "Enhancer", .) %>%
    gsub("false_ratio_p", "Promoter", .)
false_exap2_box <- ggplot(false_exap2_melt, aes(x = element, y = false_ratio, fill = element)) +
    geom_boxplot(color = "black", linewidth = 0.4, outlier.size = 0.2) +
    scale_fill_manual(values = c("Enhancer" = "#ed4343", "Promoter" = "#73b8d5")) +
    # scale_x_discrete(labels = c("Enhancer" = "Enhancer", "Promoter" = "Promoter")) +
    theme_classic() +
    ylab("Ratio") +
    scale_y_continuous(labels = scales::percent, breaks = c(0.2, 0.35, 0.5)) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2), axis.ticks = element_line(color = "black", size = 0.2))
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/false_exap2_box.pdf", false_exap2_box, height = 6, width = 6, unit = "cm")

ya_df <- read.csv("/media/Data/zhangz/chip/analysis/summary2/exap/young_anc_dna.csv", header = TRUE, sep = ",")
ya_df <- ya_df[ya_df$species %in% species, ]
y_df <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    y_e <- nrow(anno[which(anno$element == "enhancer" & anno$act_time == 0), ])
    y_p <- nrow(anno[which(anno$element == "promoter" & anno$act_time == 0), ])
    res <- data.frame(
        species = sp, y_e = y_e, y_p = y_p
    )
    return(res)
}, .parallel = TRUE)
ya_df <- merge(ya_df, y_df, by = "species")
ya_df$ye_yd <- ya_df$t3_yd_e_3 / ya_df$y_e
ya_df$yp_yd <- ya_df$t3_yd_p_3 / ya_df$y_p
ya_df$ye_ad <- ya_df$t3_ad_e_3 / ya_df$y_e
ya_df$yp_ad <- ya_df$t3_ad_p_3 / ya_df$y_p
summary(ya_df[which(ya_df$species != "Petaurus_breviceps"), c("ye_yd", "yp_yd", "ye_ad", "yp_ad")])
y_df_t1 <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    y_e_t1 <- nrow(anno[which(anno$element == "enhancer" & anno$ec_div_time == 0), ])
    y_p_t1 <- nrow(anno[which(anno$element == "promoter" & anno$ec_div_time == 0), ])
    res <- data.frame(
        species = sp, y_e_t1 = y_e_t1, y_p_t1 = y_p_t1
    )
    return(res)
}, .parallel = TRUE)
y_df_t1 <- y_df_t1[, -1]

ya_df <- merge(ya_df, y_df_t1, by = "species")
ya_df$t1_ye_yd <- ya_df$t1_yd_e_3 / ya_df$y_e_t1
ya_df$t1_yp_yd <- ya_df$t1_yd_p_3 / ya_df$y_p_t1
ya_df$t1_ye_ad <- ya_df$t1_ad_e_3 / ya_df$y_e_t1
ya_df$t1_yp_ad <- ya_df$t1_ad_p_3 / ya_df$y_p_t1
summary(ya_df[which(ya_df$species != "Petaurus_breviceps"), c("t1_ye_yd", "t1_yp_yd", "t1_ye_ad", "t1_yp_ad")])
ya_melt <- reshape2::melt(ya_df, id.vars = "species", variable.name = "t_ay_ep_m", value.name = "cre_num")
ya_melt <- separate(ya_melt, t_ay_ep_m, into = c("tissue_num", "DNA_age", "element", "method"), sep = "_")
ya_melt <- ya_melt[ya_melt$species != "Petaurus_breviceps", ]
compars <- list(c("yd", "ad"))
ya_box_p_1t <- ggplot(
    ya_melt[which(ya_melt$element == "p" & ya_melt$method == "3" & ya_melt$tissue_num == "t1"), ],
    aes(x = DNA_age, y = cre_num / 1000, color = DNA_age)
) +
    geom_boxplot(fill = "white", outlier.size = 0.2, linewidth = 0.2) +
    geom_line(aes(group = species), color = "grey", linewidth = 0.1) +
    stat_compare_means(aes(x = DNA_age), comparisons = compars, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_color_manual(values = c(ad = "#02263e", yd = "#73b8d5")) +
    theme_classic() +
    scale_x_discrete(limits = c("yd", "ad"), labels = c(yd = "Young DNA", ad = "Ancestral DNA")) +
    ylab("Promoter number (*1000)") +
    ggtitle("Recently evolved promoter") +
    theme(text = element_text(size = 6), axis.title.x = element_blank()) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2), axis.ticks = element_line(color = "black", linewidth = 0.2))
ya_box_p_3t <- ggplot(
    ya_melt[which(ya_melt$element == "p" & ya_melt$method == "3" & ya_melt$tissue_num == "t3"), ],
    aes(x = DNA_age, y = cre_num / 1000, color = DNA_age)
) +
    geom_boxplot(fill = "white", outlier.size = 0.2, linewidth = 0.3) +
    geom_line(aes(group = species), color = "grey", linewidth = 0.1) +
    stat_compare_means(aes(x = DNA_age), comparisons = compars, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_color_manual(values = c(ad = "#02263e", yd = "#73b8d5")) +
    theme_classic() +
    scale_x_discrete(limits = c("yd", "ad"), labels = c(yd = "Young DNA", ad = "Ancestral DNA")) +
    ylab("Promoter number (*1000)") +
    ggtitle("Recently evolved promoter") +
    theme(text = element_text(size = 6), axis.title.x = element_blank()) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2), axis.ticks = element_line(color = "black", linewidth = 0.2))
ya_box_e_1t <- ggplot(
    ya_melt[which(ya_melt$element == "e" & ya_melt$method == "3" & ya_melt$tissue_num == "t1"), ],
    aes(x = DNA_age, y = cre_num / 1000, color = DNA_age)
) +
    geom_boxplot(fill = "white", outlier.size = 0.2, linewidth = 0.2) +
    geom_line(aes(group = species), color = "grey", linewidth = 0.1) +
    stat_compare_means(aes(x = DNA_age), comparisons = compars, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_color_manual(values = c(ad = "#c01d2e", yd = "#ef4343")) +
    theme_classic() +
    scale_x_discrete(limits = c("yd", "ad"), labels = c(yd = "Young DNA", ad = "Ancestral DNA")) +
    ylab("Enhancer number (*1000)") +
    ggtitle("Recently evolved enhancer") +
    theme(text = element_text(size = 6), axis.title.x = element_blank()) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2), axis.ticks = element_line(color = "black", linewidth = 0.2))
ya_box_e_3t <- ggplot(
    ya_melt[which(ya_melt$element == "e" & ya_melt$method == "3" & ya_melt$tissue_num == "t3"), ],
    aes(x = DNA_age, y = cre_num / 1000, color = DNA_age)
) +
    geom_boxplot(fill = "white", outlier.size = 0.2, linewidth = 0.3) +
    geom_line(aes(group = species), color = "grey", linewidth = 0.1) +
    stat_compare_means(aes(x = DNA_age), comparisons = compars, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_color_manual(values = c(ad = "#c01d2e", yd = "#ef4343")) +
    theme_classic() +
    scale_x_discrete(limits = c("yd", "ad"), labels = c(yd = "Young DNA", ad = "Ancestral DNA")) +
    ylab("Enhancer number (*1000)") +
    ggtitle("Recently evolved enhancer") +
    theme(text = element_text(size = 6), axis.title.x = element_blank()) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2), axis.ticks = element_line(color = "black", linewidth = 0.2))

ya_box <- plot_grid(ya_box_p_1t, ya_box_e_1t, ya_box_p_3t, ya_box_e_3t, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/young_anc_dna.pdf", ya_box, width = 20, height = 20, units = "cm")

recent_cre <- plot_grid(ya_box_e_3t, ya_box_p_3t, nrow = 1, ncol = 2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/recent_cre.pdf", recent_cre, width = 10, height = 5, units = "cm")

recent_cre_1t <- plot_grid(ya_box_e_1t, ya_box_p_1t, nrow = 1, ncol = 2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/recent_cre_1t.pdf", recent_cre_1t, width = 10, height = 5, units = "cm")

rep_ratio <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    anno <- anno[, c(colnames(anno)[1:44], "s0_0", "s0_1", "s1_0", "s1_1", "act_time", "act_sp")]
    # rep_class <- paste(anno$repeat_class, collapse = ";") %>%
    #     strsplit(";") %>%
    #     unlist() %>%
    #     unique()
    # rep_class <- rep_class[!rep_class %in% c("", "Low_complexity", "Simple_repeat", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", "srpRNA", "RNA", "DNA", "LTR", "LINE", "SINE", "Unknown")]
    rep_ele <- mdply(expand.grid(tis = unique(anno$tissue), ele = unique(anno$element)), function(tis, ele) {
        ele_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele), ])
        rep_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$repeat_class != ""), ])
        rep_ratio <- rep_num / ele_num
        yy_ele_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time == 0 & anno$align_div_time == 0), ])
        yy_rep_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time == 0 & anno$align_div_time == 0 & anno$repeat_class != ""), ])
        yy_rep_ratio <- yy_rep_num / yy_ele_num
        aa_ele_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time >= 80 & anno$act_time <= 100 & anno$align_div_time >= 80 & anno$align_div_time <= 100), ])
        aa_rep_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time >= 80 & anno$act_time <= 100 & anno$align_div_time >= 80 & anno$align_div_time <= 100 & anno$repeat_class != ""), ])
        aa_rep_ratio <- aa_rep_num / aa_ele_num
        res <- data.frame(
            tissue = tis, element = ele, ele_num = ele_num, rep_num = rep_num, rep_ratio = rep_ratio,
            yy_ele_num = yy_ele_num, yy_rep_num = yy_rep_num, yy_rep_ratio = yy_rep_ratio,
            aa_ele_num = aa_ele_num, aa_rep_num = aa_rep_num, aa_rep_ratio = aa_rep_ratio
        )
        return(res)
    }, .parallel = TRUE)
    rep_ele$species <- sp
    return(rep_ele)
}, .parallel = TRUE)
rep_ratio <- rep_ratio[, -1]
write.csv(rep_ratio, "/media/Data/zhangz/chip/analysis/summary2/exap/rep_ratio.csv", row.names = FALSE)
rep_ratio <- read.csv("/media/Data/zhangz/chip/analysis/summary2/exap/rep_ratio.csv", header = TRUE, sep = ",")

rep_ratio <- rep_ratio[which(rep_ratio$species != "Petaurus_breviceps"), ]
## draw two plots, both are paired boxplot
e_yy_ratio_df <- rep_ratio[which(rep_ratio$element == "enhancer"), c("species", "tissue", "element", "yy_rep_ratio", "rep_ratio")] %>%
    gather(key = "ratio_type", value = "ratio", yy_rep_ratio, rep_ratio)
e_yy_ratio_df$ratio_type <- factor(e_yy_ratio_df$ratio_type, levels = c("yy_rep_ratio", "rep_ratio"))
comps <- list(c("yy_rep_ratio", "rep_ratio"))
e_yy_ratio_p <- ggplot(e_yy_ratio_df, aes(x = tissue, y = ratio, fill = ratio_type)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.2) +
    stat_compare_means(comparisons = list(c("yy_rep_ratio", "rep_ratio")), method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_fill_manual(values = c("yy_rep_ratio" = "#ef4343", "rep_ratio" = "grey60")) +
    # stat_compare_means(comparisons = comps, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    geom_signif(comparisons = comps, map_signif_level = TRUE, textsize = 3) +
    theme_classic() +
    ylab("Ratio") +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", linewidth = 0.2))

p_yy_ratio_df <- rep_ratio[which(rep_ratio$element == "promoter"), c("species", "tissue", "element", "yy_rep_ratio", "rep_ratio")] %>%
    gather(key = "ratio_type", value = "ratio", yy_rep_ratio, rep_ratio)
p_yy_ratio_df$ratio_type <- factor(p_yy_ratio_df$ratio_type, levels = c("yy_rep_ratio", "rep_ratio"))

p_yy_ratio_p <- ggplot(p_yy_ratio_df, aes(x = tissue, y = ratio, fill = ratio_type)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.2) +
    stat_compare_means(comparisons = list(c("yy_rep_ratio", "rep_ratio")), method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_fill_manual(values = c("yy_rep_ratio" = "#73b8d5", "rep_ratio" = "grey60")) +
    theme_classic() +
    ylab("Ratio") +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", linewidth = 0.2))

e_aa_ratio_df <- rep_ratio[which(rep_ratio$element == "enhancer"), c("species", "tissue", "element", "aa_rep_ratio", "rep_ratio")] %>%
    gather(key = "ratio_type", value = "ratio", aa_rep_ratio, rep_ratio)
e_aa_ratio_df$ratio_type <- factor(e_aa_ratio_df$ratio_type, levels = c("aa_rep_ratio", "rep_ratio"))
e_aa_ratio_p <- ggplot(e_aa_ratio_df, aes(x = tissue, y = ratio, fill = ratio_type)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.2) +
    stat_compare_means(aes(group = ratio_type)) +
    scale_fill_manual(values = c("aa_rep_ratio" = "#c01d2e", "rep_ratio" = "grey60")) +
    theme_classic() +
    ylab("Ratio") +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", linewidth = 0.2))

p_aa_ratio_df <- rep_ratio[which(rep_ratio$element == "promoter"), c("species", "tissue", "element", "aa_rep_ratio", "rep_ratio")] %>%
    gather(key = "ratio_type", value = "ratio", aa_rep_ratio, rep_ratio)


e_ya_ratio_df <- rep_ratio[which(rep_ratio$element == "enhancer"), c("species", "tissue", "element", "yy_rep_ratio", "aa_rep_ratio", "rep_ratio")] %>%
    gather(key = "ratio_type", value = "ratio", yy_rep_ratio, aa_rep_ratio, rep_ratio)
e_ya_ratio_df$ratio_type <- factor(e_ya_ratio_df$ratio_type, levels = c("yy_rep_ratio", "aa_rep_ratio", "rep_ratio"))
comps <- list(c("yy_rep_ratio", "rep_ratio"), c("aa_rep_ratio", "rep_ratio"))
e_ya_ratio_p <- ggplot(e_ya_ratio_df, aes(x = tissue, y = ratio, fill = ratio_type)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.05, color = "black", width = 0.5) +
    # stat_compare_means(aes(group = ratio_type), comparisons = comps, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    # stat_compare_means(aes(group = ratio_type), comparisons = list(c("yy_rep_ratio", "rep_ratio")), method = "wilcox.test", label = "p.signif", paired = TRUE) +
    stat_compare_means(aes(group = ratio_type), method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_fill_manual(
        values = c("yy_rep_ratio" = "#ef4343", "aa_rep_ratio" = "#c01d2e", "rep_ratio" = "grey80"),
        labels = c("yy_rep_ratio" = "Young enhancer", "aa_rep_ratio" = "Ancestral enhancer", "rep_ratio" = "All enhancer")
    ) +
    theme_classic() +
    ylab("Repeat element proportion") +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "inside", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", linewidth = 0.2)) +
    theme(legend.position.inside = c(0.5, 0.1)) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.direction = "horizontal") +
    theme(legend.title = element_blank())
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/e_ya_ratio_p.pdf", e_ya_ratio_p, height = 6, width = 8, unit = "cm")

p_ya_ratio_df <- rep_ratio[which(rep_ratio$element == "promoter"), c("species", "tissue", "element", "yy_rep_ratio", "aa_rep_ratio", "rep_ratio")] %>%
    gather(key = "ratio_type", value = "ratio", yy_rep_ratio, aa_rep_ratio, rep_ratio)
p_ya_ratio_df$ratio_type <- factor(p_ya_ratio_df$ratio_type, levels = c("yy_rep_ratio", "aa_rep_ratio", "rep_ratio"))
p_ya_ratio_p <- ggplot(p_ya_ratio_df, aes(x = tissue, y = ratio, fill = ratio_type)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.1, color = "black", width = 0.5) +
    stat_compare_means(aes(group = ratio_type), method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_fill_manual(
        values = c("yy_rep_ratio" = "#73b8d5", "aa_rep_ratio" = "#02263e", "rep_ratio" = "grey80"),
        labels = c("yy_rep_ratio" = "Young enhancer", "aa_rep_ratio" = "Ancestral enhancer", "rep_ratio" = "All enhancer")
    ) +
    theme_classic() +
    ylab("Repeat element proportion") +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "inside", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", linewidth = 0.2)) +
    theme(legend.position.inside = c(0.5, 0.1)) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.direction = "horizontal") +
    theme(legend.title = element_blank())
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/p_ya_ratio_p.pdf", p_ya_ratio_p, height = 6, width = 8, unit = "cm")

rep_enrich2 <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    anno <- anno[, c(colnames(anno)[1:44], "s0_0", "s0_1", "s1_0", "s1_1", "act_time", "act_sp")]
    rep_class <- paste(anno$repeat_class, collapse = ";") %>%
        strsplit(";") %>%
        unlist() %>%
        unique()
    rep_class <- rep_class[!rep_class %in% c("", "Low_complexity", "Simple_repeat", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", "srpRNA", "RNA", "DNA", "LTR", "LINE", "SINE", "Unknown")]
    young_e_row <- which(anno$act_time == 0 & anno$align_div_time == 0 & anno$element == "enhancer")
    young_p_row <- which(anno$act_time == 0 & anno$align_div_time == 0 & anno$element == "promoter")
    anc_e_row <- which(anno$act_time >= 80 & anno$act_time <= 100 & anno$align_div_time >= 80 & anno$align_div_time <= 100 & anno$element == "enhancer")
    anc_p_row <- which(anno$act_time >= 80 & anno$act_time <= 100 & anno$align_div_time >= 80 & anno$align_div_time <= 100 & anno$element == "promoter")

    e_row <- which(anno$element == "enhancer")
    p_row <- which(anno$element == "promoter")
    young_e_num <- sum(anno$act_time == 0 & anno$align_div_time == 0 & anno$element == "enhancer")
    young_p_num <- sum(anno$act_time == 0 & anno$align_div_time == 0 & anno$element == "promoter")
    anc_e_num <- sum(anno$act_time >= 80 & anno$act_time <= 100 & anno$align_div_time >= 80 & anno$align_div_time <= 100 & anno$element == "enhancer")
    anc_p_num <- sum(anno$act_time >= 80 & anno$act_time <= 100 & anno$align_div_time >= 80 & anno$align_div_time <= 100 & anno$element == "promoter")

    e_num <- sum(anno$element == "enhancer")
    p_num <- sum(anno$element == "promoter")
    rep_num <- adply(rep_class, 1, function(rep) {
        res <- data.frame(
            rep_class = rep,
            young_e = sum(grepl(rep, anno$repeat_class[young_e_row])),
            young_p = sum(grepl(rep, anno$repeat_class[young_p_row])),
            anc_e = sum(grepl(rep, anno$repeat_class[anc_e_row])),
            anc_p = sum(grepl(rep, anno$repeat_class[anc_p_row])),
            e_rep = sum(grepl(rep, anno$repeat_class[e_row])),
            p_rep = sum(grepl(rep, anno$repeat_class[p_row])),
            young_e_num = young_e_num,
            young_p_num = young_p_num,
            anc_e_num = anc_e_num,
            anc_p_num = anc_p_num,
            e_num = e_num,
            p_num = p_num
        )
        res$young_e_p <- phyper(res$young_e, res$e_rep, e_num - res$e_rep, young_e_num, lower.tail = FALSE)
        res$young_e_enrich <- (res$young_e / res$young_e_num) / (res$e_rep / res$e_num)
        res$young_p_p <- phyper(res$young_p, res$p_rep, p_num - res$p_rep, young_p_num, lower.tail = FALSE)
        res$young_p_enrich <- (res$young_p / res$young_p_num) / (res$p_rep / res$p_num)
        res$anc_e_p <- phyper(res$anc_e, res$e_rep, e_num - res$e_rep, anc_e_num, lower.tail = FALSE)
        res$anc_e_enrich <- (res$anc_e / res$anc_e_num) / (res$e_rep / res$e_num)
        res$anc_p_p <- phyper(res$anc_p, res$p_rep, p_num - res$p_rep, anc_p_num, lower.tail = FALSE)
        res$anc_p_enrich <- (res$anc_p / res$anc_p_num) / (res$p_rep / res$p_num)
        return(res)
    }, .parallel = TRUE)
    rep_num <- rep_num[, -1]
    rep_num$species <- sp
    return(rep_num)
}, .parallel = TRUE)
summary(rep_enrich2)
rep_enrich2 <- rep_enrich2[, -1]
rep_enrich2 <- rep_enrich2[, c(1, 22, 2:21)]
rep_enrich2$order <- sp_info$order[match(rep_enrich2$species, sp_info$species)]
write.csv(rep_enrich2, "/media/Data/zhangz/chip/analysis/summary2/exap/rep_enrich2.csv", row.names = FALSE)

rep_enrich2 <- read.csv("/media/Data/zhangz/chip/analysis/summary2/exap/rep_enrich2.csv", header = TRUE, sep = ",")
rep_enrich2 <- rep_enrich2[which(rep_enrich2$species != "Petaurus_breviceps"), ]
rep_enrich_sum <- rep_enrich2[, c("rep_class", "species", "order", "young_e_p", "young_e_enrich", "young_p_p", "young_p_enrich", "anc_e_p", "anc_e_enrich", "anc_p_p", "anc_p_enrich")] %>%
    gather(key = "enrich_type", value = "enrich", young_e_enrich, young_p_enrich, anc_e_enrich, anc_p_enrich) %>%
    separate(enrich_type, into = c("time", "element"), sep = "_") %>%
    mutate(element = ifelse(element == "p", "promoter", "enhancer")) %>%
    mutate(time = ifelse(time == "anc", "ancestral", "young")) %>%
    gather(key = "p_value_type", value = "p_value", young_e_p, young_p_p, anc_e_p, anc_p_p) %>%
    separate(p_value_type, into = c("time", "element"), sep = "_") %>%
    mutate(element = ifelse(element == "p", "promoter", "enhancer")) %>%
    mutate(time = ifelse(time == "anc", "ancestral", "young"))
rep_enrich_mid <- rep_enrich_sum %>%
    aggregate(enrich ~ rep_class + element + time + order, data = ., FUN = median)
rep_p_mid <- rep_enrich_sum %>%
    aggregate(p_value ~ rep_class + element + time + order, data = ., FUN = median)
rep_enrich_mid <- merge(rep_enrich_mid, rep_p_mid, by = c("rep_class", "element", "time", "order"))

rep_enrich_mid <- rep_enrich_mid[!rep_enrich_mid$rep_class %in% c("ARTEFACT", "DNA/Merlin", "DNA/TcMar-Pogo", "LINE/L1−Tx1", "LTR/DIRS", "Satellite/subtelo", "DNA/hAT-hAT19", "LINE/L1-Tx1", "Satellite/centr"), ]
e_rep_plot <- ggplot(rep_enrich_mid[which(rep_enrich_mid$element == "enhancer"), ], aes(x = order, y = rep_class, color = log2(enrich), size = -log10(p_value))) +
    geom_point() +
    scale_color_gradient(low = "white", high = "#c01d2e") +
    scale_size_continuous(range = c(1, 5)) +
    theme_bw() +
    theme(legend.position = "none") +
    ggtitle("Repeat class enrichment in enhancer") +
    theme(text = element_text(size = 8), plot.title = element_text(hjust = 0.5)) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", size = 0.2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.title = element_blank()) +
    facet_grid(. ~ time)
p_rep_plot <- ggplot(rep_enrich_mid[which(rep_enrich_mid$element == "promoter"), ], aes(x = order, y = rep_class, color = log2(enrich), size = -log10(p_value))) +
    geom_point() +
    scale_color_gradient(low = "white", high = "#02263e") +
    scale_size_continuous(range = c(1, 5)) +
    theme_bw() +
    theme(legend.position = "none") +
    ggtitle("Repeat class enrichment in promoter") +
    theme(text = element_text(size = 8), plot.title = element_text(hjust = 0.5)) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", size = 0.2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.title = element_blank()) +
    facet_grid(. ~ time)
rep_plot <- e_rep_plot %>% insert_right(p_rep_plot, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/rep_enrich2.pdf", rep_plot, width = 18, height = 18, units = "cm")

sele_rep <- c(
    "SINE/tRNA-RTE", "SINE/tRNA-Deu", "SINE/tRNA", "SINE/MIR", "SINE/ID", "SINE/B4", "SINE/B2", "SINE/Alu", "SINE/5S-Deu-L2",
    "LINE/RTE-X", "LINE/RTE-BovB", "LINE/Penelope", "LINE/L2", "LINE/L1", "LINE/I-Jockey", "LINE/Dong-R4", "LINE/CR1",
    "LTR/Gypsy", "LTR/ERVL-MaLR", "LTR/ERVL", "LTR/ERVK", "LTR/ERV1",
    "DNA/TcMar-Tigger", "DNA/TcMar-Mariner", "DNA/TcMar", "DNA/hAT-Tip100", "DNA/hAT-Charlie", "DNA/hAT-Blackjack", "DNA/hAT"
)
rep_enrich_mid$order <- factor(rep_enrich_mid$order, levels = c("Artiodactyla", "Carnivora", "Perissodactyla", "Chiroptera", "Eulipotyphla", "Lagomorpha", "Rodentia", "Primates", "Scandentia", "Hyracoidea"))
rep_enrich_mid_tmp <- rep_enrich_mid[which(rep_enrich_mid_tmp$rep_class %in% sele_rep), ]
rep_enrich_mid_tmp$rep_class <- factor(rep_enrich_mid_tmp$rep_class, levels = sele_rep)
rep_enrich_mid_tmp$p_value[rep_enrich_mid_tmp$p_value < 1e-200] <- 1e-200
e_rep_plot_sele <- ggplot(rep_enrich_mid_tmp[which(rep_enrich_mid_tmp$element == "enhancer"), ], aes(x = order, y = rep_class, color = log2(enrich), size = -log10(p_value))) +
    geom_point() +
    scale_color_gradient(low = "white", high = "#c01d2e") +
    scale_size_continuous(range = c(0.1, 5), limits = c(1.30103, 200)) + # , breaks = c(1.30103, 20, 40, 314)
    theme_bw() +
    theme(legend.position = "right") +
    ggtitle("Repeat class enrichment in enhancer") +
    theme(text = element_text(size = 8), plot.title = element_text(hjust = 0.5)) +
    theme(axis.line = element_blank()) +
    theme(axis.ticks = element_line(color = "black", size = 0.2)) +
    # theme(axis.ticks.y.right = element_line(color = "black", size = 0.2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.title = element_blank()) +
    theme(panel.grid = element_line(linetype = "dotted")) +
    facet_grid(. ~ time)
p_rep_plot_sele <- ggplot(rep_enrich_mid_tmp[which(rep_enrich_mid_tmp$element == "promoter"), ], aes(x = order, y = rep_class, color = log2(enrich), size = -log10(p_value))) +
    geom_point() +
    scale_color_gradient(low = "white", high = "#02263e") +
    scale_size_continuous(range = c(0.1, 5), limits = c(1.30103, 200)) +
    theme_bw() +
    theme(legend.position = "right") +
    ggtitle("Repeat class enrichment in promoter") +
    theme(text = element_text(size = 8), plot.title = element_text(hjust = 0.5)) +
    theme(axis.line = element_blank()) +
    theme(axis.ticks = element_line(color = "black", size = 0.2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.title = element_blank()) +
    theme(panel.grid = element_line(linetype = "dotted")) +
    facet_grid(. ~ time)
rep_plot_sele <- e_rep_plot_sele %>% insert_right(p_rep_plot_sele, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/rep_enrich2_sele.pdf", rep_plot_sele, width = 24, height = 15, units = "cm")

rep_ratio <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    anno <- anno[, c(colnames(anno)[1:44], "s0_0", "s0_1", "s1_0", "s1_1", "act_time", "act_sp")]
    # rep_class <- paste(anno$repeat_class, collapse = ";") %>%
    #     strsplit(";") %>%
    #     unlist() %>%
    #     unique()
    # rep_class <- rep_class[!rep_class %in% c("", "Low_complexity", "Simple_repeat", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", "srpRNA", "RNA", "DNA", "LTR", "LINE", "SINE", "Unknown")]
    rep_ele <- mdply(expand.grid(tis = unique(anno$tissue), ele = unique(anno$element)), function(tis, ele) {
        ele_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele), ])
        rep_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$repeat_class != ""), ])
        rep_ratio <- rep_num / ele_num
        ysya_ele_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time == 0 & anno$align_div_time == 0), ])
        ysya_rep_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time == 0 & anno$align_div_time == 0 & anno$repeat_class != ""), ])
        ysya_rep_ratio <- ysya_rep_num / ysya_ele_num
        asaa_ele_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time >= 80 & anno$act_time <= 100 & anno$align_div_time >= 80 & anno$align_div_time <= 100), ])
        asaa_rep_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time >= 80 & anno$act_time <= 100 & anno$align_div_time >= 80 & anno$align_div_time <= 100 & anno$repeat_class != ""), ])
        asaa_rep_ratio <- asaa_rep_num / asaa_ele_num
        asya_ele_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time == 0 & anno$align_div_time >= 80 & anno$align_div_time <= 100), ])
        asya_rep_num <- nrow(anno[which(anno$tissue == tis & anno$element == ele & anno$act_time == 0 & anno$align_div_time >= 80 & anno$align_div_time <= 100 & anno$repeat_class != ""), ])
        asya_rep_ratio <- asya_rep_num / asya_ele_num
        res <- data.frame(
            tissue = tis, element = ele, ele_num = ele_num, rep_num = rep_num, rep_ratio = rep_ratio,
            yy_ele_num = ysya_ele_num, yy_rep_num = ysya_rep_num, yy_rep_ratio = ysya_rep_ratio,
            aa_ele_num = asaa_ele_num, aa_rep_num = asaa_rep_num, aa_rep_ratio = asaa_rep_ratio,
            asya_ele_num = asya_ele_num, asya_rep_num = asya_rep_num, asya_rep_ratio = asya_rep_ratio
        )
        return(res)
    }, .parallel = TRUE)
    rep_ele$species <- sp
    return(rep_ele)
}, .parallel = TRUE)
rep_ratio <- rep_ratio[, -1]
write.csv(rep_ratio, "/media/Data/zhangz/chip/analysis/summary2/exap/rep_ratio.csv", row.names = FALSE)
rep_ratio <- read.csv("/media/Data/zhangz/chip/analysis/summary2/exap/rep_ratio.csv", header = TRUE, sep = ",")


rep_ratio <- rep_ratio[which(rep_ratio$species != "Petaurus_breviceps"), ]
summary(rep_ratio[rep_ratio$element == "enhancer", ])
summary(rep_ratio[rep_ratio$element == "promoter", ])
e_ya_ratio_df <- rep_ratio[which(rep_ratio$element == "enhancer"), c("species", "tissue", "element", "yy_rep_ratio", "aa_rep_ratio", "asya_rep_ratio", "rep_ratio")] %>%
    gather(key = "ratio_type", value = "ratio", yy_rep_ratio, aa_rep_ratio, asya_rep_ratio, rep_ratio)
e_ya_ratio_df$ratio_type <- factor(e_ya_ratio_df$ratio_type, levels = c("yy_rep_ratio", "aa_rep_ratio", "asya_rep_ratio", "rep_ratio"))
comps <- list(c("yy_rep_ratio", "rep_ratio"), c("aa_rep_ratio", "rep_ratio"), c("asya_rep_ratio", "rep_ratio"))
e_ya_ratio_p <- ggplot(e_ya_ratio_df, aes(x = tissue, y = ratio, fill = ratio_type)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.01, color = "black", width = 0.5) +
    # stat_compare_means(aes(group = ratio_type), comparisons = comps, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    # stat_compare_means(aes(group = ratio_type), comparisons = list(c("yy_rep_ratio", "rep_ratio")), method = "wilcox.test", label = "p.signif", paired = TRUE) +
    stat_compare_means(aes(group = ratio_type), method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_fill_manual(
        values = c("yy_rep_ratio" = "#ef4343", "aa_rep_ratio" = "#c01d2e", "asya_rep_ratio" = "#e83747", "rep_ratio" = "grey80"),
        labels = c("yy_rep_ratio" = "Young enhancer", "aa_rep_ratio" = "Ancestral enhancer", "asya_rep_ratio" = "Young enhancer in ancestral sequences", "rep_ratio" = "All enhancer")
    ) +
    theme_classic() +
    ylab("Repeat element proportion") +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "inside", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", linewidth = 0.2)) +
    theme(legend.position.inside = c(0.5, 0.1)) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.direction = "horizontal") +
    theme(legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/e_ya_ratio_p_2.pdf", e_ya_ratio_p, height = 6, width = 8, unit = "cm")

for (tis in c("Brain", "Kidney", "Liver")) {
    print(tis)
    print(
        wilcox.test(e_ya_ratio_df[which(e_ya_ratio_df$tissue == tis & e_ya_ratio_df$ratio_type == "asya_rep_ratio"), ]$ratio,
            e_ya_ratio_df[which(e_ya_ratio_df$tissue == tis & e_ya_ratio_df$ratio_type == "rep_ratio"), ]$ratio,
            paired = TRUE
        )$p.value
    )
}
for (tis in c("Brain", "Kidney", "Liver")) {
    print(tis)
    print(
        wilcox.test(p_ya_ratio_df[which(p_ya_ratio_df$tissue == tis & e_ya_ratio_df$ratio_type == "asya_rep_ratio"), ]$ratio,
            p_ya_ratio_df[which(p_ya_ratio_df$tissue == tis & p_ya_ratio_df$ratio_type == "rep_ratio"), ]$ratio,
            paired = TRUE
        )$p.value
    )
}

p_ya_ratio_df <- rep_ratio[which(rep_ratio$element == "promoter"), c("species", "tissue", "element", "yy_rep_ratio", "aa_rep_ratio", "asya_rep_ratio", "rep_ratio")] %>%
    gather(key = "ratio_type", value = "ratio", yy_rep_ratio, aa_rep_ratio, asya_rep_ratio, rep_ratio)
p_ya_ratio_df$ratio_type <- factor(p_ya_ratio_df$ratio_type, levels = c("yy_rep_ratio", "aa_rep_ratio", "asya_rep_ratio", "rep_ratio"))
p_ya_ratio_p <- ggplot(p_ya_ratio_df, aes(x = tissue, y = ratio, fill = ratio_type)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.01, color = "black", width = 0.5) +
    stat_compare_means(aes(group = ratio_type), method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_fill_manual(
        values = c("yy_rep_ratio" = "#73b8d5", "aa_rep_ratio" = "#02263e", "asya_rep_ratio" = "#457b9d", "rep_ratio" = "grey80"),
        labels = c("yy_rep_ratio" = "Young enhancer", "aa_rep_ratio" = "Ancestral enhancer", "asya_rep_ratio" = "Young promoter in ancestral sequences", "rep_ratio" = "All enhancer")
    ) +
    theme_classic() +
    ylab("Repeat element proportion") +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "inside", axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(color = "black", linewidth = 0.2)) +
    theme(legend.position.inside = c(0.5, 0.1)) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.direction = "horizontal") +
    theme(legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/p_ya_ratio_p_2.pdf", p_ya_ratio_p, height = 6, width = 8, unit = "cm")
