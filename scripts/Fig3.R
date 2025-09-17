# draw plot about gene annotation, repeat elements annotation and entropies
# load libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(tidyr)
library(aplot)
library(ggstance)
library(dplyr)
library(ggpubr)
library(ggtree)
library(ggridges)
library(ggunchained)
library(Rmisc)
library(nlme)
library(rr2)
library(ggimage)

# read data
dis_gene <- read.csv("/media/Data/zhangz/chip/analysis/summary2/each_sp/dis_gene_num.csv")
dis_gene <- dis_gene[!is.na(dis_gene$median), ]
dis_gene <- dis_gene[, -1]

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

# turn dis_gene into wide dataframe
# we need columns including species, ele, abc, median_dis, median_gene_number

dis_gene <- dis_gene[dis_gene$species %in% species, ]
vars <- data.frame(old = unique(dis_gene$variable))
vars$new <- c("all_distance", "ABC_distance", "other_distance", "all_genen", "ABC_genen", "other_genen")
for (i in 1:nrow(vars)) {
    dis_gene$variable <- gsub(paste0("^", vars$old[i]), vars$new[i], dis_gene$variable)
}

# separate variable into two columns
dis_gene <- dis_gene %>%
    separate(., variable, into = c("annotation", "var"))
dis_gene_wide <- dis_gene %>%
    spread(key = "var", value = "median")
dis_gene_wide <- dis_gene_wide[which(dis_gene_wide$annotation != "all" & dis_gene_wide$ele != "all"), ]
dis_gene_point <- ggplot(dis_gene_wide, aes(x = genen, y = log10(distance), color = ele, shape = annotation)) +
    geom_point(size = 1.7) +
    theme_classic() +
    scale_color_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    scale_y_continuous(labels = c(1, 1000, 1000000), breaks = c(0, 3, 6)) +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    labs(x = "Gene number", y = "Distance (bp)", color = "Element", shape = "Annotation method") +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.4), text = element_text(size = 8)) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1))
dis_box <- ggplot(dis_gene_wide, aes(x = log10(distance), y = annotation, fill = ele, group = interaction(annotation, ele))) +
    geom_density_ridges(linewidth = 0.1, color = "black", alpha = 0.7) +
    theme_classic() +
    scale_fill_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    # scale_color_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    # scale_y_continuous(labels = c(1, 1000, 1000000), breaks = c(0, 3, 6)) +
    labs(y = "Annotation") +
    theme(legend.position = "none", text = element_text(size = 8), axis.title.y = element_blank()) +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_flip()
genen_box <- ggplot(dis_gene_wide, aes(y = annotation, x = genen, fill = ele, group = interaction(annotation, ele))) +
    geom_density_ridges(linewidth = 0.1, color = "black", alpha = 0.7) +
    theme_classic() +
    scale_fill_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    # scale_color_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    labs(y = "Annotation") +
    theme(legend.position = "none", text = element_text(size = 8), axis.title.x = element_blank()) +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

dis_gene_plot <- dis_gene_point %>%
    insert_right(dis_box, width = 0.2) %>%
    insert_top(genen_box, height = 0.2)
# dis_gene_plot <- dis_gene_plot
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/dis_gene_plot.pdf", dis_gene_plot, width = 12, height = 8, dpi = 300, units = "cm")

res <- read.csv("/media/Data/zhangz/chip/analysis/summary2/each_sp/cre_per_gene.csv")
res_m <- reshape2::melt(res, id.vars = c("gene", "species"), measure.vars = c("all_enhancer", "all_promoter"), value.name = "count")
res_mid <- tapply(res_m$count, list(res_m$variable, res_m$species), median) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var1") %>%
    gather(key = "Var2", value = "value", -Var1)
anno_count_ratio <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/anno_count_ratio.csv", header = TRUE)
anno_count_ratio$type <- factor(anno_count_ratio$type, levels = c("1M_w_ne", "1M_wo_ne", "nearest", "abc_w_ne", "abc_wo_ne", "abc_1M_w_ne", "abc_1M_wo_ne"))
colnames(res_mid) <- c("element", "species", "median_cre_number")
res_mid$element <- gsub("all_", "", res_mid$element)
res_mid <- res_mid[res_mid$species %in% species, ]
unique(anno_count_ratio$type)
anno_tmp <- anno_count_ratio[anno_count_ratio$type %in% c("abc_1M_w_ne", "abc_1M_wo_ne", "abc_wo_ne", "abc_w_ne"), ]
anno_tmp <- anno_tmp[anno_tmp$species %in% species, ]
anno_tmp <- anno_tmp %>%
    group_by(species, element) %>%
    summarize(abc_proportion = sum(ratio))
dim(anno_tmp)
dim(res_mid)
abc_type <- c("abc_1M_w_ne", "abc_1M_wo_ne", "abc_w_ne", "abc_wo_ne")
res_mid <- merge(res_mid, anno_tmp, by = c("species", "element"))
tapply(res_mid$median_cre_number[res_mid$type %in% abc_type], res_mid$element[res_mid$type %in% abc_type], summary)
anno_point <- ggplot(res_mid, aes(x = median_cre_number, y = abc_proportion / 100, fill = element)) +
    geom_point(size = 1.7, shape = 24, color = "transparent") +
    theme_classic() +
    # scale_color_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    scale_fill_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5), limits = c(0, 0.5)) +
    scale_x_continuous(breaks = c(0, 30, 60)) +
    labs(x = "CRE number per gene", y = "ABC annotation proportion", fill = "Element") +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.4), text = element_text(size = 8)) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 0.4)) +
    theme(axis.line = element_line(color = "black", size = 0.4))
cre_num_box <- ggplot(res_mid, aes(y = element, x = median_cre_number, fill = element)) +
    geom_boxploth(size = 0.2, outlier.size = 0.1) +
    theme_dendrogram() +
    scale_fill_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    theme(legend.position = "none", text = element_text(size = 8)) +
    # labs(y = "Element") +
    theme(axis.text.y = element_text(size = 6, color = "black"), axis.line.y = element_line(color = "black", size = 0.4)) +
    theme(axis.ticks.y = element_line(color = "black", size = 0.4))
anno_box <- ggplot(res_mid, aes(y = abc_proportion / 100, x = element, fill = element)) +
    geom_boxplot(linewidth = 0.2, outlier.size = 0.1) +
    theme_dendrogram() +
    scale_fill_manual(values = c("enhancer" = "#ef4343", "promoter" = "#73b8d5")) +
    theme(legend.position = "none", text = element_text(size = 8)) +
    # labs(x = "Element") +
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 0.5)) +
    theme(axis.line.x = element_line(color = "black", size = 0.4)) +
    theme(axis.ticks.x = element_line(color = "black", size = 0.4))
anno_p <- anno_point %>%
    insert_top(cre_num_box, height = 0.1) %>%
    insert_right(anno_box, width = 0.1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/anno_point.pdf", anno_p, width = 11, height = 8, dpi = 300, units = "cm")

## violin plot of trun_ratio in rep_ele_ratio
library(Rmisc)
rep_ele_ratio <- read.csv("/media/Data/zhangz/chip/analysis/summary2/rep_ele/rep_ele_ratio.csv", header = TRUE)
rep_ele_ratio <- rep_ele_ratio[rep_ele_ratio$species %in% species, ]
rep_ratio_summary <- summarySE(rep_ele_ratio, measurevar = "trun_ratio", groupvars = c("tis", "ele"))
rep_ratio_violin <- ggplot(rep_ele_ratio, aes(x = tis, y = trun_ratio, fill = ele)) +
    geom_split_violin(color = "white", linewidth = 0.1, alpha = 0.8) +
    # geom_boxplot(width = 0.2, alpha = 0.5, position = position_dodge(0.3)) +
    geom_point(data = rep_ratio_summary, aes(x = tis, y = trun_ratio), pch = 19, position = position_dodge(0.5), size = 1) +
    geom_errorbar(data = rep_ratio_summary, aes(x = tis, ymin = trun_ratio - se, ymax = trun_ratio + se), width = 0.05, position = position_dodge(0.5)) +
    # facet_grid(tis ~ .) +
    theme_classic() +
    scale_fill_manual(values = c("#e83747", "#a9dbdc"), labels = c("enhancer" = "Enhancer", "promoter" = "Promoter")) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Repeat element proportion", fill = "") +
    theme(legend.position = "top", text = element_text(size = 8)) +
    theme(plot.title = element_blank()) +
    theme(legend.key.size = unit(0.5, "cm")) +
    theme(axis.title.x = element_blank())
ggsave("/media/Data/zhangz/chip/analysis/summary2/rep_ele/fig/rep_ratio_violin2.pdf", rep_ratio_violin, width = 6, height = 6, units = "cm")

# for supplement, compare ratio of trun_ratio and raw_ratio in rep_ele_ratio
rep_ele_ratio_trun <- melt(rep_ele_ratio[, c("species", "tis", "ele", "raw_ratio", "trun_ratio")], id.vars = c("species", "tis", "ele"), variable.name = "ratio_type", value.name = "ratio")

rep_ratio_summary2 <- summarySE(rep_ele_ratio_trun, measurevar = "ratio", groupvars = c("tis", "ele", "ratio_type"))
rep_ratio_trun_raw_violin_e <- ggplot(rep_ele_ratio_trun[rep_ele_ratio_trun$ele == "enhancer", ], aes(x = tis, y = ratio, fill = ratio_type)) +
    geom_split_violin(color = "white", linewidth = 0.1, alpha = 0.8) +
    geom_point(data = rep_ratio_summary2[rep_ratio_summary2$ele == "enhancer", ], aes(x = tis, y = ratio), pch = 19, position = position_dodge(0.5), size = 0.2) +
    geom_errorbar(data = rep_ratio_summary2[rep_ratio_summary2$ele == "enhancer", ], aes(x = tis, ymin = ratio - se, ymax = ratio + se), width = 0.05, position = position_dodge(0.5)) +
    theme_classic() +
    scale_fill_manual(values = c("raw_ratio" = "#ef4343", "trun_ratio" = "#c01d2e"), labels = c("raw_ratio" = "All enhancers", "trun_ratio" = "Filtered enhancers")) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Repeat element proportion", fill = "") +
    theme(legend.position = "top", text = element_text(size = 8)) +
    theme(plot.title = element_blank()) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(axis.title.x = element_blank()) +
    theme(legend.box.background = element_rect(color = "transparent", size = 0.1))
rep_ratio_trun_raw_violine_p <- ggplot(rep_ele_ratio_trun[rep_ele_ratio_trun$ele == "promoter", ], aes(x = tis, y = ratio, fill = ratio_type)) +
    geom_split_violin(color = "white", linewidth = 0.1, alpha = 0.8) +
    geom_point(data = rep_ratio_summary2[rep_ratio_summary2$ele == "promoter", ], aes(x = tis, y = ratio), pch = 19, position = position_dodge(0.5), size = 0.2) +
    geom_errorbar(data = rep_ratio_summary2[rep_ratio_summary2$ele == "promoter", ], aes(x = tis, ymin = ratio - se, ymax = ratio + se), width = 0.05, position = position_dodge(0.5)) +
    theme_classic() +
    scale_fill_manual(values = c("raw_ratio" = "#73b8d5", "trun_ratio" = "#02263e"), labels = c("raw_ratio" = "All promoters", "trun_ratio" = "Filtered promoters")) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Repeat element proportion", fill = "") +
    theme(legend.position = "top", text = element_text(size = 8)) +
    theme(plot.title = element_blank()) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(axis.title.x = element_blank()) +
    theme(legend.box.background = element_rect(color = "transparent", size = 0.1))
rep_raw_ratio_p <- rep_ratio_trun_raw_violin_e %>% insert_bottom(rep_ratio_trun_raw_violine_p, height = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/rep_ele/fig/rep_ratio_trun_raw_violin.pdf", rep_raw_ratio_p, width = 11, height = 8, units = "cm")

# plot of repeat element proportion in each species

# prepare data
# genome_stats <- data.frame(
#     species = c("Atelerix_albiventris", "Petaurus_breviceps", "Mus_musculus"),
#     scaffold_n50 = c(0.1, 0.1, 0.1), genome_size = c(0.1, 0.1, 0.1),
#     contig_n50 = c(0.1, 0.1, 0.1), repeat_ratio = c(0.1, 0.1, 0.1),
#     gene_num = c(0.1, 0.1, 0.1)
# )
assembly <- read.csv("/media/Data/zhangz/chip/genomes/quality/report.tsv", header = TRUE, sep = "\t")
for (sp in species) {
    align_cov[align_cov$Genome == sp, "N50"] <- assembly[assembly$Assembly == "N50", sp]
    align_cov[align_cov$Genome == sp, "N_per_100_kbp"] <- assembly[assembly$Assembly == "# N's per 100 kbp", sp]
}
species <- c(
    "Ovis_aries",
    "Bos_taurus", "Neophocaena_asiaeorientalis",
    "Sus_scrofa", "Lama_glama", "Mustela_putorius",
    "Canis_lupus", "Felis_catus", "Equus_asinus",
    "Equus_caballus", "Rhinolophus_pusillus", "Myotis_chinensis",
    "Atelerix_albiventris", "Mus_musculus", "Rattus_norvegicus",
    "Cavia_porcellus", "Oryctolagus_cuniculus", "Macaca_fascicularis",
    "Rhinopithecus_roxellana", "Tupaia_belangeri",
    "Procavia_capensis", "Petaurus_breviceps"
)

## busco plot
busco <- data.frame(
    species = rep(c("Atelerix_albiventris", "Petaurus_breviceps"), 6),
    cat = c("C", "C", "S", "S", "D", "D", "F", "F", "M", "M", "n", "n"),
    value = c(72.2, 84.8, 70.4, 81.8, 1.8, 3.0, 16.4, 7.3, 11.4, 7.9, 3354, 3354)
)

library(ggstance)
# stacked bar plot
busco$cat <- factor(busco$cat, levels = c("M", "F", "D", "S", "C"))
busco_bar <- ggplot(busco[busco$cat %in% c("S", "D", "F", "M"), ], aes(x = value, y = species, fill = cat)) +
    geom_barh(position = "stackv", stat = "identity", linewidth = 0.5, color = "black", width = 0.5) +
    theme_classic() +
    labs(x = "proportions", fill = "BUSCO category") +
    scale_fill_manual(
        values = c("S" = "#b8f4b8", "D" = "#D3E9C5", "F" = "#5C9388", "M" = "#FFC8BE"),
        labels = c("S" = "Complex (C) and single copy (S)", "D" = "Complex (C) and single copy (S)", "F" = "Fragmented (F)", "M" = "Missing (M)")
    ) +
    scale_x_continuous(labels = scales::percent, n.breaks = 3) +
    scale_y_discrete(labels = c("Atelerix_albiventris" = "Four-toed hedgehog", "Petaurus_breviceps" = "Sugar glider")) +
    annotate(geom = "text", x = 0.1, y = "Atelerix_albiventris", label = "C: 72.2% (S: 70.4%, D: 1.8%), F: 16.4%, M: 11.4%, n: 3354", hjust = 0) +
    annotate(geom = "text", x = 0.1, y = "Petaurus_breviceps", label = "C: 84.8% (S: 81.8%, D: 3.0%), F: 7.3%, M: 7.9%, n: 3354", hjust = 0) +
    theme(legend.position = "top", axis.title.y = element_blank()) +
    theme(text = element_text(size = 8), axis.text = element_text(size = 10), axis.title = element_text(size = 10))
# guides(fill = guide_legend(nrow = 2))
busco_bar
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/busco_bar.pdf", busco_bar, width = 9, height = 6)

busco_bar2 <- ggplot(busco[busco$cat %in% c("S", "D", "F", "M"), ], aes(x = species, y = value, fill = cat)) +
    geom_bar(position = "stack", stat = "identity", linewidth = 0.5, color = "black") +
    theme_classic() +
    labs(y = "proportions", fill = "BUSCO category") +
    scale_fill_manual(
        values = c("S" = "#b8f4b8", "D" = "#D3E9C5", "F" = "#5C9388", "M" = "#FFC8BE"),
        labels = c("S" = "Complex (C) and single copy (S)", "D" = "Complex (C) and single copy (S)", "F" = "Fragmented (F)", "M" = "Missing (M)")
    ) +
    scale_x_discrete(labels = c("Atelerix_albiventris" = "Four-toed hedgehog", "Petaurus_breviceps" = "Sugar glider")) +
    annotate(geom = "text", x = "Atelerix_albiventris", y = 10, label = "C: 72.2% (S: 70.4%, D: 1.8%), F: 16.4%, M: 11.4%, n: 3354", hjust = 0) +
    annotate(geom = "text", x = "Petaurus_breviceps", y = 10, label = "C: 84.8% (S: 81.8%, D: 3.0%), F: 7.3%, M: 7.9%, n: 3354", hjust = 0) +
    theme(legend.position = "top", axis.title.x = element_blank()) +
    theme(text = element_text(size = 8))



n50_genen <- data.frame(
    species = rep(c("Atelerix_albiventris", "Petaurus_breviceps"), 2),
    term = c("N50", "N50", "Gene number", "Gene number"),
    value = c(115005015, 363524205, 25035, 23815)
)

n50_genen$value_million <- n50_genen$value / 1e6
n50_bar <- ggplot(n50_genen[n50_genen$term == "N50", ], aes(x = value_million, y = species)) +
    geom_barh(position = "stackv", stat = "identity", size = 0.5, fill = "#bac6ed") +
    theme_classic() +
    labs(x = "Contig N50 (Million)") +
    # scale_fill_manual(values = c("#b8f4b8", "#b0e9e9")) +
    scale_y_discrete(labels = c("Atelerix_albiventris" = "Four-toed hedgehog", "Petaurus_breviceps" = "Sugar glider")) +
    scale_x_continuous(labels = scales::comma) +
    annotate(geom = "text", x = 10, y = "Atelerix_albiventris", label = "115.005 M", hjust = 0) +
    annotate(geom = "text", x = 10, y = "Petaurus_breviceps", label = "363.524 M", hjust = 0) +
    theme(legend.position = "top", axis.title.y = element_blank()) +
    theme(text = element_text(size = 8))
n50_bar
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/n50_bar.pdf", n50_bar, width = 6, height = 2)

n50_genen$value_kilo <- n50_genen$value / 1e3
genen_bar <- ggplot(n50_genen[n50_genen$term == "Gene number", ], aes(x = value_kilo, y = species)) +
    geom_barh(position = "stackv", stat = "identity", width = 0.5, fill = "#A8B0C3", color = "black", linewidth = 0.5) +
    theme_classic() +
    labs(x = "Non-redundant gene number (Thousand)") +
    # scale_fill_manual(values = c("#b8f4b8", "#b0e9e9")) +
    scale_y_discrete(labels = c("Atelerix_albiventris" = "Four-toed hedgehog", "Petaurus_breviceps" = "Sugar glider")) +
    scale_x_continuous(labels = scales::comma, n.breaks = 3) +
    annotate(geom = "text", x = 1, y = "Atelerix_albiventris", label = "25,035", hjust = 0) +
    annotate(geom = "text", x = 1, y = "Petaurus_breviceps", label = "23,815", hjust = 0) +
    theme(legend.position = "top", axis.title.y = element_blank()) +
    theme(text = element_text(size = 8), axis.text = element_text(size = 10), axis.title = element_text(size = 10))
genen_bar
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/genen_bar.pdf", genen_bar, width = 9, height = 6)


# functions
se <- function(x, model) {
    vi <- vcov(model)[1, 1] + x * vcov(model)[1, 2] * 1 + (1 * vcov(model)[2, 1] + x * vcov(model)[2, 2]) * x
    se <- sqrt(vi)
    return(se)
}

pgls_pmax_model <- function(var1, var2, data_df, tree, model, lg = TRUE) {
    data_df <- data_df[!is.na(data_df[[var1]]) & !is.na(data_df[[var2]]), ]
    spp <- data_df$species
    tree <- keep.tip(tree, spp)
    p_list <- adply(spp, 1, function(sp) {
        data_sp <- data_df[data_df$species != sp, ]
        tree_sp <- drop.tip(tree, sp)
        if (lg == TRUE) {
            formula_str <- paste("log(", var2, ") ~ log(", var1, ")", sep = "")
        } else {
            formula_str <- paste(var2, " ~ ", var1, sep = "")
        }
        if (model == "BM") {
            res <- gls(as.formula(formula_str), correlation = corBrownian(phy = tree_sp), data = data_sp, method = "ML", na.action = na.omit)
        } else if (model == "OU") {
            res <- gls(as.formula(formula_str), correlation = corMartins(1, phy = tree_sp), data = data_sp, method = "ML", na.action = na.omit)
        }
        res_summary <- summary(res)
        return(data.frame(species = sp, p_value = res_summary$tTable[2, 4]))
    })
    return(c(max(p_list$p_value), p_list$species[which.max(p_list$p_value)]))
}

div_coverage_pred <- data.frame(x = align_cov$div_time, pred = )
coverage_n_point <- ggplot(align_cov, aes(x = N_per_100_kbp, y = coverage)) +
    geom_point() +
    theme_classic() +
    xlab("Number of N's per 100 kbp") +
    ylab("Genome alignment coverage") +
    theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70")) +
    geom_smooth(method = "lm", se = FALSE, color = "red")
ggsave("/media/Data/zhangz/chip/analysis/summary2/align_coverage.pdf", coverage_n_point, width = 7, height = 7, units = "cm")

tre <- read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
tre$tip.label <- gsub("Macaca_mulatta", "Macaca_fascicularis", tre$tip.label)
tre <- ape::keep.tip(tre, gsub("Macaca_mulatta", "Macaca_fascicularis", species))
setwd("/media/Data/zhangz/chip/genomes/halcoverage")
align_cov <- read.csv("/media/Data/zhangz/chip/genomes/halcoverage/coverage.log", header = TRUE, sep = ",")
align_cov <- align_cov[align_cov$Genome %in% gsub("Macaca_mulatta", "Macaca_fascicularis", species), ]

# sum row
align_cov$sum <- rowSums(align_cov[, -1])
align_cov$coverage <- align_cov$sum / max(align_cov$sum)
# align_cov <- align_cov[align_cov$Genome %in% species, ]
div_time <- ape::cophenetic.phylo(tre)
for (sp in gsub("Macaca_mulatta", "Macaca_fascicularis", species)) {
    align_cov[align_cov$Genome == sp, "div_time"] <- div_time[sp, "Mus_musculus"] / 2
}

# colors
order_colors <- c("Artiodactyla" = "#3682be", "Carnivora" = "#45a776", "Perissodactyla" = "#f05330", "Chiroptera" = "#eed777", "Eulipotyphla" = "#38cb7d", "Rodentia" = "#334f65", "Lagomorpha" = "#ddae33", "Primates" = "#b3974e", "Scandentia" = "#844bb3", "Hyracoidea" = "#93c555", "Diprotodontia" = "#5f6694")

# assign orders
sp_info <- read.csv("/media/Data/zhangz/chip/scripts/info/species.csv", header = TRUE)
align_cov$Order <- sp_info$order[match(align_cov$Genome, sp_info$species)]
align_cov$Order[align_cov$Genome == "Macaca_fascicularis"] <- "Primates"

align_cov$svg_path <- paste0("/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/", align_cov$Genome, ".svg")
align_cov$svg_path[which(align_cov$Genome %in% c("Myotis_chinensis", "Myotis_ricketti", "Rhinolophus_ferrumequinum", "Rhinolophus_pusillus"))] <- "/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/Hipposideros_larvatus.svg"
align_cov$svg_path[align_cov$Genome == "Macaca_fascicularis"] <- "/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/Macaca_mulatta.svg"
for (svg in align_cov$svg_path) {
    if (!file.exists(svg)) {
        print(svg)
    }
}

div_coverage_gls <- gls(coverage ~ div_time, data = align_cov, method = "ML")
div_coverage_gls
div_coverage_r2 <- R2_lik(mod = div_coverage_gls)
div_coverage_point <- ggplot(align_cov[align_cov$Genome != "Mus_musculus", ], aes(x = div_time, y = coverage)) +
    geom_point(aes(fill = Order, color = Order), size = 1, shape = 24) +
    geom_smooth(method = "lm", color = "#ef4343", linewidth = 0.4) +
    theme_classic() +
    scale_color_manual(values = order_colors) +
    scale_fill_manual(values = order_colors) +
    scale_x_continuous(breaks = c(10, 80, 160)) +
    geom_image(image = align_cov$svg_path[align_cov$Genome != "Mus_musculus"], size = 0.05, height = 0.03, alpha = 0.5) +
    xlab("Species divergence time (from mouse)") +
    ylab("Genome alignment coverage") +
    theme(text = element_text(size = 8), se = TRUE) +
    annotate("text", x = 80, y = 0.7, label = paste(
        "R^2 = ", round(div_coverage_r2, 2),
        ", p value =", round(summary(div_coverage_gls)$tTable[2, 4], 2)
    ), size = 2) +
    theme(
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8),
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.line = element_line(colour = "black"), legend.position = "none"
    )

ggsave("/media/Data/zhangz/chip/analysis/summary2/align_coverage_div.pdf", div_coverage_point, width = 7, height = 7, units = "cm")

assembly <- read.csv("/media/Data/zhangz/chip/genomes/quality/report.tsv", header = TRUE, sep = "\t")
for (sp in gsub("Macaca_mulatta", "Macaca_fascicularis", species)) {
    align_cov[align_cov$Genome == sp, "N50"] <- assembly[assembly$Assembly == "N50", sp]
    align_cov[align_cov$Genome == sp, "N_per_100_kbp"] <- assembly[assembly$Assembly == "# N's per 100 kbp", sp]
}

# , "Rattus_norvegicus", "Petaurus_breviceps"
cov_n50_point1 <- ggplot(align_cov[!align_cov$Genome %in% c("Mus_musculus"), ], aes(x = N50, y = coverage, color = Order, fill = Order)) +
    geom_point(shape = 24, size = 1) +
    theme_classic() +
    scale_color_manual(values = order_colors) +
    scale_fill_manual(values = order_colors) +
    theme(legend.position = "none") +
    geom_image(image = align_cov$svg_path[!align_cov$Genome %in% c("Mus_musculus")], size = 0.05, height = 0.03, alpha = 0.7) +
    scale_x_continuous(breaks = c(0, 200000000, 400000000), labels = scales::comma, limits = c(0, 400000000))
cov_n50_point2 <- cov_n50_point1 +
    scale_y_continuous(limits = c(0.15, 0.35), n.breaks = 3, labels = scales::percent) +
    scale_x_continuous(limits = c(0, 150000000), n.breaks = 3, labels = scales::comma) +
    theme(axis.title = element_blank())
library(patchwork)
cov_n50_point <- cov_n50_point1 +
    inset_element(cov_n50_point2, left = 0.4, bottom = 0.4, right = 1, top = 1)

div_coverage_gls <- gls(coverage ~ div_time, data = align_cov, method = "ML")
# 计算拟合模型的残差
residuals_fit <- residuals(div_coverage_gls)

# 将残差添加到原始数据框中
align_cov$residuals <- residuals_fit

# 打印残差
print(align_cov$residuals)
cov_res_n50_point <- ggplot(align_cov[!align_cov$Genome %in% c("Mus_musculus", "Petaurus_breviceps"), ], aes(x = N50, y = residuals)) +
    geom_point(aes(color = Order, fill = Order), shape = 24, size = 2) +
    theme_classic() +
    scale_color_manual(values = order_colors) +
    scale_fill_manual(values = order_colors) +
    theme(legend.position = "none") +
    labs(x = "Contig N50 (million)", y = "Genome alignment residuals") +
    geom_image(image = align_cov$svg_path[!align_cov$Genome %in% c("Mus_musculus", "Petaurus_breviceps")], size = 0.05, height = 0.05, alpha = 0.5) +
    scale_x_continuous(breaks = c(0, 100000000, 200000000), labels = scales::comma, limits = c(0, 200000000)) +
    scale_y_continuous(breaks = c(-0.1, 0, 0.1), limits = c(-0.2, 0.2)) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), axis.line = element_line(color = "black", size = 0.4), axis.ticks = element_line(color = "black", size = 0.4))
ggsave("/media/Data/zhangz/chip/analysis/summary2/align_coverage_residuals_n50.pdf", cov_res_n50_point, width = 5, height = 5)

cov_n50_point <- ggplot(align_cov[!align_cov$Genome %in% c("Mus_musculus", "Petaurus_breviceps"), ], aes(x = N50, y = coverage, color = Order, fill = Order)) +
    geom_point(shape = 24, size = 1) +
    theme_classic() +
    scale_color_manual(values = order_colors) +
    scale_fill_manual(values = order_colors) +
    theme(legend.position = "none") +
    geom_image(image = align_cov$svg_path[!align_cov$Genome %in% c("Mus_musculus", "Petaurus_breviceps")], size = 0.05, height = 0.03, alpha = 0.7) +
    scale_x_continuous(breaks = c(0, 75000000, 150000000), labels = scales::comma, limits = c(0, 150000000))
ggsave("/media/Data/zhangz/chip/analysis/summary2/align_coverage_n50.pdf", cov_n50_point, width = 7, height = 7)

coverage_n50_gls_ou <- gls(coverage ~ N50, data = align_cov, method = "ML", correlation = corMartins(1, phy = tre))
coverage_n50_gls_ou_r2 <- R2_lik(mod = coverage_n50_gls_ou)
coverage_n50_gls_bm <- gls(coverage ~ N50, data = align_cov, method = "ML", correlation = corBrownian(1, phy = tre))
coverage_n50_gls_bm_r2 <- R2_lik(mod = coverage_n50_gls_bm)
cor.test(pic(align_cov$coverage, phy = tre), pic(align_cov$N50, phy = tre), method = "spearman")
coverage_n_gls_ou <- gls(coverage ~ N_per_100_kbp, data = align_cov, method = "ML", correlation = corMartins(1, phy = tre))
coverage_n_gls_ou_r2 <- R2_lik(mod = coverage_n_gls_ou)
coverage_n_gls_bm <- gls(coverage ~ N_per_100_kbp, data = align_cov, method = "ML", correlation = corBrownian(1, phy = tre))
coverage_n_gls_bm_r2 <- R2_lik(mod = coverage_n_gls_bm)
log_coverage_n_gls_ou <- gls(coverage ~ 1 / N_per_100_kbp, data = align_cov, method = "ML", correlation = corMartins(1, phy = tre))

ggplot(pic_matrix, aes(x = N_per_100_kbp, y = coverage)) +
    geom_point(size = 1, shape = 24)

power_model <- function(x, a, b) {
    a / x^b
}

start_params <- list(a = 1, b = 1)

fit <- gnls(coverage ~ power_model(N_per_100_kbp, a, b), data = align_cov, start = start_params, correlation = corMartins(1, phy = tre))
summary(fit)

cor.test(pic(align_cov$coverage, phy = tre), pic(align_cov$N_per_100_kbp, phy = tre), method = "spearman")
cor.test(pic(align_cov$coverage, phy = tre), pic(align_cov$N_per_100_kbp, phy = tre), method = "pearson")
cor.test(align_cov$coverage, align_cov$N_per_100_kbp, method = "spearman")

# correlation heat map of coverage, N50, N_per_100_kbp
pic_matrix <- data.frame(coverage = pic(align_cov$coverage, phy = tre), N50 = pic(align_cov$N50, phy = tre), N_per_100_kbp = pic(align_cov$N_per_100_kbp, phy = tre))
cor_matrix <- cor(pic_matrix, use = "pairwise.complete.obs", method = "spearman")
library(Hmisc)
corr_p <- rcorr(as.matrix(pic_matrix), type = "spearman")
library(corrplot)
library(grid)
library(cowplot)
library(ggsci)
grid.newpage()
pdf("/media/Data/zhangz/chip/analysis/summary2/align_coverage_cor_heatmap.pdf")
corrplot(cor_matrix,
    method = "ellipse",
    type = "upper", tl.col = "black", tl.cex = 0.8, tl.srt = 45, tl.pos = "lt",
    p.mat = corr_p$P, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig",
    col = colorRampPalette(c("#1d355b", "white", "#c01d2e"))(200)
)
corrplot(cor_matrix,
    method = "number",
    type = "lower", add = TRUE, tl.col = "n", tl.cex = 0.8, tl.pos = "n",
    p.mat = corr_p$P, sig.level = 0.05,
    col = colorRampPalette(c("#1d355b", "white", "#c01d2e"))(200)
)
dev.off()

## still need a distribution plot of cre alignability
sp <- "Mus_musculus"
ele_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_align22.csv"), header = TRUE, sep = ",")

# bottom part of enhancer distribution
align_dis_enhancer_plot_bottom <- ggplot(ele_df[which(ele_df$align22 > 0 & ele_df$element == "enhancer"), ], aes(x = align22)) +
    geom_histogram(binwidth = 1, fill = "#ef4343", color = "black", linewidth = 0.01) +
    theme_classic() +
    scale_x_continuous(breaks = c(1, 10, 21)) +
    scale_y_continuous(n.breaks = 3) +
    coord_cartesian(ylim = c(0, 8000)) +
    xlab(paste0("Number of species that CREs sequence in\n", gsub("_", " ", sp), " can be aligned to")) +
    theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70"))
# top part of enhancer distribution
align_dis_enhancer_plot_top <- ggplot(ele_df[which(ele_df$align22 > 0 & ele_df$element == "enhancer"), ], aes(x = align22)) +
    geom_histogram(binwidth = 1, fill = "#ef4343", color = "black", linewidth = 0.01) +
    theme_classic() +
    scale_x_continuous(breaks = c(1, 10, 21)) +
    scale_y_continuous(n.breaks = 2) +
    coord_cartesian(ylim = c(20000, 21000)) +
    xlab(paste0("Number of species that CREs sequence in\n", gsub("_", " ", sp), " can be aligned to")) +
    theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70")) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
align_dis_enhancer_plot <- plot_grid(align_dis_enhancer_plot_top, align_dis_enhancer_plot_bottom, nrow = 2, align = "v", axis = "l", rel_heights = c(0.5, 2))

# ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_enhancer_align_dis.pdf"), align_dis_enhancer_plot, width = 7, height = 7, units = "cm")
align_dis_promoter_plot_bottom <- ggplot(ele_df[which(ele_df$align22 > 0 & ele_df$element == "promoter"), ], aes(x = align22)) +
    geom_histogram(binwidth = 1, fill = "#73b8d5", color = "black", size = 0.01) +
    theme_classic() +
    scale_x_continuous(breaks = c(1, 10, 21), position = "top") +
    scale_y_reverse(n.breaks = 3) +
    coord_cartesian(ylim = c(4000, 0)) +
    xlab(paste0("Number of species that CREs sequence in\n", gsub("_", " ", sp), " can be aligned to")) +
    theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70"))
align_dis_promoter_plot_top <- ggplot(ele_df[which(ele_df$align22 > 0 & ele_df$element == "promoter"), ], aes(x = align22)) +
    geom_histogram(binwidth = 1, fill = "#73b8d5", color = "black", size = 0.01) +
    theme_classic() +
    scale_x_continuous(breaks = c(1, 10, 21), position = "top") +
    scale_y_reverse(n.breaks = 2) +
    coord_cartesian(ylim = c(7000, 6000)) +
    xlab(paste0("Number of species that CREs sequence in\n", gsub("_", " ", sp), " can be aligned to")) +
    theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70")) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
    )
align_dis_promoter_plot <- plot_grid(align_dis_promoter_plot_bottom, align_dis_promoter_plot_top, nrow = 2, align = "v", axis = "l", rel_heights = c(2, 0.5))
# ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_promoter_align_dis.pdf"), align_dis_promoter_plot, width = 7, height = 7, units = "cm")
align_dis_plot <- plot_grid(align_dis_enhancer_plot_top, align_dis_enhancer_plot_bottom, align_dis_promoter_plot_bottom, align_dis_promoter_plot_top, nrow = 4, ncol = 1, align = "v", axis = "l", rel_heights = c(0.5, 2, 2, 0.5))
ggsave("/media/Data/zhangz/chip/analysis/summary2/mouse_align_div.pdf", align_dis_plot, width = 7, height = 10, units = "cm")

## we need to show Atelerix_albiventris because the genome quality is not so good as others
sp <- "Atelerix_albiventris"
ele_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_align22.csv"), header = TRUE, sep = ",")

# enhancer distribution
align_dis_enhancer_plot <- ggplot(ele_df[which(ele_df$align22 > 0 & ele_df$element == "enhancer"), ], aes(x = align22)) +
    geom_histogram(binwidth = 1, fill = "#ef4343", color = "black", linewidth = 0.01) +
    theme_classic() +
    scale_x_continuous(breaks = c(1, 10, 21)) +
    scale_y_continuous(n.breaks = 3) +
    xlab(paste0("Number of species that CREs sequence in\n", gsub("_", " ", sp), " can be aligned to")) +
    theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70"))

# promoter distribution
align_dis_promoter_plot <- ggplot(ele_df[which(ele_df$align22 > 0 & ele_df$element == "promoter"), ], aes(x = align22)) +
    geom_histogram(binwidth = 1, fill = "#73b8d5", color = "black", size = 0.01) +
    theme_classic() +
    scale_x_continuous(breaks = c(1, 10, 21), position = "top") +
    scale_y_reverse(n.breaks = 3) +
    xlab(paste0("Number of species that CREs sequence in\n", gsub("_", " ", sp), " can be aligned to")) +
    theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70"))

align_dis_plot <- plot_grid(align_dis_enhancer_plot, align_dis_promoter_plot, nrow = 2, align = "v", axis = "l", rel_heights = c(1, 1))
ggsave("/media/Data/zhangz/chip/analysis/summary2/hedgehog_align_dis.pdf", align_dis_plot, width = 7, height = 9, units = "cm")
