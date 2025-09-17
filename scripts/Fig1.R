## 　draw plots of Fig1

## 　set up the environment
library(plyr)
library(treedataverse)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(aplot)
library(reshape2)
library(tibble)
library(ggridges)
library(tidyverse)
library(doMC)
doMC::registerDoMC(4)

## load the data and base information including tree
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
tissue_color <- c("Brain" = "#FFD2A8", "Kidney" = "#84C990", "Liver" = "#7E8BB4")

all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/sum_all/all_element.csv", header = TRUE)

tree <- read.newick("/media/Data/zhangz/chip/scripts2/info/sps.nwk", node.label = "label")
# tree <- keep.tip(tree, c("Ovis_aries", "Equus_caballus", "Rhinolophus_pusillus", "Rhinopithecus_roxellana", "Sus_scrofa", "Felis_catus"))
# tt_tre <- ggtree(tree, layout = "rectangular", size = 0.95) +
#     theme_tree2() +
#     scale_x_continuous(limits = c(-0.1, 94.1), breaks = c(0, 94), labels = c(94, 0), expand = c(0, 0))
# ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/tree/minitree.pdf", tt_tre, height = 12, width = 2, units = "cm")
base_tp2 <- ggtree(keep.tip(tree, species), layout = "rectangular", size = 0.95) +
    theme_tree2() +
    scale_x_continuous(limits = c(-0.1, 160.1), breaks = c(0, 15, 94, 136.92, 160), labels = c(160, 145, 66, 23.08, 0), expand = c(0, 0))
sp_df <- data.frame(species = species, label = gsub("_", " ", species))
sp_p <- ggplot(sp_df, aes(x = 0, y = species)) +
    geom_text(aes(label = label), size = 4, fontface = "italic") +
    theme_void()

## draw plots of CRE length and fc
all <- all[, c(14, 1:13, 15:26)]
all <- all[all$species %in% species, ]
all$tissue <- factor(all$tissue, levels = c("Brain", "Kidney", "Liver"))
all$motif_density <- all$motif_num / all$length
all_enh <- all[all$element == "enhancer", ]
all_pro <- all[all$element == "promoter", ]
all_enh$species <- factor(all_enh$species, levels = tree$tip.label)

enh_length <- ggplot(all_enh[all_enh$length < 1000, ], aes(x = length, y = species, fill = tissue)) +
    geom_density_ridges(alpha = 0.7, scale = 1, rel_min_height = 0.01, linewidth = 0.2) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_blank(),
        legend.position = "none", plot.title = element_text(size = 8, hjust = 0.5),
        axis.line.x = element_line(colour = "grey70")
    ) +
    labs(title = "Enhance\nlength") +
    scale_fill_manual(values = tissue_color) +
    theme(axis.line.x = element_line(colour = "black")) +
    scale_x_continuous(limits = c(0, 1000), breaks = c(0, 500, 1000), labels = c(0, 500, 1000), expand = c(0, 0))
pro_length <- ggplot(all_pro[all_pro$length < 1000, ], aes(x = length, y = species, fill = tissue)) +
    geom_density_ridges(alpha = 0.7, scale = 1, rel_min_height = 0.01, linewidth = 0.2) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        axis.line.x = element_line(colour = "grey70")
    ) +
    labs(title = "Promoter\nlength", fill = "Tissue") +
    scale_fill_manual(values = tissue_color) +
    theme(axis.line.x = element_line(colour = "black")) +
    scale_x_continuous(limits = c(0, 1000), breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
enh_fc <- ggplot(all_enh[all_enh$mean_fc < 10, ], aes(x = mean_fc, y = species, fill = tissue)) +
    geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01, linewidth = 0.2) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5), legend.position = "none",
        axis.line.x = element_line(colour = "grey70")
    ) +
    ggtitle("Enhancer\nfoldchange") +
    scale_fill_manual(values = tissue_color) +
    theme(axis.line.x = element_line(colour = "black")) +
    scale_x_continuous(limits = c(0, 10), breaks = c(0, 5, 10), labels = c(0, 5, 10))
pro_fc <- ggplot(all_pro[all_pro$mean_fc < 10, ], aes(x = mean_fc, y = species, fill = tissue)) +
    geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01, linewidth = 0.2) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5), legend.position = "none",
        axis.line.x = element_line(colour = "grey70")
    ) +
    labs(title = "Promoter\nfoldchange", fill = "Tissue") +
    scale_fill_manual(values = tissue_color) +
    theme(axis.line.x = element_line(colour = "black")) +
    scale_x_continuous(limits = c(0, 10), breaks = c(0, 5, 10), labels = c(0, 5, 10))

## draw sub plots of pleiotropy and specificity
specificity_cat <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/specificity_cat.csv", header = TRUE)
specificity_cat_m <- melt(specificity_cat, id.vars = "species")
colnames(specificity_cat_m) <- c("species", "cat", "proportion")
specificity_cat_m <- specificity_cat_m %>%
    separate(cat, into = c("sps", "ts", "ele"), sep = "_") %>%
    unite(sp_t_sp, sps, ts, sep = "_")
specificity_e_plot <- ggplot(specificity_cat_m[specificity_cat_m$ele == "e", ], aes(x = proportion, y = species, fill = sp_t_sp)) +
    geom_barh(stat = "identity", position = "stack") +
    scale_fill_manual(
        values = c("evol_pleio" = "#b23f3f", "evol_ts" = "#ce6868", "sp_ts" = "#ffbebe", "sp_nts" = "#ffd7d7"),
        labels = c("evol_pleio" = "evolutionary pleiotropy", "evol_ts" = "evolutionary tissue-specificity", "sp_ts" = "species & tissue-specific", "sp_nts" = "species-specific pleiotropy")
    ) +
    theme_void() +
    theme(
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        legend.title = element_text(size = 8), plot.title = element_text(size = 8, hjust = 0.5),
        axis.text.x = element_text(size = 8)
    ) +
    scale_x_continuous(n.breaks = 3, labels = scales::percent) +
    labs(title = "Enhancer specificity", fill = "Enhancer specificity")
specificity_p_plot <- ggplot(specificity_cat_m[specificity_cat_m$ele == "p", ], aes(x = proportion, y = species, fill = sp_t_sp)) +
    geom_barh(stat = "identity", position = "stack") +
    scale_fill_manual(
        values = c("evol_pleio" = "#3a6c82", "evol_ts" = "#548195", "sp_ts" = "#99b7c5", "sp_nts" = "#c0d3dc"),
        labels = c("evol_pleio" = "evolutionary pleiotropy", "evol_ts" = "evolutionary tissue-specificity", "sp_ts" = "species & tissue-specific", "sp_nts" = "species-specific pleiotropy")
    ) +
    theme_void() +
    theme(
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        legend.title = element_text(size = 8), plot.title = element_text(size = 8, hjust = 0.5),
        axis.text.x = element_text(size = 8)
    ) +
    scale_x_continuous(n.breaks = 3, labels = scales::percent) +
    labs(title = "Promoter specificity", fill = "Promoter specificity")

## arrange the plots
tree_ridge_plot <- sp_p %>%
    insert_right(specificity_e_plot, width = 1) %>%
    insert_right(specificity_p_plot, width = 1) %>%
    insert_left(base_tp2, width = 1) %>%
    insert_right(enh_length, width = 0.5) %>%
    insert_right(pro_length, width = 0.5) %>%
    insert_right(enh_fc, width = 0.5) %>%
    insert_right(pro_fc, width = 0.5)
tree_ridge_plott <- print(tree_ridge_plot) & theme(legend.position = "bottom") &
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))

## save the whole plot
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/tree/ridge_tree3.pdf", tree_ridge_plott, width = 14, height = 9, dpi = 600)
