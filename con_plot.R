## we have used /media/Data/zhangz/chip/scripts2/con_tissues.py
## to calculate how many elements are EC in different tissues in
## other species, for example ortholog mouse brain enhancer function
## in Ovis aries kidney and so on.
## this script was aimed to draw a heatmap to show the results.

## load libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(doMC)
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
doMC::registerDoMC(cores = 10)

# maybe more

# read tree
sp_tre <- ape::read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
sp_tre <- drop.tip(sp_tre, c("Neophocaena_asiaeorientalis", "Rhinopithecus_roxellana"))

## read in data
func <- read.csv("/media/Data/zhangz/chip/analysis/conservation3/intersect_count.csv",
    header = TRUE, sep = ","
)
sp_info <- read.csv("/media/Data/zhangz/chip/scripts/info/species.csv", header = TRUE, sep = ",")
# assign order according to sp_info
func$order1 <- sp_info$order[match(func$sp1, sp_info$species)]
func$order2 <- sp_info$order[match(func$sp2, sp_info$species)]
func$st1 <- interaction(func$sp1, func$tis1)
func$st2 <- interaction(func$sp2, func$tis2)
func$type <- "EC"

## i want to use interaction of sp1 and tis1 as x-axis, sp2 and tis2 as y-axis
## and the count value as fill color, facet by ele
## both axis should be grouped by sp and the plot will be a heatmap with tree plot
## on the left and top side.
# 检查树是否是等时树
if (!is.ultrametric(sp_tre)) {
    # 如果树不是等时树，将它转换为等时树
    sp_tre <- chronos(sp_tre)
}

# 现在你可以将树转换为层次聚类
hclus <- as.hclust(sp_tre)
hclus <- as.hclust(sp_tre)
vtre <- ggtree(sp_tre, layout = "rectangular") %<+% func[, c("sp2", "tis2")] + scale_y_continuous(expand = c(0, 0)) + geom_tiplab()
htre <- ggtree(sp_tre, layout = "rectangular") %<+% func[, c("sp1", "tis1")] + coord_flip() + scale_x_continuous(expand = c(0, 0)) + geom_tiplab()
ggsave("/media/Data/zhangz/chip/analysis/conservation3/vtre.pdf", vtre, width = 10, height = 20, units = "cm")
ggsave("/media/Data/zhangz/chip/analysis/conservation3/htre.pdf", htre, width = 20, height = 10, units = "cm")
st1 <- get_taxa_name(htre)
st2 <- get_taxa_name(vtre)
func$sp1 <- factor(func$sp1, levels = rev(st1))
func$sp2 <- factor(func$sp2, levels = rev(st2))

ec_heat_map_e <- ggplot(func[func$ele == "enhancer", ], aes(x = interaction(tis1, sp1), y = interaction(tis2, sp2))) +
    geom_tile(aes(fill = count), color = "white") +
    scale_fill_gradient(low = "white", high = "#c45058") +
    guides(x = "axis_nested", y = "axis_nested") +
    # facet_grid(ele ~ .) +
    theme_void() +
    theme(legend.position = "right")
ec_heat_map_p <- ggplot(func[func$ele == "promoter", ], aes(x = interaction(tis1, sp1), y = interaction(tis2, sp2))) +
    geom_tile(aes(fill = (count)), color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    guides(x = "axis_nested", y = "axis_nested") +
    theme_void() +
    theme(legend.position = "right")
# scale_x_dendrogram(hclust = hclus) +
# scale_y_dendrogram(hclust = hclus)
# theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     strip.text.x = element_text(size = 8),
#     strip.text.y = element_text(size = 8),
#     strip.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.line = element_line(color = "black"),
#     legend.position = "none",
# )
# tre_ec_heat <- cowplot::plot_grid(vtre, ec_heat_map, NULL, htre, ncol = 2, rel_widths = c(0.2, 0.8), rel_heights = c(0.8, 0.2))

# test <- ec_heat_map_e + scale_x_dendrogram(hclust = hclus) + scale_y_dendrogram(hclust = hclus)

# strip of tis1
tissue_color <- c("Brain" = "#628F91", "Kidney" = "#9DCD8B", "Liver" = "#DBE79D")
tis_type1 <- ggplot(func, aes(x = interaction(tis1, sp1), y = type, fill = tis1)) +
    geom_tile() +
    theme_void() +
    scale_fill_manual(values = tissue_color) +
    theme(legend.position = "none")
order_colors <- c(
    "Artiodactyla" = "#6C9F9F", "Carnivora" = "#C56161", "Chiroptera" = "#C58E61", "Diprotodontia" = "#D09090",
    "Eulipotyphla" = "#3A7676", "Hyracoidea" = "#8c564b", "Lagomorpha" = "#4D9D4D", "Perissodactyla" = "#7BBD7B",
    "Primates" = "#F0A2A2", "Rodentia" = "#F0C8A7", "Scandentia" = "#87C787", "#FF55a3"
)

tis_type2 <- ggplot(func, aes(x = type, y = interaction(tis2, sp2), fill = tis2)) +
    geom_tile() +
    theme_void() +
    scale_fill_manual(values = tissue_color) #+
# theme(legend.position = "none")
order1 <- ggplot(func, aes(x = interaction(tis1, sp1), y = type, fill = order1)) +
    geom_tile() +
    theme_void() +
    scale_fill_manual(values = order_colors) #+
# theme(legend.position = "none")
order2 <- ggplot(func, aes(x = type, y = interaction(tis2, sp2), fill = order2)) +
    geom_tile() +
    theme_void() +
    scale_fill_manual(values = order_colors) +
    theme(legend.position = "none")
ec_tree_e <- ec_heat_map_e %>%
    insert_right(tis_type2, width = 0.02) %>%
    insert_top(tis_type1, height = 0.02) %>%
    insert_left(order2, width = 0.02) %>%
    insert_bottom(order1, height = 0.02)
ggsave("/media/Data/zhangz/chip/analysis/conservation3/ec_heatmap_e.pdf", ec_tree_e, width = 20, height = 17, units = "cm")
ec_tree_p <- ec_heat_map_p %>%
    insert_right(tis_type2, width = 0.02) %>%
    insert_top(tis_type1, height = 0.02) %>%
    insert_left(order2, width = 0.02) %>%
    insert_bottom(order1, height = 0.02)
ggsave("/media/Data/zhangz/chip/analysis/conservation3/ec_heatmap_p.pdf", ec_tree_p, width = 20, height = 17, units = "cm")
#  %>%
# insert_left(vtre, width = 0.1) %>%
# insert_bottom(htre, height = 0.1)
#     p <- cowplot::plot_grid(vtre, (ec_tree_e), NULL, htre, ncol = 2, rel_widths = c(0.1, 0.9), rel_heights = c(0.9, 0.1))
# test <- cowplot::plot_grid(vtre, ec_heat_map_e, NULL, htre, ncol = 2, rel_widths = c(0.2, 0.8), rel_heights = c(0.8, 0.2)) %>%
#     insert_left(order2, width = 0.05) %>%
#     insert_bottom(order1, height = 0.05) %>%
#     insert_right(tis_type2, width = 0.05) %>%
#     insert_top(tis_type1, height = 0.05)


## for condition of func[, c("sp1", "tis1", "sp2", "ele")],
## /media/Data/zhangz/chip/analysis/{sp1}/compare2/{tis1}/{ele}/{sp1}2{sp2}_{tis1}_{ele}_halper.bed4
## check the number of elements can be found in the other species
## use these numbers to normalize the count value
## and draw a heatmap to show the conservation of elements in different species

func2 <- mdply(func[, c("sp1", "tis1", "sp2", "ele")], function(sp1, tis1, sp2, ele) {
    file <- sprintf("/media/Data/zhangz/chip/analysis/%s/compare2/%s/%s/%s2%s_%s_%s_halper.bed4", sp1, tis1, ele, sp1, sp2, tis1, ele)
    if (file.exists(file)) {
        count <- readLines(file) %>% length()
    } else {
        count <- 0
    }
    return(data.frame(sp1 = sp1, tis1 = tis1, sp2 = sp2, ele = ele, count = count))
}, .parallel = TRUE)
colnames(func2)[5] <- "gc_count"

# add gc_count to func
func$gc_count <- func2$gc_count[match(paste(func$sp1, func$tis1, func$sp2, func$ele), paste(func2$sp1, func2$tis1, func2$sp2, func2$ele))]
func$ec_ratio <- func$count / func$gc_count
write.csv(func, "/media/Data/zhangz/chip/analysis/conservation3/intersect_count.csv", row.names = FALSE)
# draw heatmap
ec_ratio_heat_e <- ggplot(func[func$ele == "enhancer", ], aes(x = interaction(tis1, sp1), y = interaction(tis2, sp2))) +
    geom_tile(aes(fill = ec_ratio), color = "white") +
    scale_fill_gradient(low = "white", high = "#c45058") +
    guides(x = "axis_nested", y = "axis_nested") +
    theme_void() +
    theme(legend.position = "right")
ec_ratio_heat_p <- ggplot(func[func$ele == "promoter", ], aes(x = interaction(tis1, sp1), y = interaction(tis2, sp2))) +
    geom_tile(aes(fill = ec_ratio), color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    guides(x = "axis_nested", y = "axis_nested") +
    theme_void() +
    theme(legend.position = "right")
ec_ratio_tree_e <- ec_ratio_heat_e %>%
    insert_right(tis_type2, width = 0.02) %>%
    insert_top(tis_type1, height = 0.02) %>%
    insert_left(order2, width = 0.02) %>%
    insert_bottom(order1, height = 0.02)
ggsave("/media/Data/zhangz/chip/analysis/conservation3/ec_ratio_heatmap_e.pdf", ec_ratio_tree_e, width = 20, height = 17, units = "cm")
ec_ratio_tree_p <- ec_ratio_heat_p %>%
    insert_right(tis_type2, width = 0.02) %>%
    insert_top(tis_type1, height = 0.02) %>%
    insert_left(order2, width = 0.02) %>%
    insert_bottom(order1, height = 0.02)
ggsave("/media/Data/zhangz/chip/analysis/conservation3/ec_ratio_heatmap_p.pdf", ec_ratio_tree_p, width = 20, height = 17, units = "cm")

# box plot of ec_ratio with outlier point, grouped by func$tissue, facet by func$ele
func$tissue <- (func$tis1 == func$tis2)
ec_violin <- ggplot(func, aes(x = tissue, y = ec_ratio, color = tissue)) +
    # geom_violin() +
    geom_boxplot(width = 0.4, outliers = TRUE) +
    # geom_jitter() +
    geom_hline(yintercept = median(func$ec_ratio[func$tissue == FALSE]), linetype = "dashed", color = "red") +
    facet_grid(ele ~ ., scale = "free") +
    theme_classic() +
    scale_color_manual(values = c("TRUE" = "dark blue", "FALSE" = "grey")) +
    stat_compare_means(method = "wilcox.test", label.x = 1.3) + # label = "p.format",
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.x = element_text(size = 8), strip.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black"), legend.position = "top")
ggsave("/media/Data/zhangz/chip/analysis/conservation3/ec_violin.pdf", ec_violin, width = 10, height = 10, units = "cm")


## try to summarize whether each element in each species and tissue is activited in every other species and tissue
### first, combine rep annotation and adjusted distance of filtered elements
species <- sp_info$species
species <- species[!species %in% c("Hipposideros_larvatus", "Rhinolophus_ferrumequinum", "Myotis_ricketti")]
plyr::a_ply(species, 1, function(sp) {
    rep_anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno_rep.csv"), header = TRUE, sep = ",")
    dis_anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno_new.csv"), header = TRUE, sep = ",")
    dis_col <- which(colnames(rep_anno) %in% c("distanceToTSS", "overlap_distance"))
    rep_anno <- rep_anno[, -dis_col]
    rep_anno <- dplyr::arrange(rep_anno, unique_id, peak)
    dis_anno <- dplyr::arrange(dis_anno, unique_id, peak)
    # unique_id were not unique, try to make it unique
    rep_anno$unique_id <- paste(rep_anno$unique_id, rep_anno$peak, sep = "_")
    dis_anno$unique_id <- paste(dis_anno$unique_id, dis_anno$peak, sep = "_")
    rep_anno <- merge(rep_anno, dis_anno[, c("unique_id", "distanceToTSS", "overlap_distance")], by = "unique_id")
    write.csv(rep_anno, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), row.names = FALSE)
})
plyr::a_ply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    if (ncol(anno) < 100) {
        ec_cols <- expand.grid(species = species[species != sp], tissue = c("Brain", "Kidney", "Liver"))
        ec_cols <- paste(ec_cols$species, ec_cols$tissue, sep = "_")
        ec_cols <- ec_cols[!ec_cols %in% c("Rhinopithecus_roxellana_Brain", "Neophocaena_asiaeorientalis_Kidney", "Neophocaena_asiaeorientalis_Liver")]
        for (col in ec_cols) {
            anno <- dplyr::mutate(anno, !!col := 0)
        }
        condts <- expand.grid(tis1 = unique(anno$tissue), ele = unique(anno$element), sp2 = species[species != sp], tis2 = c("Brain", "Kidney", "Liver"))
        for (i in 1:nrow(condts)) {
            tis1 <- condts$tis1[i]
            ele <- condts$ele[i]
            sp2 <- condts$sp2[i]
            tis2 <- condts$tis2[i]
            print(paste(sp, tis1, ele, sp2, tis2))
            intersect_file <- paste0("/media/Data/zhangz/chip/analysis/", sp, "/compare2/", tis1, "/", ele, "/", sp, "2", sp2, "_", tis1, "_", tis2, "_", ele, "_intersect.bed")
            if (file.exists(intersect_file) && file.size(intersect_file) > 0) {
                itsect <- read.csv(intersect_file, header = FALSE, sep = "\t")
                colnames(itsect) <- c("chr", "start", "end", "peak")
                anno[which(anno$peak %in% itsect$peak & anno$tissue == tis1 & anno$element == ele), paste(sp2, tis2, sep = "_")] <- 1
            }
        }
        write.csv(anno, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), row.names = FALSE)
    }
})
library(ape)
tre <- ape::read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
time_matrix <- cophenetic(tre)
time_matrix[sp, sp2] / 2
count_ec <- function(row, sps) {
    res <- data.frame(s0_0 = 0, s0_1 = 0, s1_0 = 0, s1_1 = 0, act_time = 0)
    tis <- row$tissue
    dt <- tiss[which(tiss != tis)]
    act_time <- c()
    act_sp <- 0
    for (s in sps) {
        # print(s)
        if (!s %in% c("Neophocaena_asiaeorientalis", "Rhinopithecus_roxellana")) {
            if (any(row[[paste(s, tis, sep = "_")]] == 1) & (any(row[[paste(s, dt[1], sep = "_")]] == 1) | any(row[[paste(s, dt[2], sep = "_")]] == 1))) {
                res$s1_1 <- res$s1_1 + 1
                act_sp <- act_sp + 1
            } else if (any(row[[paste(s, tis, sep = "_")]] == 1) & all(row[[paste(s, dt[1], sep = "_")]] == 0) & all(row[[paste(s, dt[2], sep = "_")]] == 0)) {
                res$s1_0 <- res$s1_0 + 1
                act_sp <- act_sp + 1
            } else if (all(row[[paste(s, tis, sep = "_")]] == 0) & (any(row[[paste(s, dt[1], sep = "_")]] == 1) | any(row[[paste(s, dt[2], sep = "_")]] == 1))) {
                res$s0_1 <- res$s0_1 + 1
                act_sp <- act_sp + 1
            } else {
                res$s0_0 <- res$s0_0 + 1
            }
        } else if ((s == "Neophocaena_asiaeorientalis" & tis == "Brain") | (s == "Rhinopithecus_roxellana" & tis != "Brain")) {
            if (any(row[[paste(s, tis, sep = "_")]] == 1)) {
                res$s1_0 <- res$s1_0 + 1
                act_sp <- act_sp + 1
            } else {
                res$s0_0 <- res$s0_0 + 1
            }
        }
        # else if ((s == "Neophocaena_asiaeorientalis" & tis != "Brain") | (s == "Rhinopithecus_roxellana" & tis == "Brain")) {
        #     if (row[[paste(s, tis, sep = "_")]] == 1) {
        #         res$s0_1 <- res$s1_1 + 1
        #     } else {
        #         res$s0_0 <- res$s0_0 + 1
        #     }
        # }
        if (s == "Neophocaena_asiaeorientalis") {
            tt <- c("Brain")
        } else if (s == "Rhinopithecus_roxellana") {
            tt <- c("Kidney", "Liver")
        } else {
            tt <- c("Brain", "Kidney", "Liver")
        }
        t_col <- paste(s, tt, sep = "_")
        if (any(row[t_col] == 1)) {
            act_time <- c(act_time, time_matrix[sp, s] / 2)
            print(s)
        }
    }
    res$act_time <- max(act_time)
    res$act_time <- ifelse(res$act_time == -Inf, 0, res$act_time)
    res$act_sp <- act_sp
    return(res)
}
library(doMC)
doMC::registerDoMC(cores = 6)
# plyr::a_ply(species, 1, function(sp) {
for (sp in species) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    ec_cols <- expand.grid(species = species[species != sp], tissue = c("Brain", "Kidney", "Liver"))
    ec_cols <- paste(ec_cols$species, ec_cols$tissue, sep = "_")
    ec_cols <- ec_cols[!ec_cols %in% c("Rhinopithecus_roxellana_Brain", "Neophocaena_asiaeorientalis_Kidney", "Neophocaena_asiaeorientalis_Liver")]
    tiss <- c("Brain", "Kidney", "Liver")
    if (any(c("s0_0", "s0_1", "s1_0", "s1_1", "act_time") %in% colnames(anno))) {
        anno <- dplyr::select(anno, -c("s0_0", "s0_1", "s1_0", "s1_1", "act_time", "act_sp"))
    }
    for (col in c("s0_0", "s0_1", "s1_0", "s1_1", "act_time", "act_sp")) {
        anno <- dplyr::mutate(anno, !!col := 0)
    }
    anno <- plyr::adply(anno, 1, count_ec, sps = species[species != sp], .parallel = TRUE)
    write.csv(anno, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), row.names = FALSE)
} # , .parallel = TRUE)

# ## a function to find the oldest active species of each element
# find_oldest <- function(row, sps) {
#     act_time <- c()
#     for (s in sps) {
#         if (s == "Neophocaena_asiaeorientalis") {
#             tt <- c("Brain")
#         } else if (s == "Rhinopithecus_roxellana") {
#             tt <- c("Kidney", "Liver")
#         } else {
#             tt <- c("Brain", "Kidney", "Liver")
#         }
#         t_col <- paste(s, tt, sep = "_")
#         if (any(row[t_col] == 1)) {
#             act_time <- c(act_time, time_matrix[sp, s] / 2)
#         }
#     }
#     return(max(act_time))
# }
summary(anno[which(anno$align > 0), c("s1_1", "s1_0", "s0_1", "s0_0", "ec_div_time", "act_time", "align_div_time")])
head(anno[which(anno$act_time == -Inf), c("s1_1", "s1_0", "s0_1", "s0_0", "act_time")])
head(anno[is.na(anno$ec_div_time), c("s1_1", "s1_0", "s0_1", "s0_0", "act_time", "ec_div_time")])
# change -Inf to 0
anno$act_time[which(anno$act_time == -Inf)] <- 0
library(plyr)
# a_ply(species, 1, function(sp) {
#     anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
#     anno$act_time <- plyr::adply(anno, 1, find_oldest, sps = species[species != sp], .parallel = TRUE)[, ""]
#     write.csv(act_time, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), row.names = FALSE)
# }, .parallel = TRUE)

## first, find the recently evolved elements, which means the elements that are not active in other species,
## in my data, there act_time will be 0. then split them into two groups, one is the elements in ancestral DNA
## sequences, the other is the elements in derived DNA sequences. the threshold is 100 ma for now.
## then, calculate the enrichment of rep elements in each group.
## i dont know how to define recently evolved elements in cell report, so firstly define as species-specific
species <- species[!species %in% c("Hipposideros_larvatus", "Rhinolophus_ferrumequinum", "Myotis_ricketti")]
ya_df <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    anno$act_time[which(anno$act_time == -Inf)] <- 0
    anno$ec_div_time[which(is.na(anno$ec_div_time))] <- 0
    write.csv(anno, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), row.names = FALSE)
    t3_ad_p_1 <- nrow(anno[which(anno$act_time == 0 & anno$align_div_time > 100 & anno$element == "promoter"), ]) # ancestral DNA, three tissue
    t3_yd_p_1 <- nrow(anno[which(anno$act_time == 0 & anno$align_div_time <= 40 & anno$element == "promoter"), ]) # young DNA, three tissue
    t3_ad_e_1 <- nrow(anno[which(anno$act_time == 0 & anno$align_div_time > 100 & anno$element == "enhancer"), ]) # ancestral DNA, three tissue
    t3_yd_e_1 <- nrow(anno[which(anno$act_time == 0 & anno$align_div_time <= 40 & anno$element == "enhancer"), ]) # young DNA, three tissue
    t1_ad_p_1 <- nrow(anno[which(anno$ec_div_time == 0 & anno$align_div_time > 100 & anno$element == "promoter"), ])
    t1_yd_p_1 <- nrow(anno[which(anno$ec_div_time == 0 & anno$align_div_time <= 40 & anno$element == "promoter"), ])
    t1_ad_e_1 <- nrow(anno[which(anno$ec_div_time == 0 & anno$align_div_time > 100 & anno$element == "enhancer"), ])
    t1_yd_e_1 <- nrow(anno[which(anno$ec_div_time == 0 & anno$align_div_time <= 40 & anno$element == "enhancer"), ])

    ## then try to use 40 ma as threshold
    t3_ad_p_2 <- nrow(anno[which(anno$act_time <= 40 & anno$align_div_time > 80 & anno$element == "promoter"), ]) # ancestral DNA, three tissue
    t3_yd_p_2 <- nrow(anno[which(anno$act_time <= 40 & anno$align_div_time <= 40 & anno$element == "promoter"), ]) # young DNA, three tissue
    t3_ad_e_2 <- nrow(anno[which(anno$act_time <= 40 & anno$align_div_time > 80 & anno$element == "enhancer"), ]) # ancestral DNA, three tissue
    t3_yd_e_2 <- nrow(anno[which(anno$act_time <= 40 & anno$align_div_time <= 40 & anno$element == "enhancer"), ]) # young DNA, three tissue
    t1_ad_p_2 <- nrow(anno[which(anno$ec_div_time <= 40 & anno$align_div_time > 80 & anno$element == "promoter"), ])
    t1_yd_p_2 <- nrow(anno[which(anno$ec_div_time <= 40 & anno$align_div_time <= 40 & anno$element == "promoter"), ])
    t1_ad_e_2 <- nrow(anno[which(anno$ec_div_time <= 40 & anno$align_div_time > 80 & anno$element == "enhancer"), ])
    t1_yd_e_2 <- nrow(anno[which(anno$ec_div_time <= 40 & anno$align_div_time <= 40 & anno$element == "enhancer"), ])

    ## recently evolved elements = act_time == 0, young DNA = align_div_time <= 40, ancestral DNA = align_div_time > 80
    t3_ad_p_3 <- nrow(anno[which(anno$act_time == 0 & anno$align_div_time > 80 & anno$element == "promoter"), ]) # ancestral DNA, three tissue
    t3_yd_p_3 <- nrow(anno[which(anno$act_time == 0 & anno$align_div_time <= 40 & anno$element == "promoter"), ]) # young DNA, three tissue
    t3_ad_e_3 <- nrow(anno[which(anno$act_time == 0 & anno$align_div_time > 80 & anno$element == "enhancer"), ]) # ancestral DNA, three tissue
    t3_yd_e_3 <- nrow(anno[which(anno$act_time == 0 & anno$align_div_time <= 40 & anno$element == "enhancer"), ]) # young DNA, three tissue
    t1_ad_p_3 <- nrow(anno[which(anno$ec_div_time == 0 & anno$align_div_time > 80 & anno$element == "promoter"), ])
    t1_yd_p_3 <- nrow(anno[which(anno$ec_div_time == 0 & anno$align_div_time <= 40 & anno$element == "promoter"), ])
    t1_ad_e_3 <- nrow(anno[which(anno$ec_div_time == 0 & anno$align_div_time > 80 & anno$element == "enhancer"), ])
    t1_yd_e_3 <- nrow(anno[which(anno$ec_div_time == 0 & anno$align_div_time <= 40 & anno$element == "enhancer"), ])

    res <- data.frame(
        species = sp, t3_ad_p_1 = t3_ad_p_1, t3_yd_p_1 = t3_yd_p_1, t3_ad_e_1 = t3_ad_e_1, t3_yd_e_1 = t3_yd_e_1,
        t1_ad_p_1 = t1_ad_p_1, t1_yd_p_1 = t1_yd_p_1, t1_ad_e_1 = t1_ad_e_1, t1_yd_e_1 = t1_yd_e_1,
        t3_ad_p_2 = t3_ad_p_2, t3_yd_p_2 = t3_yd_p_2, t3_ad_e_2 = t3_ad_e_2, t3_yd_e_2 = t3_yd_e_2,
        t1_ad_p_2 = t1_ad_p_2, t1_yd_p_2 = t1_yd_p_2, t1_ad_e_2 = t1_ad_e_2, t1_yd_e_2 = t1_yd_e_2,
        t3_ad_p_3 = t3_ad_p_3, t3_yd_p_3 = t3_yd_p_3, t3_ad_e_3 = t3_ad_e_3, t3_yd_e_3 = t3_yd_e_3,
        t1_ad_p_3 = t1_ad_p_3, t1_yd_p_3 = t1_yd_p_3, t1_ad_e_3 = t1_ad_e_3, t1_yd_e_3 = t1_yd_e_3
    )
    return(res)
}, .parallel = TRUE)
summary(ya_df)
ya_df[which(ya_df$species %in% c("Macaca_mulatta", "Mus_musculus", "Bos_taurus", "Canis_lupus")), ]
ya_df[which.max(ya_df$t1_yd_e_1), ]
ya_df <- ya_df[, -1]
write.csv(ya_df, "/media/Data/zhangz/chip/analysis/summary2/exap/young_anc_dna.csv", row.names = FALSE)

## draw box plot
## melt
ya_melt <- reshape2::melt(ya_df, id.vars = "species", variable.name = "t_ay_ep_m", value.name = "cre_num")
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(aplot)
ya_melt <- separate(ya_melt, t_ay_ep_m, into = c("tissue_num", "DNA_age", "element", "method"), sep = "_")
ya_melt <- ya_melt[ya_melt$species != "Petaurus_breviceps", ]
compars <- list(c("yd", "ad"))
ya_box_p_1t <- ggplot(
    ya_melt[which(ya_melt$element == "p" & ya_melt$method == "3" & ya_melt$tissue_num == "t1"), ],
    aes(x = DNA_age, y = cre_num / 1000, color = DNA_age)
) +
    geom_boxplot(fill = "white", outlier.size = 0.2) +
    geom_line(aes(group = species), color = "grey", linewidth = 0.1) +
    stat_compare_means(aes(x = DNA_age), comparisons = compars, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_color_manual(values = c(ad = "#02263e", yd = "#73b8d5")) +
    theme_light() +
    scale_x_discrete(limits = c("yd", "ad"), labels = c(yd = "Young DNA", ad = "Ancestral DNA")) +
    ylab("") +
    xlab("") +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none")
ya_box_p_3t <- ggplot(
    ya_melt[which(ya_melt$element == "p" & ya_melt$method == "3" & ya_melt$tissue_num == "t3"), ],
    aes(x = DNA_age, y = cre_num / 1000, color = DNA_age)
) +
    geom_boxplot(fill = "white", outlier.size = 0.2) +
    geom_line(aes(group = species), color = "grey", linewidth = 0.2) +
    stat_compare_means(aes(x = DNA_age), comparisons = compars, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_color_manual(values = c(ad = "#02263e", yd = "#73b8d5")) +
    theme_light() +
    scale_x_discrete(limits = c("yd", "ad"), labels = c(yd = "Young DNA", ad = "Ancestral DNA")) +
    ylab("") +
    xlab("") +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none")
ya_box_e_1t <- ggplot(
    ya_melt[which(ya_melt$element == "e" & ya_melt$method == "3" & ya_melt$tissue_num == "t1"), ],
    aes(x = DNA_age, y = cre_num, color = DNA_age)
) +
    geom_boxplot(fill = "white", outlier.size = 0.2) +
    geom_line(aes(group = species), color = "grey", linewidth = 0.3) +
    stat_compare_means(aes(x = DNA_age), comparisons = compars, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_color_manual(values = c(ad = "#c01d2e", yd = "#ef4343")) +
    theme_light() +
    scale_x_discrete(limits = c("yd", "ad"), labels = c(yd = "Young DNA", ad = "Ancestral DNA")) +
    ylab("") +
    xlab("") +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none")
ya_box_e_3t <- ggplot(
    ya_melt[which(ya_melt$element == "e" & ya_melt$method == "3" & ya_melt$tissue_num == "t3"), ],
    aes(x = DNA_age, y = cre_num, color = DNA_age)
) +
    geom_boxplot(fill = "white", outlier.size = 0.2) +
    geom_line(aes(group = species), color = "grey", linewidth = 0.4) +
    stat_compare_means(aes(x = DNA_age), comparisons = compars, method = "wilcox.test", label = "p.signif", paired = TRUE) +
    scale_color_manual(values = c(ad = "#c01d2e", yd = "#ef4343")) +
    theme_light() +
    scale_x_discrete(limits = c("yd", "ad"), labels = c(yd = "Young DNA", ad = "Ancestral DNA")) +
    ylab("") +
    xlab("") +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none")
library(cowplot)
ya_box <- plot_grid(ya_box_p_1t, ya_box_e_1t, ya_box_p_3t, ya_box_e_3t, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/young_anc_dna.pdf", ya_box, width = 20, height = 20, units = "cm")

## next, by dividing recently evolved elements into two groups, young DNA and ancestral DNA,
## calculate rep elements enrichment use rep number in cres as background


## first, get the rep number in cres
rep_enrich <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    anno <- anno[, c(colnames(anno)[1:44], "s0_0", "s0_1", "s1_0", "s1_1", "act_time", "act_sp")]
    rep_class <- paste(anno$repeat_class, collapse = ";") %>%
        strsplit(";") %>%
        unlist() %>%
        unique()
    rep_class <- rep_class[!rep_class %in% c("", "Low_complexity", "Simple_repeat", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", "srpRNA", "RNA", "DNA", "LTR", "LINE", "SINE", "Unknown")]
    rec_young_row_e <- which(anno$act_time == 0 & anno$align_div_time <= 40 & anno$element == "enhancer")
    rec_young_row_p <- which(anno$act_time == 0 & anno$align_div_time <= 40 & anno$element == "promoter")
    rec_anc_row_e <- which(anno$act_time == 0 & anno$align_div_time > 80 & anno$element == "enhancer")
    rec_anc_row_p <- which(anno$act_time == 0 & anno$align_div_time > 80 & anno$element == "promoter")
    cons_young_row_e <- which(anno$act_time > 0 & anno$align_div_time <= 40 & anno$element == "enhancer")
    cons_young_row_p <- which(anno$act_time > 0 & anno$align_div_time <= 40 & anno$element == "promoter")
    cons_anc_row_e <- which(anno$act_time > 0 & anno$align_div_time > 80 & anno$element == "enhancer")
    cons_anc_row_p <- which(anno$act_time > 0 & anno$align_div_time > 80 & anno$element == "promoter")
    e_row <- which(anno$element == "enhancer")
    p_row <- which(anno$element == "promoter")
    rec_young_e_num <- sum(anno$act_time == 0 & anno$align_div_time <= 40 & anno$element == "enhancer")
    rec_young_p_num <- sum(anno$act_time == 0 & anno$align_div_time <= 40 & anno$element == "promoter")
    rec_anc_e_num <- sum(anno$act_time == 0 & anno$align_div_time > 80 & anno$element == "enhancer")
    rec_anc_p_num <- sum(anno$act_time == 0 & anno$align_div_time > 80 & anno$element == "promoter")
    cons_young_e_num <- sum(anno$act_time > 0 & anno$align_div_time <= 40 & anno$element == "enhancer")
    cons_young_p_num <- sum(anno$act_time > 0 & anno$align_div_time <= 40 & anno$element == "promoter")
    cons_anc_e_num <- sum(anno$act_time > 0 & anno$align_div_time > 80 & anno$element == "enhancer")
    cons_anc_p_num <- sum(anno$act_time > 0 & anno$align_div_time > 80 & anno$element == "promoter")
    e_num <- sum(anno$element == "enhancer")
    p_num <- sum(anno$element == "promoter")
    rep_num <- adply(rep_class, 1, function(rep) {
        res <- data.frame(
            rep_class = rep,
            rec_young_e = sum(grepl(rep, anno$repeat_class[rec_young_row_e])),
            rec_young_p = sum(grepl(rep, anno$repeat_class[rec_young_row_p])),
            rec_anc_e = sum(grepl(rep, anno$repeat_class[rec_anc_row_e])),
            rec_anc_p = sum(grepl(rep, anno$repeat_class[rec_anc_row_p])),
            cons_young_e = sum(grepl(rep, anno$repeat_class[cons_young_row_e])),
            cons_young_p = sum(grepl(rep, anno$repeat_class[cons_young_row_p])),
            cons_anc_e = sum(grepl(rep, anno$repeat_class[cons_anc_row_e])),
            cons_anc_p = sum(grepl(rep, anno$repeat_class[cons_anc_row_p])),
            e_rep = sum(grepl(rep, anno$repeat_class[e_row])),
            p_rep = sum(grepl(rep, anno$repeat_class[p_row])),
            rec_young_e_num = rec_young_e_num,
            rec_young_p_num = rec_young_p_num,
            rec_anc_e_num = rec_anc_e_num,
            rec_anc_p_num = rec_anc_p_num,
            cons_young_e_num = cons_young_e_num,
            cons_young_p_num = cons_young_p_num,
            cons_anc_e_num = cons_anc_e_num,
            cons_anc_p_num = cons_anc_p_num,
            e_num = e_num,
            p_num = p_num
        )
        res$rec_young_e_p <- phyper(res$rec_young_e, res$e_rep, e_num - res$e_rep, rec_young_e_num, lower.tail = FALSE)
        res$rec_young_e_enrich <- (res$rec_young_e / res$rec_young_e_num) / (res$e_rep / res$e_num)
        res$rec_young_p_p <- phyper(res$rec_young_p, res$p_rep, p_num - res$p_rep, rec_young_p_num, lower.tail = FALSE)
        res$rec_young_p_enrich <- (res$rec_young_p / res$rec_young_p_num) / (res$p_rep / res$p_num)
        res$rec_anc_e_p <- phyper(res$rec_anc_e, res$e_rep, e_num - res$e_rep, rec_anc_e_num, lower.tail = FALSE)
        res$rec_anc_e_enrich <- (res$rec_anc_e / res$rec_anc_e_num) / (res$e_rep / res$e_num)
        res$rec_anc_p_p <- phyper(res$rec_anc_p, res$p_rep, p_num - res$p_rep, rec_anc_p_num, lower.tail = FALSE)
        res$rec_anc_p_enrich <- (res$rec_anc_p / res$rec_anc_p_num) / (res$p_rep / res$p_num)
        res$cons_young_e_p <- phyper(res$cons_young_e, res$e_rep, e_num - res$e_rep, cons_young_e_num, lower.tail = FALSE)
        res$cons_young_e_enrich <- (res$cons_young_e / res$cons_young_e_num) / (res$e_rep / res$e_num)
        res$cons_young_p_p <- phyper(res$cons_young_p, res$p_rep, p_num - res$p_rep, cons_young_p_num, lower.tail = FALSE)
        res$cons_young_p_enrich <- (res$cons_young_p / res$cons_young_p_num) / (res$p_rep / res$p_num)
        res$cons_anc_e_p <- phyper(res$cons_anc_e, res$e_rep, e_num - res$e_rep, cons_anc_e_num, lower.tail = FALSE)
        res$cons_anc_e_enrich <- (res$cons_anc_e / res$cons_anc_e_num) / (res$e_rep / res$e_num)
        res$cons_anc_p_p <- phyper(res$cons_anc_p, res$p_rep, p_num - res$p_rep, cons_anc_p_num, lower.tail = FALSE)
        res$cons_anc_p_enrich <- (res$cons_anc_p / res$cons_anc_p_num) / (res$p_rep / res$p_num)
        return(res)
    }, .parallel = TRUE)
    rep_num <- rep_num[, -1]
    rep_num$species <- sp
    return(rep_num)
}, .parallel = TRUE)
summary(rep_enrich)
rep_enrich <- rep_enrich[, -1]
rep_enrich <- rep_enrich[, c(1, 38, 2:37)]
write.csv(rep_enrich, "/media/Data/zhangz/chip/analysis/summary2/exap/rep_enrich.csv", row.names = FALSE)

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
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("ancestral DNA")
rec_e_rep_p <- rec_young_e_rep_p %>% insert_right(rec_anc_e_rep_p, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/rec_e_rep_p.pdf", rec_e_rep_p, height = 12, width = 18, unit = "cm")

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
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("ancestral DNA")
rec_p_rep_p <- rec_young_p_rep_p %>% insert_right(rec_anc_p_rep_p, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/rec_p_rep_p.pdf", rec_p_rep_p, height = 12, width = 18, unit = "cm")

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
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("ancestral DNA")
cons_e_rep_p <- cons_young_e_rep_p %>% insert_right(cons_anc_e_rep_p, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/cons_e_rep_p.pdf", cons_e_rep_p, height = 12, width = 18, unit = "cm")

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
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("ancestral DNA")
cons_p_rep_p <- cons_young_p_rep_p %>% insert_right(cons_anc_p_rep_p, width = 1)
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/cons_p_rep_p.pdf", cons_p_rep_p, height = 12, width = 18, unit = "cm")

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
pseudo_ts[which(pseudo_ts$pseudo_ratio_e == 0), ]
## boxplot of pseudo ratio in enhancer and promoter
pseudo_ts_melt <- reshape2::melt(pseudo_ts[, c("species", "pseudo_ratio_e", "pseudo_ratio_p")], id.vars = "species", variable.name = "element", value.name = "pseudo_ratio")
pseudo_ts_melt$element <- pseudo_ts_melt$element %>%
    gsub("pseudo_ratio_e", "enhancer", .) %>%
    gsub("pseudo_ratio_p", "promoter", .)
pseudo_ts_box <- ggplot(pseudo_ts_melt, aes(x = element, y = pseudo_ratio, color = element)) +
    geom_boxplot(fill = "white", outlier.size = 0.2) +
    # geom_line(aes(group = species), color = "grey", linewidth = 0.2) +
    scale_color_manual(values = c("enhancer" = "#c01d2e", "promoter" = "#033250")) +
    theme_light() +
    ylab("Evolutionary pleiotropy ratio") +
    xlab("") +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none")
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/pseudo_ts_box.pdf", pseudo_ts_box, height = 6, width = 6, unit = "cm")

# how many cres were considered as exaptation by 1 tissue but pleiotropy by 3 tissues
## the main aim of below if to calculate the ratio of evolutionary pleiotropy in enhancer and promoter
evol_ts <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    evol_pleio_e <- nrow(anno[which(!anno$ele_pattern %in% c("B", "K", "L") & (anno$s0_1 != 0 | anno$s1_1 != 0) & anno$element == "enhancer"), ])
    pleio_e_sp <- nrow(anno[which(!anno$ele_pattern %in% c("B", "K", "L") & anno$element == "enhancer"), ])
    evol_pleio_p <- nrow(anno[which(!anno$ele_pattern %in% c("B", "K", "L") & (anno$s0_1 != 0 | anno$s1_1 != 0) & anno$element == "promoter"), ])
    pleio_p_sp <- nrow(anno[which(!anno$ele_pattern %in% c("B", "K", "L") & anno$element == "promoter"), ])
    # exaptation is recently evolved elements from ancestral DNA,
    pseudo_exapt_e <- nrow(anno[which(anno$act_time > 0 & anno$align_div_time > 80 & anno$ec_div_time == 0 & anno$element == "enhancer"), ])
    pseudo_exapt_p <- nrow(anno[which(anno$act_time > 0 & anno$align_div_time > 80 & anno$ec_div_time == 0 & anno$element == "promoter"), ])
    # pleio_e <- nrow(anno[which(anno$pleio_sp > 0 & anno$element == "enhancer"), ])
    evol_ratio_e <- evol_pleio_e / pleio_e_sp
    evol_ratio_p <- evol_pleio_p / pleio_p_sp
    pseudo_exapt_e_ratio <- pseudo_exapt_e / nrow(anno[which(anno$element == "enhancer" & anno$align_div_time > 80), ])
    pseudo_exapt_p_ratio <- pseudo_exapt_p / nrow(anno[which(anno$element == "promoter" & anno$align_div_time > 80), ])
    pseudo_sp_e <- nrow(anno[which(anno$act_time > 0 & anno$ec_div_time == 0 & anno$element == "enhancer"), ])
    pseudo_sp_p <- nrow(anno[which(anno$act_time > 0 & anno$ec_div_time == 0 & anno$element == "promoter"), ])
    pseudo_sp_e_ratio <- pseudo_sp_e / nrow(anno[which(anno$element == "enhancer" & anno$ec_div_time == 0), ])
    pseudo_sp_p_ratio <- pseudo_sp_p / nrow(anno[which(anno$element == "promoter" & anno$ec_div_time == 0), ])
    pleio_e <- nrow(anno[which(anno$pleio_sp > 0 & anno$element == "enhancer"), ])
    pleio_p <- nrow(anno[which(anno$pleio_sp > 0 & anno$element == "promoter"), ])
    pleio_e_ratio <- pleio_e / nrow(anno[which(anno$element == "enhancer"), ])
    pleio_p_ratio <- pleio_p / nrow(anno[which(anno$element == "promoter"), ])
    pleio_e_sp_evo <- nrow(anno[which(anno$pleio_sp > 0 & !anno$ele_pattern %in% c("B", "K", "L") & anno$element == "enhancer"), ])
    pleio_p_sp_evo <- nrow(anno[which(anno$pleio_sp > 0 & !anno$ele_pattern %in% c("B", "K", "L") & anno$element == "promoter"), ])
    pleio_e_sp_evo_ratio <- pleio_e_sp_evo / nrow(anno[which(!anno$ele_pattern %in% c("B", "K", "L") & anno$element == "enhancer"), ])
    pleio_p_sp_evo_ratio <- pleio_p_sp_evo / nrow(anno[which(!anno$ele_pattern %in% c("B", "K", "L") & anno$element == "promoter"), ])
    res <- data.frame(
        species = sp, evol_pleio_e = evol_pleio_e, pleio_e_sp = pleio_e_sp, evol_ratio_e = evol_ratio_e,
        evol_pleio_p = evol_pleio_p, pleio_p_sp = pleio_p_sp, evol_ratio_p = evol_ratio_p,
        pseudo_exapt_e = pseudo_exapt_e, pseudo_exapt_p = pseudo_exapt_p,
        pseudo_exapt_e_ratio = pseudo_exapt_e_ratio, pseudo_exapt_p_ratio = pseudo_exapt_p_ratio,
        pseudo_sp_e = pseudo_sp_e, pseudo_sp_p = pseudo_sp_p,
        pseudo_sp_e_ratio = pseudo_sp_e_ratio, pseudo_sp_p_ratio = pseudo_sp_p_ratio,
        pleio_e = pleio_e, pleio_p = pleio_p, pleio_e_ratio = pleio_e_ratio, pleio_p_ratio = pleio_p_ratio,
        pleio_e_sp_evo = pleio_e_sp_evo, pleio_p_sp_evo = pleio_p_sp_evo,
        pleio_e_sp_evo_ratio = pleio_e_sp_evo_ratio, pleio_p_sp_evo_ratio = pleio_p_sp_evo_ratio
    )
    return(res)
}, .parallel = TRUE)
evol_ts <- evol_ts[, -1]
summary(evol_ts)
## boxplot of evolutionary pleiotropy ratio in enhancer and promoter
evol_ts_melt <- reshape2::melt(evol_ts[, c("species", "evol_ratio_e", "evol_ratio_p")], id.vars = "species", variable.name = "element", value.name = "evol_ratio")
evol_ts_melt$element <- evol_ts_melt$element %>%
    gsub("evol_ratio_e", "enhancer", .) %>%
    gsub("evol_ratio_p", "promoter", .)
evol_ts_box <- ggplot(evol_ts_melt, aes(x = element, y = evol_ratio, color = element)) +
    geom_boxplot(fill = "white", outlier.size = 0.2) +
    # geom_line(aes(group = species), color = "grey", linewidth = 0.2) +
    scale_color_manual(values = c("enhancer" = "#c01d2e", "promoter" = "#033250")) +
    theme_light() +
    ylab("Evolutionary pleiotropy ratio") +
    xlab("") +
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none")
ggsave("/media/Data/zhangz/chip/analysis/summary2/exap/evol_ts_box.pdf", evol_ts_box, height = 6, width = 6, unit = "cm")

pleio_judge <- function(row, sps) {
    cols <- expand.grid(species = sps, tissue = c("Brain", "Kidney", "Liver"))
    cols <- paste0(cols$species, "_", cols$tissue)
    cols <- cols[cols %in% colnames(row)]
    pleio_sp <- 0
    pleio_time <- c()
    for (ss in sps) {
        if (sum(row[, cols[grepl(ss, cols)]]) > 1) {
            pleio_sp <- pleio_sp + 1
            pleio_time <- c(pleio_time, time_matrix[sp, ss] / 2)
        }
    }
    pleio_time <- max(pleio_time)
    return(data.frame(pleio_sp = pleio_sp, pleio_time = pleio_time))
}
## how many pleiotropy elements were also pleiotropy in at least one other species
pleio_ts <- adply(species, 1, function(sp) {
    anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    sps <- setdiff(species, sp)
    anno <- adply(anno, 1, pleio_judge, sps, .parallel = TRUE)
    write.csv(anno, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), row.names = FALSE)
    pleio_e <- nrow(anno[which(anno$pleio_sp > 0 & anno$element == "enhancer"), ])
    pleio_p <- nrow(anno[which(anno$pleio_sp > 0 & anno$element == "promoter"), ])
    pleio_e_ratio <- pleio_e / nrow(anno[which(anno$element == "enhancer"), ])
    pleio_p_ratio <- pleio_p / nrow(anno[which(anno$element == "promoter"), ])
    res <- data.frame(
        species = sp, pleio_e = pleio_e, pleio_p = pleio_p, pleio_e_ratio = pleio_e_ratio, pleio_p_ratio = pleio_p_ratio
    )
    return(res)
}, .parallel = TRUE)
pleio_ts <- pleio_ts[, -1]
summary(pleio_ts)
