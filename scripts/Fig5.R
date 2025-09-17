# draw plots for Fig5, start from raw data process
# 2025-06-29
# load libraries
library(magrittr)
library(ggplot2)
library(Biostrings)
library(reshape2)
library(doMC)
library(plyr)
library(dplyr)
library(tidyr)
library(ggpubr)
# library(ggbreak)
registerDoMC(cores = 4)

# --------------------------------------A--------------------------------------------------------------
## panel A, basic stats of cres involved in mcl and in top groups

# read cre ids used for mcl
# cre_ids <- read.csv("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl_input.tab",
#     header = FALSE, sep = '\t', row.names = 1
# )
# colnames(cre_ids) <- 'id'
# cre_ids$species <- sapply(strsplit(as.character(cre_ids$id), "_"), function(x) paste(x[1:2], collapse = '_'))
# cre_ids <- aaply(cre_ids, 1, function(x) {
#     eles <- strsplit(as.character(x["id"]), "_")[[1]]
#     x["species"] <- paste(eles[1:2], collapse = "_")
#     x["cre"] <- ifelse(eles[3] == 'pro', 'promoter', 'enhancer')
#     return(x)
# }, .parallel = TRUE)
# use awk command to do this
# awk -F'\t' '{split($2, a, "_");print a[1]"_"a[2] "\t" a[3]}' all_cre_blast_para_mcl_input.tab > all_cre_blast_para_mcl_input_table.txt
cre_ids <- read.csv("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl_input_table.txt",
    header = FALSE, sep = "\t"
)
colnames(cre_ids) <- c("species", "cre")
# count how many enh and pro in each species
cre_ids <- cre_ids %>%
    dplyr::group_by(species, cre) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>% # 添加.groups参数避免警告
    # pivot_wider(
    #     names_from = cre,
    #     values_from = count,
    #     values_fill = list(count = 0) # 处理缺失值
    # ) %>%
    arrange(species) # 按物种排序
top_cre_ids <- read.csv("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/top_cres_table.tsv",
    header = FALSE, sep = "\t"
)
species_order <- c(
    "Ovis_aries", "Bos_taurus", "Neophocaena_asiaeorientalis",
    "Sus_scrofa", "Lama_glama", "Mustela_putorius",
    "Canis_lupus", "Felis_catus", "Equus_asinus",
    "Equus_caballus", "Myotis_chinensis", "Rhinolophus_pusillus",
    "Atelerix_albiventris", "Mus_musculus", "Rattus_norvegicus",
    "Cavia_porcellus", "Oryctolagus_cuniculus", "Rhinopithecus_roxellana",
    "Macaca_mulatta", "Tupaia_belangeri", "Procavia_capensis", "Petaurus_breviceps"
)
colnames(top_cre_ids) <- c("species", "cre")
top_cre_ids <- top_cre_ids %>%
    dplyr::group_by(species, cre) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>% # 添加.groups参数避免警告
    # pivot_wider(
    #     names_from = cre,
    #     values_from = count,
    #     values_fill = list(count = 0) # 处理缺失值
    # ) %>%
    arrange(species) # 按物种排序
# concat cre_ids and top_cre_ids
cre_ids$group <- "all"
top_cre_ids$group <- "top"
cre_ids <- rbind(cre_ids, top_cre_ids)
cre_ids$count[cre_ids$cre == "pro"] <- -cre_ids$count[cre_ids$cre == "pro"]
# two-way (enh and pro) col plot of how many cres involved in mcl and in top groups
num_col_p <- ggplot(cre_ids, aes(x = count, y = species, fill = cre)) +
    geom_bar(stat = "identity", position = "identity", width = 0.8, alpha = 0.7) +
    # geom_bar(data = cre_ids[cre_ids$group == "top", ], stat = "identity", position = "identity", width = 0.5, alpha = 0.5) +
    scale_fill_manual(values = c("enh" = "#ef4343", "pro" = "#73b8d5")) +
    labs(x = "Number of CRES", y = "Species", fill = "CRE Type") +
    scale_x_continuous(labels = function(x) abs(x)) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10, face = "italic"),
    ) +
    facet_wrap(~group, scales = "free_x") # add facet wrap to separate all and top groups
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/fig5/top_cre.pdf", num_col_p, height = 4, width = 5)



# 从/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl_input.tab取前2-10行，
# 每一列保存为一个文件，其中每一行是制表符分割的，将每个元素转换成一行
# bash /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/split_count.sh

all_num_p <- ggplot(cre_ids[cre_ids$group == "all", ], aes(x = count, y = species, fill = cre)) +
    geom_bar(stat = "identity", position = "identity", width = 0.8, alpha = 0.7) +
    scale_fill_manual(values = c("enh" = "#ef4343", "pro" = "#73b8d5")) +
    labs(x = "Number of CRES", y = "Species", fill = "CRE Type") +
    scale_x_continuous(labels = function(x) abs(x)) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x), limits = rev(species_order)) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10, face = "italic"),
        axis.line = element_line(linewidth = 0.5, color = "grey90"),
    )
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/fig5/all_cre.pdf", all_num_p, height = 8, width = 6)


grouped_ids <- adply(1:20, 1, function(i) {
    # read No.i group and count ele number
    group_ids <- read.csv(
        paste0("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/", i, "_temp.tabel"),
        header = FALSE, sep = "\t"
    )
    colnames(group_ids) <- c("species", "cre")
    group_ids <- group_ids %>%
        dplyr::group_by(species, cre) %>%
        dplyr::summarise(count = n(), .groups = "drop") %>% # 添加.groups参数避免警告
        arrange(species) # 按物种排序
    group_ids$group <- i
    group_ids$count[group_ids$cre == "pro"] <- -group_ids$count[group_ids$cre == "pro"]
    write.csv(group_ids,
        paste0("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/group_", i, "_cres_count.csv"),
        row.names = FALSE
    )
    return(group_ids)
}, .parallel = TRUE)
grouped_ids <- grouped_ids[, -1]
grouped_ids <- adply(1:20, 1, function(i) {
    # read No.i group and count ele number
    group_ids <- read.csv(
        paste0("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/group_", i, "_cres_count.csv"),
        header = TRUE, sep = ","
    )
    group_ids$group <- i
    return(group_ids)
}, .parallel = TRUE)
grouped_ids <- grouped_ids[, -1]

label_fun <- function(variable, value) {
    # create a label for facet_wrap
    return(paste0("Group ", value))
}

cluster_source <- ggplot(grouped_ids, aes(x = count, y = species, fill = cre)) +
    geom_bar(stat = "identity", position = "identity", width = 0.8, alpha = 0.7) +
    scale_fill_manual(values = c("enh" = "#ef4343", "pro" = "#73b8d5")) +
    labs(x = "Number of CRES", y = "Species", fill = "CRE Type") +
    scale_x_continuous(labels = function(x) abs(x), n.breaks = 3) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x), limits = rev(species_order)) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10, face = "italic"),
        axis.line = element_line(linewidth = 0.5, color = "grey90"),
    ) +
    facet_wrap(~group, scales = "free_x", nrow = 2, labeller = label_fun)
ggsave("/media/Data/zhangz/chip/analysis/summary2/fig/fig5/mcl_cluster_source.pdf", cluster_source, height = 8, width = 10)


library(data.table)

mcl_identity <- fread(
    "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl_input.txt",
    select = 3
)
# draw a distribution density of mcl_identity
mcl_identity_p <- ggplot(mcl_identity, aes(x = V3)) +
    geom_density(fill = "#ef4343", alpha = 0.7) +
    labs(x = "Blast identity (%)", y = "Density") +
    scale_x_continuous(limits = c(50, 100), n.breaks = 3) +
    scale_y_continuous(n.breaks = 3) +
    theme_classic() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
ggsave("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/mcl_identity_distribution.pdf", mcl_identity_p, width = 4, height = 3, dpi = 300)



mcl_identity_all <- fread(
    cmd = "awk 'BEGIN {srand(); rate=0.01} rand() < rate' /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/all_cre_blast.out",
    select = 3
)
mcl_identity_all_p <- ggplot(mcl_identity_all, aes(x = V3)) +
    geom_density(fill = "#ef4343", alpha = 0.7) +
    labs(x = "Identity (%)", y = "Density") +
    scale_x_continuous(limits = c(0, 100), n.breaks = 5) +
    scale_y_continuous(n.breaks = 5) +
    theme_classic() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
ggsave("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/mcl_identity_distribution_all.pdf", mcl_identity_all_p, width = 8, height = 6, dpi = 300)


# --------------------------------------B--------------------------------------------------------------


## read the mcl result by lines
# mcl_res <- "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_blast_mcl.mcl"
## only for test, use results of single species
mcl_res <- "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl.mcl"
mcl_res <- readLines(mcl_res)
## seperate each line by \t and count how many elements in each line
mcl_res <- lapply(mcl_res, function(x) {
    return(strsplit(x, "\t")[[1]])
})
mcl_count <- lapply(mcl_res, function(x) {
    return(length(x))
}) %>% unlist()
sum(mcl_count)

## accumulation of mcl_count
mcl_count_acc <- cumsum(mcl_count)
mcl_count_acc <- data.frame(rank = seq_along(mcl_count), accumulate = mcl_count_acc)
mcl_count_acc <- rbind(data.frame(rank = 0, accumulate = 0), mcl_count_acc)
## how many clusters account for 80% of all cres
ye <- sum(mcl_count) * 0.8
cluster_num <- mcl_count_acc[which(mcl_count_acc$accumulate >= ye), "rank"][1]
mcl_accumu_p <- ggplot(mcl_count_acc, aes(x = rank, y = accumulate)) +
    # geom_line() +
    geom_smooth(method = "gam") +
    # geom_point() +
    ## add line of 80% of total number of cres
    annotate("segment", x = 0, xend = cluster_num, y = ye, yend = ye, color = "#ef4343") +
    annotate("text", x = 3, y = ye + 200, label = paste0("80% of total number of cres: ", ye), hjust = 0, vjust = 1, color = "red", size = 8) +
    annotate("segment", x = cluster_num, xend = cluster_num, y = 0, yend = ye, color = "#ef4343") +
    annotate("text", x = cluster_num + 4, y = 80, label = paste0("cluster number: ", cluster_num), hjust = 0, vjust = 1, color = "red", size = 8) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 600000, 1200000)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 6042856), n.breaks = 3) +
    theme_classic() +
    labs(x = "Cluster rank", y = "Accumulated number of CREs") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
ggsave("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/mcl_cluster_accumulate.pdf", mcl_accumu_p, width = 4, height = 3, dpi = 300)

# -----------------------------------------------------------------------------C--------------------------------------------------------------
## result of the second dual-luciferase reporter assay
# read the result of dual-luciferase reporter assay
dual_luc <- read.csv("/media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/res/dual_luci_res.csv", header = TRUE)
luc_plot <- ggplot(dual_luc, aes(x = group_c, y = Standardized_F.R, fill = group_c)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.4) +
    geom_violin(aes(fill = group_c), alpha = 0.3, scale = "width", adjust = 1.5, width = 0.6) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    labs(y = "Relative enhancer activity", x = "") +
    scale_fill_manual(values = c("#ef4343", "grey10", "grey70")) +
    stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(c("exp", "neg"))) +
    theme_classic() +
    scale_x_discrete(
        limits = c("pos", "exp", "neg"),
        labels = c("Positive control", "De novo CREs", "Negative control")
    ) +
    theme(
        axis.text.x = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
    )
ggsave("/media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/res/dual_luci_res.pdf", luc_plot, width = 6, height = 4, dpi = 600)

random_overlap <- read.csv("/media/Data/zhangz/chip/analysis/summary2/blast_all/pattern_search/bos_random/random_overlap.txt", header = TRUE)
#       X42
#  Min.   : 0.00
#  1st Qu.:37.00
#  Median :37.00
#  Mean   :34.78
#  3rd Qu.:37.00
#  Max.   :53.00
colnames(random_overlap) <- c("ratio")
random_overlap$ratio <- random_overlap$ratio / 1000
random_overlap$type <- "random"
bos_ratio <- data.frame(ratio = c(180 / 2281, 96 / 1242), type = "bos")
# [1] 0.07891276 0.07729469
mus_ratio <- data.frame(ratio = c(6 / 278, 6 / 200), type = "mus")
# [1] 0.02158273 0.03000000
overlap_ratio <- rbind(random_overlap, bos_ratio)
overlap_p <- ggplot(overlap_ratio, aes(x = type, y = ratio)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.4) +
    geom_violin(aes(fill = type), alpha = 0.3, scale = "width", adjust = 1.5, width = 0.6) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    labs(y = "Overlap ratio", x = "") +
    stat_compare_means(
        method = "wilcox",
        label = "p.signif",
        comparisons = list(c("random", "bos"))
    ) +
    scale_y_continuous(n.breaks = 3) +
    theme_classic() +
    scale_x_discrete(labels = c("Pattern", "Random")) +
    theme(legend.position = 'none')
ggsave("/media/Data/zhangz/chip/analysis/summary2/blast_all/pattern_search/bos_random/random_overlap.pdf", overlap_p, width = 4, height = 3, dpi = 600)

mus_random <- read.csv("/media/Data/zhangz/chip/analysis/summary2/blast_all/pattern_search/mus_random/random_overlap.txt", header = TRUE)
colnames(mus_random) <- c("ratio")
mus_random$ratio <- mus_random$ratio / 1000
mus_random$type <- "random"
mus_ratio <- data.frame(ratio = c(6 / 278, 6 / 200), type = "mus")
mus_ratio <- rbind(mus_ratio, mus_random)
mus_overlap_p <- ggplot(mus_ratio, aes(x = type, y = ratio)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.4) +
    geom_violin(aes(fill = type), alpha = 0.3, scale = "width", adjust = 1.5, width = 0.6) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    labs(y = "Overlap ratio", x = "") +
    stat_compare_means(
        method = "wilcox",
        label = "p.signif",
        comparisons = list(c("random", "mus"))
    )
ggsave("/media/Data/zhangz/chip/analysis/summary2/blast_all/pattern_search/mus_random/random_overlap.pdf", mus_overlap_p, width = 6, height = 4, dpi = 600)
