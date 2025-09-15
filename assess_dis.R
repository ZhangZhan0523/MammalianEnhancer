## assess the influence of different distance threshold on the annotation
library(ggplot2)
library(reshape2)
library(stringr)
library(magrittr)

# sp <- "Mus_musculus"
diss <- c("1M", "500k", "200k", "100k", "50k")
sp_all <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_all_1M.csv"), header = TRUE)
d_cnames <- list(
    "1M" = "anno_gene_1M", "500k" = "anno_gene_500k", "200k" = "anno_gene_200k",
    "100k" = "anno_gene_100k", "50k" = "anno_gene_50k"
)
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
species_nots <- c(
    "Neophocaena_asiaeorientalis", "Rhinopithecus_roxellana",
    "Rhinolophus_ferrumequinum", "Myotis_ricketti"
)
species_ts <- species[!species %in% species_nots]

setwd("/media/Data/zhangz/chip/analysis/summary2/distance_assess")
# first is the number of annotated genes to the each cre and as a whole
gene_dis <- data.frame()
for (sp in species_ts) {
    sp_all <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_all_1M.csv"), header = TRUE)
    tmp_df <- data.frame(species = sp)
    for (dis in diss) {
        cname <- d_cnames[[dis]]
        for (ele in c("enhancer", "promoter")) {
            genes <- unique(str_split(paste(sp_all[[cname]][sp_all$element == ele], collapse = ";"), ";")[[1]])
            genes <- genes[genes != "NULL" & genes != ""]
            tmp_df[paste(ele, dis, sep = "_")] <- length(genes)
        }
    }
    # remove "NULL", "", NA
    sp_all$g_per_e_1M <- sapply(sp_all$anno_gene_1M, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    sp_all$g_per_e_500k <- sapply(sp_all$anno_gene_500k, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    sp_all$g_per_e_200k <- sapply(sp_all$anno_gene_200k, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    sp_all$g_per_e_100k <- sapply(sp_all$anno_gene_100k, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    sp_all$g_per_e_50k <- sapply(sp_all$anno_gene_50k, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    tmp_df$g_per_e_50k_mean <- mean(sp_all$g_per_e_50k[sp_all$g_per_e_50k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_50k_median <- median(sp_all$g_per_e_50k[sp_all$g_per_e_50k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_50k_max <- max(sp_all$g_per_e_50k[sp_all$g_per_e_50k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_100k_mean <- mean(sp_all$g_per_e_100k[sp_all$g_per_e_100k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_100k_median <- median(sp_all$g_per_e_100k[sp_all$g_per_e_100k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_100k_max <- max(sp_all$g_per_e_100k[sp_all$g_per_e_100k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_200k_mean <- mean(sp_all$g_per_e_200k[sp_all$g_per_e_200k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_200k_median <- median(sp_all$g_per_e_200k[sp_all$g_per_e_200k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_200k_max <- max(sp_all$g_per_e_200k[sp_all$g_per_e_200k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_500k_mean <- mean(sp_all$g_per_e_500k[sp_all$g_per_e_500k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_500k_median <- median(sp_all$g_per_e_500k[sp_all$g_per_e_500k != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_500k_max <- max(sp_all$g_per_e_500k[sp_all$g_per_e_500k != 0 & sp_all$element == "enhancer"])

    # tmp_df$g_per_e_1M_min <- min(sp_all$g_per_e_1M)
    tmp_df$g_per_e_1M_mean <- mean(sp_all$g_per_e_1M[sp_all$g_per_e_1M != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_1M_median <- median(sp_all$g_per_e_1M[sp_all$g_per_e_1M != 0 & sp_all$element == "enhancer"])
    tmp_df$g_per_e_1M_max <- max(sp_all$g_per_e_1M[sp_all$g_per_e_1M != 0 & sp_all$element == "enhancer"])
    # tmp_df$g_per_e_500k_min <- min(sp_all$g_per_e_500k)
    # tmp_df$g_per_e_200k_min <- min(sp_all$g_per_e_200k)
    tmp_df$g_per_p_50k_mean <- mean(sp_all$g_per_e_50k[sp_all$g_per_e_50k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_50k_median <- median(sp_all$g_per_e_50k[sp_all$g_per_e_50k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_50k_max <- max(sp_all$g_per_e_50k[sp_all$g_per_e_50k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_100k_mean <- mean(sp_all$g_per_e_100k[sp_all$g_per_e_100k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_100k_median <- median(sp_all$g_per_e_100k[sp_all$g_per_e_100k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_100k_max <- max(sp_all$g_per_e_100k[sp_all$g_per_e_100k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_200k_mean <- mean(sp_all$g_per_e_200k[sp_all$g_per_e_200k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_200k_median <- median(sp_all$g_per_e_200k[sp_all$g_per_e_200k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_200k_max <- max(sp_all$g_per_e_200k[sp_all$g_per_e_200k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_500k_mean <- mean(sp_all$g_per_e_500k[sp_all$g_per_e_500k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_500k_median <- median(sp_all$g_per_e_500k[sp_all$g_per_e_500k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_500k_max <- max(sp_all$g_per_e_500k[sp_all$g_per_e_500k != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_1M_mean <- mean(sp_all$g_per_e_1M[sp_all$g_per_e_1M != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_1M_median <- median(sp_all$g_per_e_1M[sp_all$g_per_e_1M != 0 & sp_all$element == "promoter"])
    tmp_df$g_per_p_1M_max <- max(sp_all$g_per_e_1M[sp_all$g_per_e_1M != 0 & sp_all$element == "promoter"])
    gene_dis <- rbind(gene_dis, tmp_df)
}
gene_dis_melt <- melt(gene_dis, id.vars = "species")

e_num_p <- ggplot(gene_dis_melt[gene_dis_melt$variable %in% c("enhancer_1M", "enhancer_500k", "enhancer_200k", "enhancer_100k", "enhancer_50k"), ], aes(x = variable, y = value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Number of annotated genes in enhancer", x = "distance", y = "Number of genes") +
    scale_x_discrete(limits = c("enhancer_50k", "enhancer_100k", "enhancer_200k", "enhancer_500k", "enhancer_1M"),
    labels = c("enhancer_50k" = "50k", "enhancer_100k" = "100k", "enhancer_200k" = "200k", "enhancer_500k" = "500k", "enhancer_1M" = "1M")) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(plot.title = element_text(size = 30))
ggsave("e_num_p.png", e_num_p, width = 15, height = 10)
p_num_p <- ggplot(gene_dis_melt[gene_dis_melt$variable %in% c("promoter_1M", "promoter_500k", "promoter_200k", "promoter_100k", "promoter_50k"), ], aes(x = variable, y = value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Number of annotated genes in promoter", x = "Species", y = "Number of genes") +
    scale_x_discrete(
        limits = c("promoter_50k", "promoter_100k", "promoter_200k", "promoter_500k", "promoter_1M"),
        labels = c("promoter_50k" = "50k", "promoter_100k" = "100k", "promoter_200k" = "200k", "promoter_500k" = "500k", "promoter_1M" = "1M")
    ) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(plot.title = element_text(size = 30))
ggsave("p_num_p.png", p_num_p, width = 15, height = 10)
g_per_e_mean_p <- ggplot(gene_dis_melt[gene_dis_melt$variable %in% c("g_per_e_1M_mean", "g_per_e_500k_mean", "g_per_e_200k_mean", "g_per_e_100k_mean", "g_per_e_50k_mean"), ], aes(x = variable, y = value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mean number of annotated genes per enhancer", x = "Species", y = "Number of genes") +
    scale_x_discrete(
        limits = c("g_per_e_50k_mean", "g_per_e_100k_mean", "g_per_e_200k_mean", "g_per_e_500k_mean", "g_per_e_1M_mean"),
        labels = c("g_per_e_50k_mean" = "50k", "g_per_e_100k_mean" = "100k", "g_per_e_200k_mean" = "200k", "g_per_e_500k_mean" = "500k", "g_per_e_1M_mean" = "1M")
    ) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(plot.title = element_text(size = 30))
ggsave("g_per_e_mean_p.png", g_per_e_mean_p, width = 15, height = 10)
g_per_e_median_p <- ggplot(gene_dis_melt[gene_dis_melt$variable %in% c("g_per_e_1M_median", "g_per_e_500k_median", "g_per_e_200k_median", "g_per_e_100k_median", "g_per_e_50k_median"), ], aes(x = variable, y = value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Median number of annotated genes per enhancer", x = "Species", y = "Number of genes") +
    scale_x_discrete(limits = c("g_per_e_50k_median", "g_per_e_100k_median", "g_per_e_200k_median", "g_per_e_500k_median", "g_per_e_1M_median"),
    labels = c("g_per_e_50k_median" = "50k", "g_per_e_100k_median" = "100k", "g_per_e_200k_median" = "200k", "g_per_e_500k_median" = "500k", "g_per_e_1M_median" = "1M")) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(plot.title = element_text(size = 30))
ggsave("g_per_e_median_p.png", g_per_e_median_p, width = 15, height = 10)
g_per_e_max_p <- ggplot(gene_dis_melt[gene_dis_melt$variable %in% c("g_per_e_1M_max", "g_per_e_500k_max", "g_per_e_200k_max", "g_per_e_100k_max", "g_per_e_50k_max"), ], aes(x = variable, y = value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Max number of annotated genes per enhancer", x = "Species", y = "Number of genes") +
    scale_x_discrete(limits = c("g_per_e_50k_max", "g_per_e_100k_max", "g_per_e_200k_max", "g_per_e_500k_max", "g_per_e_1M_max"),
    labels = c("g_per_e_50k_max" = "50k", "g_per_e_100k_max" = "100k", "g_per_e_200k_max" = "200k", "g_per_e_500k_max" = "500k", "g_per_e_1M_max" = "1M")) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(plot.title = element_text(size = 30))
ggsave("g_per_e_max_p.png", g_per_e_max_p, width = 15, height = 10)
g_per_p_mean_p <- ggplot(gene_dis_melt[gene_dis_melt$variable %in% c("g_per_p_1M_mean", "g_per_p_500k_mean", "g_per_p_200k_mean", "g_per_p_100k_mean", "g_per_p_50k_mean"), ], aes(x = variable, y = value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mean number of annotated genes per promoter", x = "Species", y = "Number of genes") +
    scale_x_discrete(limits = c("g_per_p_50k_mean", "g_per_p_100k_mean", "g_per_p_200k_mean", "g_per_p_500k_mean", "g_per_p_1M_mean"),
    labels = c("g_per_p_50k_mean" = "50k", "g_per_p_100k_mean" = "100k", "g_per_p_200k_mean" = "200k", "g_per_p_500k_mean" = "500k", "g_per_p_1M_mean" = "1M")) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(plot.title = element_text(size = 30))
ggsave("g_per_p_mean_p.png", g_per_p_mean_p, width = 15, height = 10)
g_per_p_median_p <- ggplot(gene_dis_melt[gene_dis_melt$variable %in% c("g_per_p_1M_median", "g_per_p_500k_median", "g_per_p_200k_median", "g_per_p_100k_median", "g_per_p_50k_median"), ], aes(x = variable, y = value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Median number of annotated genes per promoter", x = "Species", y = "Number of genes") +
    scale_x_discrete(limits = c("g_per_p_50k_median", "g_per_p_100k_median", "g_per_p_200k_median", "g_per_p_500k_median", "g_per_p_1M_median"),
    labels = c("g_per_p_50k_median" = "50k", "g_per_p_100k_median" = "100k", "g_per_p_200k_median" = "200k", "g_per_p_500k_median" = "500k", "g_per_p_1M_median" = "1M")) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(plot.title = element_text(size = 30))
ggsave("g_per_p_median_p.png", g_per_p_median_p, width = 15, height = 10)
g_per_p_max_p <- ggplot(gene_dis_melt[gene_dis_melt$variable %in% c("g_per_p_1M_max", "g_per_p_500k_max", "g_per_p_200k_max", "g_per_p_100k_max", "g_per_p_50k_max"), ], aes(x = variable, y = value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Max number of annotated genes per promoter", x = "Species", y = "Number of genes") +
    scale_x_discrete(limits = c("g_per_p_50k_max", "g_per_p_100k_max", "g_per_p_200k_max", "g_per_p_500k_max", "g_per_p_1M_max"),
    labels = c("g_per_p_50k_max" = "50k", "g_per_p_100k_max" = "100k", "g_per_p_200k_max" = "200k", "g_per_p_500k_max" = "500k", "g_per_p_1M_max" = "1M")) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(plot.title = element_text(size = 30))
ggsave("g_per_p_max_p.png", g_per_p_max_p, width = 15, height = 10)
e_distances <- str_split(paste(sp_all[["anno_distance_1M"]][sp_all$element == "enhancer"], collapse = ";"), ";")[[1]]
p_distances <- str_split(paste(sp_all[["anno_distance_1M"]][sp_all$element == "promoter"], collapse = ";"), ";")[[1]]
e_dis <- data.frame(species = sp, ele = "enhancer", distance = e_distances)
p_dis <- data.frame(species = sp, ele = "promoter", distance = p_distances)
distances_df <- rbind(e_dis, p_dis)
distances_df$distance <- as.numeric(distances_df$distance)
distances_df <- distances_df[!is.na(distances_df$distance), ]
ggplot(distances_df[distances_df$distance != 0, ], aes(x = distance / 1000, fill = ele)) +
    geom_histogram(binwidth = 10) +
    geom_vline(xintercept = c(-100, -50, 0, 50, 100), linetype = "dashed", color = "black") +
    theme_bw() +
    labs(title = "Distribution of distances to annotated gene up to 1Mb", x = "Distance (kb)", y = "Density") +
    scale_x_continuous(
        limits = c(-1000, 1200),
        breaks = c(-1000, -500, -200, -100, -50, 0, 50, 100, 200, 500, 1000, 1200)
    ) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 25)) +
    theme(axis.title.x = element_text(size = 25)) +
    theme(axis.title.y = element_text(size = 25)) +
    theme(plot.title = element_text(size = 25))
ggsave("distances_1Mb.png", width = 15, height = 10)
# plot histogram of distances to annotated genes
ggplot(distances_df, aes(x = log10(abs(distance) + 1), fill = ele)) +
    geom_density(alpha = 0.5, ) +
    theme_bw() +
    labs(title = "Distribution of log10(distances+1) to annotated gene up to 1Mb", x = "log10(distances+1)", y = "Density") +
    # scale_x_continuous(limits = c(0, 1200),
    # breaks = c(0, 50, 100, 200, 500, 1000, 1200)) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 25)) +
    theme(axis.title.x = element_text(size = 25)) +
    theme(axis.title.y = element_text(size = 25)) +
    theme(plot.title = element_text(size = 25))
## plot empirical cumulative distribution curve
ggplot(distances_df[distances_df$distance != 0, ], aes(x = abs(distance), color = ele)) +
    stat_ecdf(geom = "smooth", pad = FALSE) +
    theme_bw() +
    labs(title = "Empirical cumulative distribution of distances to annotated gene up to 1Mb", x = "Distance", y = "Density") +
    scale_x_continuous(
        limits = c(0, 1200),
        breaks = c(0, 50, 100, 200, 500, 1000, 1200)
    ) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 25)) +
    theme(axis.title.x = element_text(size = 25)) +
    theme(axis.title.y = element_text(size = 25)) +
    theme(plot.title = element_text(size = 25))
for (sp in species_ts) {
    sp_all <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_all_1M.csv"), header = TRUE)
    # remove "NULL", "", NA
    sp_all$g_per_e_1M <- sapply(sp_all$anno_gene_1M, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    sp_all$g_per_e_500k <- sapply(sp_all$anno_gene_500k, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    sp_all$g_per_e_200k <- sapply(sp_all$anno_gene_200k, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    sp_all$g_per_e_100k <- sapply(sp_all$anno_gene_100k, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    sp_all$g_per_e_50k <- sapply(sp_all$anno_gene_50k, function(x) {
        unique(str_split(x, ";")[[1]] %>%
            .[. != "NULL" & . != "" & !is.na(.)] %>% length())
    })
    write.csv(sp_all, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_all_1M_num.csv"), row.names = FALSE)
    e_distances <- str_split(paste(sp_all[["anno_distance_1M"]][sp_all$element == "enhancer"], collapse = ";"), ";")[[1]]
    p_distances <- str_split(paste(sp_all[["anno_distance_1M"]][sp_all$element == "promoter"], collapse = ";"), ";")[[1]]
    e_dis <- data.frame(species = sp, ele = "enhancer", distance = e_distances)
    p_dis <- data.frame(species = sp, ele = "promoter", distance = p_distances)
    distances_df <- rbind(e_dis, p_dis)
    distances_df$distance <- as.numeric(distances_df$distance)
    distances_df <- distances_df[!is.na(distances_df$distance), ]
    p <- ggplot(distances_df[distances_df$distance != 0, ], aes(x = distance / 1000, fill = ele)) +
        geom_histogram(binwidth = 10) +
        geom_vline(xintercept = c(-100, -50, 0, 50, 100), linetype = "dashed", color = "black") +
        theme_bw() +
        labs(title = "Distribution of distances to annotated gene up to 1Mb (omit 0)", x = "Distance (kb)", y = "count") +
        scale_x_continuous(
            limits = c(-1000, 1200),
            breaks = c(-1000, -500, -200, -100, -50, 0, 50, 100, 200, 500, 1000, 1200)
        ) +
        theme(legend.position = "top") +
        theme(legend.title = element_blank()) +
        theme(legend.text = element_text(size = 20)) +
        theme(axis.text.x = element_text(size = 15)) +
        theme(axis.text.y = element_text(size = 25)) +
        theme(axis.title.x = element_text(size = 25)) +
        theme(axis.title.y = element_text(size = 25)) +
        theme(plot.title = element_text(size = 25))
    ggsave(paste0(sp, "_distances_1Mb.png"), p, width = 15, height = 10)
}
ggplot(distances_df[distances_df$distance != 0, ], aes(x = abs(distance) / 1000, fill = ele)) +
    geom_histogram(alpha = 0.5, binwidth = 10) +
    theme_bw() +
    labs(title = "Distribution of log10(distances+1) to annotated gene up to 1Mb", x = "log10(distances+1)", y = "Density") +
    scale_x_continuous(
        limits = c(0, 1200),
        breaks = c(0, 50, 100, 200, 500, 1000, 1200)
    ) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 25)) +
    theme(axis.title.x = element_text(size = 25)) +
    theme(axis.title.y = element_text(size = 25)) +
    theme(plot.title = element_text(size = 25))
ec <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", tis, "_", ele, "_ec_group.csv"), header = TRUE)
diff_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_diff.csv"), header = TRUE)
egi_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_egi.csv"), header = TRUE)
tapply(egi_df$gene, list(egi_df$note, egi_df$element), length)
egi_count <- data.frame()
for (sp in species_ts) {
    egi_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_egi.csv"), header = TRUE)
    count <- tapply(egi_df$gene, list(egi_df$note, egi_df$element), length)
    count <- data.frame(count)
    count$note <- rownames(count)
    count <- melt(count, id.vars = "note", variable.name = "element", value.name = "count")
    count$species <- sp
    egi_count <- rbind(egi_count, count)
}
egi_count$note <- factor(egi_count$note, levels = rev(c(
    "50k_nearest", "50k",
    "100k_nearest", "100k", "200k_nearest", "200k", "500k_nearest", "500k",
    "1M_nearest", "1M", "nearest"
)))
# stack bar plot
library(ggplot2)
library(tidyr)
stack_dis_p <- ggplot(drop_na(egi_count[egi_count$note != "nearest", ]), aes(x = species, y = count, fill = note)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Number of EGI (element-gene interaction) in each species", x = "Species", y = "Number of EGIs") +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 25)) +
    theme(axis.title.x = element_text(size = 25)) +
    theme(axis.title.y = element_text(size = 25)) +
    theme(plot.title = element_text(size = 25)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(
        "50k_nearest" = "#73b8d5", "50k" = "#73b8d5",
        "100k_nearest" = "#02263e", "100k" = "#02263e",
        "200k_nearest" = "#ed9797", "200k" = "#ed9797",
        "500k_nearest" = "#ef4343", "500k" = "#ef4343",
        "1M_nearest" = "#c01d2e", "1M" = "#c01d2e"
    )) +
    facet_wrap(~element)
ggsave("egi_count_stack.png", width = 15, height = 10)
sp <- "Mus_musculus"
egi_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_egi.csv"), header = TRUE)
ec <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", tis, "_", ele, "_ec_group.csv"), header = TRUE)
high_ec <- ec[ec$count > 10, ]
high_ec_egi <- egi_df[egi_df$consrv >= 10, ]
egi_con_df <- data.frame()
sp_target <- species_ts[species_ts != sp]
high_ec_egi <- data.frame()
for (sp in species_ts) {
    egi_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_egi.csv"), header = TRUE)
    high_ec_egi <- rbind(high_ec_egi, egi_df[egi_df$consrv >= 10, ])
}
# dis_note <- c("50k_nearest", "50k", "100k_nearest", "100k", "200k_nearest", "200k", "500k_nearest", "500k", "1M_nearest", "1M")
# for (tis in unique(high_ec_egi$tissue)) {
#     print(tis)
#     for (ele in unique(high_ec_egi$element)) {
#         print(ele)
#         ec <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", tis, "_", ele, "_ec_group.csv"), header = TRUE)
#         for (peak in unique(high_ec_egi$peak[high_ec_egi$tissue == tis & high_ec_egi$element == ele])) {
#             target_peaks <- ec[ec[[sp]] == peak, sp_target]
#             target_peaks <- target_peaks[, target_peaks != ""]
#             for (sp_t in colnames(target_peaks)) {
#                 target_peak <- target_peaks[[sp_t]] %>%
#                     str_split(., ",")
#                 target_peak <- unlist(target_peak)
#                 target_egi <- high_ec_egi[high_ec_egi$tissue == tis & high_ec_egi$species == sp_t &
#                     high_ec_egi$element == ele & high_ec_egi$peak %in% target_peak & high_ec_egi$note %in% dis_note, ]
#                 source_egi <- high_ec_egi[high_ec_egi$tissue == tis & high_ec_egi$species == sp &
#                     high_ec_egi$element == ele & high_ec_egi$peak == peak & high_ec_egi$note %in% dis_note, ]
#                 og_s <- unique(source_egi$gene)
#                 og_t <- unique(target_egi$gene)
#                 same_og <- intersect(unique(source_egi$og), unique(target_egi$og))
#                 same_ratio <- ifelse(nrow(target_egi) != 0,
#                     (length(same_og) / sqrt(length(og_s) * length(same_og))), NA)
#                 tmp <- data.frame(
#                     peak_s = peak, sp_s = sp, peak_t = target_peak, sp_t = sp_t,
#                     og_s = length(og_s), og_t = length(og_t),
#                     same_og = length(same_og),
#                     same_og_ratio = same_ratio,
#                     tissue = tis, element = ele
#                 )
#                 egi_con_df <- rbind(egi_con_df, tmp)
#             }
#         }
#     }
# }
## too slow
anno_rate_50k <- read.csv("/media/Data/zhangz/chip/analysis/summary2/sum_all/anno_rate_50k.csv", header = TRUE)
