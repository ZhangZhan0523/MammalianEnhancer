library(plyr)
library(treedataverse)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(magrittr)
# library(ggvenn)
# library("ggVennDiagram")
# library("VennDiagram")
library(ggtext)
library(aplot)
library(reshape2)
library(tibble)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggnewscale)
library(stringr)
library(DOSE)
library(tidyr)
library(doMC)
library(ggpubr)
library(ggimage)

doMC::registerDoMC(6)
# load("/media/Data/zhangz/chip/analysis/summary2/abc_1M/egi_con_assess.RData")
setwd("/media/Data/zhangz/chip/analysis/summary2/abc_1M/")
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

ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)
og_mark <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_mark.csv", header = TRUE)
load("/media/Data/zhangz/chip/analysis/expression_mean/lwq/31mm_GEX.RData")

td <- 1
egi_abc_df <- adply(species, 1, function(sp) read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_abc_egi_con_t", td, ".csv"), header = TRUE))
head(egi_abc_df[grep(";", egi_abc_df$same_og_name), ])
highec_egi_df <- egi_abc_df[egi_abc_df$consrv_s >= 10 & !is.na(egi_abc_df$consrv_s), ]
summary(abs(highec_egi_df$distance_diff[!is.na(highec_egi_df$distance_diff) & highec_egi_df$element == "enhancer"]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#       0       0     777  146104  149688 1998196
length(which(abs(highec_egi_df$distance_diff[!is.na(highec_egi_df$distance_diff) & highec_egi_df$element == "enhancer"]) > 1000)) / nrow(highec_egi_df[!is.na(highec_egi_df$distance_diff) & highec_egi_df$element == "enhancer", ])
# 0.4985215
# r$> summary(highec_egi_df$ne_same[highec_egi_df$element == "enhancer"])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.0000  0.0000  0.2078  0.0000  1.0000
# summary(highec_egi_df$same_og[highec_egi_df$element == "enhancer"])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  1.0000  1.0000  0.8437  1.0000  1.0000
og_s_num <- egi_abc_df[, c("ec_id", "peak_s", "sp_s", "og_s", "same_og", "unique_id_s", "distance_s", "note_abc_s", "align_s", "consrv_s", "tissue", "element", "ne_same")]
og_s_num <- og_s_num[!duplicated(og_s_num), ]
## plot distribution of og_s in og_s_num
og_s_num$og_s <- factor(og_s_num$og_s)
og_s_num$og_s <- as.numeric(og_s_num$og_s)
nrow(og_s_num[which(og_s_num$element == "enhancer" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv) & og_s_num$og_s == 1), ]) / nrow(og_s_num[which(og_s_num$element == "enhancer" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ])
# [1] 0.143102 14.3%的高保守元件enhancer只注释一个OG
nrow(og_s_num[which(og_s_num$element == "promoter" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv) & og_s_num$og_s == 1), ]) / nrow(og_s_num[which(og_s_num$element == "promoter" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ])
# [1] 0.2357815 23.6%的高保守元件promoter只注释一个OG
og_s_num_dis <- ggplot(og_s_num[which(og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ], aes(og_s, fill = element)) +
    geom_bar() +
    geom_vline(xintercept = mean(og_s_num[which(og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ]$og_s), linetype = "dashed", color = "black") +
    geom_vline(xintercept = median(og_s_num[which(og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ]$og_s), linetype = "dotted", color = "dark blue") +
    scale_fill_manual(values = c("enhancer" = "#c01d2e", "promoter" = "#033250")) +
    theme_classic() +
    facet_wrap(element ~ .) +
    theme(
        axis.text.x = element_text(hjust = 1),
        text = element_text(size = 30),
        legend.position = "none"
    ) +
    labs(title = "Distribution of annotated OGs by ABC and 1Mts", x = "OGs", y = "Number")
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_s_num_dis.pdf", og_s_num_dis, width = 15, height = 10, dpi = 300)
og_s_num_dis_enh <- ggplot(og_s_num[which(og_s_num$element == "enhancer" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ], aes(og_s)) +
    geom_bar(fill = "#c01d2e") +
    geom_vline(xintercept = mean(og_s_num[which(og_s_num$element == "enhancer" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ]$og_s), linetype = "dashed", color = "black") +
    geom_vline(xintercept = median(og_s_num[which(og_s_num$element == "enhancer" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ]$og_s), linetype = "dotted", color = "dark blue") +
    theme_classic() +
    theme(
        axis.text.x = element_text(hjust = 1),
        text = element_text(size = 30),
        legend.position = "none"
    ) +
    labs(title = "Distribution of annotated OGs of enhancers by ABC and 1Mts", x = "OGs", y = "Number")
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_s_num_dis_enh.pdf", og_s_num_dis_enh, width = 15, height = 10, dpi = 300)
# og_s_num_dis_pro <- ggplot(og_s_num[which(og_s_num$element == "promoter" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)),], aes(og_s)) +
#     geom_bar(fill = "#033250") +
#     geom_vline(xintercept = mean(og_s_num[which(og_s_num$element == "promoter" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)),]$og_s), linetype = "dashed", color = "black") +
#     geom_vline(xintercept = median(og_s_num[which(og_s_num$element == "promoter" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)),]$og_s), linetype = "dotted", color = "dark blue") +
#     theme_classic() +
#     theme(axis.text.x = element_text(hjust = 1),
#     text = element_text(size = 30),
#     legend.position = "none") +
#     labs(title = "Distribution of annotated OGs of promoters by ABC and 1Mts", x = "OGs", y = "Number")
# ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_s_num_dis_pro.pdf", og_s_num_dis_pro, width = 15, height = 10, dpi = 300)
summary(og_s_num[which(og_s_num$element == "promoter" & og_s_num$consrv >= 10 & !is.na(og_s_num$consrv)), ]$og_s)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     1.0     2.0     8.0    12.3    18.0    81.0

unique_id_dfs <- adply(unique(egi_abc_df$sp_s[!is.na(egi_abc_df$consrv_s)]), 1, function(sp) read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_abc_egi_con_t", td, "_unique_id.csv"), header = TRUE))
unique_id_dfs2 <- adply(unique_id_dfs, 1, function(row) {
    tmp <- egi_abc_df[egi_abc_df$sp_s == row$sp & egi_abc_df$unique_id_s == row$unique_id, ]
    ec_sp <- length(unique(tmp$sp_t))
    rc_ec_ratio <- row$same_sp_num / ec_sp
    return(cbind(row, ec_sp = ec_sp, rc_ec_ratio = rc_ec_ratio))
})
high_unique_id_df <- unique_id_dfs[unique_id_dfs$consrv >= 10, ]
nrow(high_unique_id_df[high_unique_id_df$element == "enhancer" & high_unique_id_df$same_og_num <= 3, ]) / nrow(high_unique_id_df[high_unique_id_df$element == "enhancer", ])
# 0.4341825, 43.4%的物种间高保守元件调控的基因在不同物种中只有不超过3个相同OG
nrow(high_unique_id_df[high_unique_id_df$element == "enhancer" & high_unique_id_df$same_sp_num >= 5, ]) / nrow(high_unique_id_df[high_unique_id_df$element == "enhancer", ])
# 0.6030397, 60.3%的物种间高保守元件调控的基因在5个及以上物种中有相同og
nrow(high_unique_id_df[high_unique_id_df$element == "promoter" & high_unique_id_df$same_og_num <= 3, ]) / nrow(high_unique_id_df[high_unique_id_df$element == "promoter", ])
# 0.5660593
nrow(high_unique_id_df[high_unique_id_df$element == "promoter" & high_unique_id_df$same_sp_num >= 5, ]) / nrow(high_unique_id_df[high_unique_id_df$element == "promoter", ])
# 0.6165033
nrow(unique_id_dfs[unique_id_dfs$element == "enhancer" & unique_id_dfs$same_og_num <= 3, ]) / nrow(unique_id_dfs[unique_id_dfs$element == "enhancer", ])
# 0.8852242

# ecdf plot of same_og_num in unique_id_dfs
same_og_num_ecdf1 <- ggplot(unique_id_dfs, aes(same_og_num, color = element)) +
    stat_ecdf(geom = "step", size = 3) + # , position = "jitter"
    scale_color_manual(values = c("#c01d2e", "#033250")) +
    geom_vline(xintercept = 4, linetype = "dashed", color = "black") +
    theme_classic() +
    scale_x_continuous(breaks = c(0, 3, 6, 10)) +
    theme(
        axis.text.x = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2)
    ) +
    labs(title = "ECDF of same OG number of conserved elements", x = "Same OG number", y = "Ratio")
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/same_og_num_ecdf1.pdf", same_og_num_ecdf1, width = 15, height = 10, dpi = 300)

# ecdf plot of same_og_num in high_unique_id_df
same_og_num_ecdf <- ggplot(high_unique_id_df, aes(same_og_num, color = element)) +
    stat_ecdf(geom = "step", size = 3, position = "jitter") +
    scale_color_manual(values = c("#c01d2e", "#033250")) +
    geom_vline(xintercept = 7, linetype = "dashed", color = "black") +
    theme_classic() +
    scale_x_continuous(breaks = c(0, 3, 6, 10)) +
    theme(
        axis.text.x = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2)
    ) +
    labs(title = "ECDF of same OG number of high conserved elements", x = "Same OG number", y = "Ratio")
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/same_og_num_ecdf.png", same_og_num_ecdf, width = 15, height = 10, dpi = 300)
tapply(egi_abc_df$unique_id_s, egi_abc_df$sp_s, function(x) length(unique(x)))
# stat_unique_id_sp <- function(sp) {
#     egi_con_df <- egi_abc_df[egi_abc_df$sp_s == sp, ]
#     resd <- adply(unique(egi_con_df$unique_id_s), 1, function(id) {
#         # for (id in unique(egi_con_df$unique_id_s)) {
#         tmp <- egi_con_df[egi_con_df$unique_id_s == id, ]
#         same_ogs <- unique(tmp$same_og_name[!is.na(tmp$same_og_name)])
#         same_og_freq <- table(tmp$same_og_name[!is.na(tmp$same_og_name)])
#         most_same_og <- names(same_og_freq)[which.max(same_og_freq)]
#         most_same_ratio <- same_og_freq[which.max(same_og_freq)] / sum(same_og_freq)
#         same_og_num <- length(same_ogs)
#         same_sp <- unique(tmp$sp_t[tmp$same_og != 0])
#         same_sp_num <- length(same_sp)
#         diff_sp <- unique(tmp$sp_t[tmp$same_og == 0])
#         diff_sp_num <- length(diff_sp)
#         distance_diff <- paste(tmp$sp_t, tmp$distance_diff[!is.na(tmp$same_og_name)], sep = ":")
#         same_og_pair <- paste(tmp$sp_t, tmp$same_og_name, sep = ":")
#         res <- data.frame(
#             unique_id = id,
#             peak = paste(unique(tmp$peak_s), collapse = ";"),
#             sp = unique(tmp$sp_s),
#             tissue = unique(tmp$tissue),
#             element = unique(tmp$element),
#             align = paste(unique(tmp$align_s), collapse = ";"),
#             # consrv = unique(tmp$consrv_s),
#             distance = paste(unique(tmp$distance_s[!is.na(tmp$same_og_name)]), collapse = ";"),
#             distance_diff = paste(distance_diff, collapse = ";"),
#             same_og_pair = paste(same_og_pair, collapse = ";"),
#             same_og_num = same_og_num,
#             same_ogs = paste(same_ogs, collapse = ";"),
#             most_same_og = most_same_og,
#             most_same_ratio = most_same_ratio,
#             same_sp_num = same_sp_num,
#             same_sp = paste(same_sp, collapse = ";"),
#             diff_sp_num = diff_sp_num,
#             diff_sp = paste(diff_sp, collapse = ";")
#         )
#         # }
#         return(res)
#     })
#     return(resd)
# }
# unique_id_df <- adply(species, 1, stat_unique_id_sp, .parallel = TRUE)
# write.csv(unique_id_df, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/unique_id_df.csv", row.names = FALSE)
unique_id_df <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/unique_id_df.csv", header = TRUE)
tapply(unique_id_dfs$same_og_num, list(unique_id_dfs$sp, unique_id_dfs$same_og_num), length)
# plot distribution of distance_diff in egi_abc_df
dif_distance_enh <- as.numeric(str_split(paste(egi_abc_df[egi_abc_df$element == "enhancer", ]$distance_diff, collapse = ";"), ";")[[1]])
diff_distance_enh <- data.frame(enhancer = dif_distance_enh, n = 1:length(dif_distance_enh))
dif_distance_pro <- as.numeric(str_split(paste(egi_abc_df[egi_abc_df$element == "promoter", ]$distance_diff, collapse = ";"), ";")[[1]])
diff_distance_pro <- data.frame(promoter = dif_distance_pro, n = 1:length(dif_distance_pro))
diff_distance_p <- ggplot(diff_distance_enh[!is.na(diff_distance_enh$enhancer), ], aes(x = log10(abs(enhancer) + 1))) +
    geom_histogram(aes(y = after_stat(count / sum(count))), fill = "#c01d2e", color = "black", position = "dodge", bins = 18, na.rm = TRUE) +
    # geom_histogram(data = diff_distance_pro[!is.na(diff_distance_pro$promoter), ], aes(x = log10(abs(promoter) + 1), y = after_stat(count / sum(count)), fill = "#73b8d5", color = "black", alpha = 0.5), position = "dodge", bins = 18, na.rm = TRUE) +
    theme_classic() +
    # scale_fill_manual(values = c("#c01d2e")) +
    theme(
        axis.text.x = element_text(hjust = 1),
        text = element_text(size = 30),
        legend.position = "none"
    ) +
    labs(title = "Distribution of distance difference of \nsame OGs between species pair (enhancers)", x = "Distance difference (log10(+1))", y = "Ratio")
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/diff_distance_p.png", diff_distance_p, width = 15, height = 10, dpi = 300)

summary(egi_abc_df[which(egi_abc_df$element == "enhancer" & egi_abc_df$distance_diff > 1000), "same_og"])

## 每个物种中有多少基因受保守元件调控，每个基因被多少EC cre调控
same_og_cre_num <- adply(unique(egi_abc_df$sp_s), 1, function(sp) {
    tmp <- egi_abc_df[egi_abc_df$sp_s == sp, ]
    same_og_name <- unique(tmp$same_og_name)
    same_og_name <- same_og_name[!is.na(same_og_name)]
    cre_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x])))
    high_ec_cre_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x & tmp$consrv_s >= 10 & !is.na(tmp$consrv_s)])))
    enh_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x & tmp$element == "enhancer"])))
    high_ec_enh_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x & tmp$element == "enhancer" & tmp$consrv_s >= 10 & !is.na(tmp$consrv_s)])))
    pro_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x & tmp$element == "promoter"])))
    high_ec_pro_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x & tmp$element == "promoter" & tmp$consrv_s >= 10 & !is.na(tmp$consrv_s)])))
    res <- data.frame(
        sp = sp, same_og_name = same_og_name, cre_num = cre_num, enh_num = enh_num, pro_num = pro_num, high_ec_cre_num = high_ec_cre_num, high_ec_enh_num = high_ec_enh_num, high_ec_pro_num = high_ec_pro_num
    )
    return(res)
})

same_og_cre_num <- same_og_cre_num[, -1]
write.csv(same_og_cre_num, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/same_og_cre_num.csv", row.names = FALSE)

high_same_og_cre_num <- adply(unique(highec_egi_df$sp_s), 1, function(sp) {
    tmp <- highec_egi_df[highec_egi_df$sp_s == sp, ]
    same_og_name <- unique(tmp$same_og_name)
    same_og_name <- same_og_name[!is.na(same_og_name)]
    cre_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x])))
    enh_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x & tmp$element == "enhancer"])))
    pro_num <- sapply(same_og_name, function(x) length(unique(tmp$unique_id_s[tmp$same_og_name == x & tmp$element == "promoter"])))
    res <- data.frame(
        sp = sp, same_og_name = same_og_name, cre_num = cre_num, enh_num = enh_num, pro_num = pro_num
    )
    return(res)
})

high_same_og_cre_num <- high_same_og_cre_num[, -1]
write.csv(high_same_og_cre_num, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/high_same_og_cre_num.csv", row.names = FALSE)

# 匹配调控熵
ent1_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ent_ouwie/longevity/ent1_all.csv", header = TRUE)
ent2_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ent_ouwie/longevity/ent2_all.csv", header = TRUE)
same_og_cre_ent <- adply(same_og_cre_num, 1, function(row) {
    same_og_name <- row$same_og_name
    cre_num <- row$cre_num
    enh_num <- row$enh_num
    pro_num <- row$pro_num
    ent1 <- sapply(same_og_name, function(x) {
        tmp <- ent1_all[ent1_all$OGID == x, row$sp]
        if (length(tmp) == 0) {
            return(NA)
        }
        return(tmp)
    })
    ent2 <- sapply(same_og_name, function(x) {
        tmp <- ent2_all[ent2_all$OGID == x, row$sp]
        if (length(tmp) == 0) {
            return(NA)
        }
        return(tmp)
    })
    res <- data.frame(
        sp = row$sp, same_og_name = same_og_name, cre_num = cre_num, enh_num = enh_num, pro_num = pro_num,
        ent1 = ent1, ent2 = ent2
    )
    return(res)
})
# sp <- "Ovis_aries"
a_ply(species, 1, function(sp) {
    ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
    ent <- merge(ent, same_og_cre_num[same_og_cre_num$sp == sp, c("same_og_name", "cre_num", "enh_num", "pro_num")],
        by.x = "OGID", by.y = "same_og_name", all.x = TRUE
    )
    ent <- replace(ent, is.na(ent), 0)
    cor.test(ent$ent1[ent$ent2 != 0], ent$cre_num[ent$ent2 != 0], method = "spearman")
    # whether ogs with EC element have higher entropy
    t.test(ent$ent1[ent$cre_num != 0 & ent$ent2 != 0], ent$ent1[ent$cre_num == 0 & ent$ent2 != 0], alternative = "greater")
    chisq.test(matrix(c(
        sum(ent$cre_num != 0 & ent$ent2 != 0), sum(ent$cre_num == 0 & ent$ent2 != 0),
        sum(ent$cre_num != 0 & ent$ent2 == 0), sum(ent$cre_num == 0 & ent$ent2 == 0)
    ), nrow = 2), correct = FALSE)
    # violin plot of ent1, grouped by whether cre_num is 0, add p value
    pairs <- list(c("TRUE", "FALSE"))
    ent_violin <- ggplot(ent[ent$ent2 != 0, ], aes(x = factor(cre_num != 0), y = ent1, fill = factor(cre_num != 0))) +
        geom_violin(alpha = 0.5) +
        geom_boxplot(width = 0.1) +
        stat_compare_means(comparisons = pairs, method = "wilcox.test") +
        scale_fill_manual(values = c("TRUE" = "#c01d2e", "FALSE" = "#033250")) +
        theme_classic() +
        theme(
            axis.text.x = element_text(hjust = 0.5),
            text = element_text(size = 30),
            legend.position = "none"
        ) +
        labs(title = paste0("Entropy distribution of OGs with and without EC element \nin ", sp), x = "Whether with EC element", y = "Entropy")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/ent_violin_", sp, ".png"), ent_violin, width = 15, height = 10, dpi = 300)
})
same_og_freq <- table(same_og_cre_num$same_og_name[!is.na(same_og_cre_num$same_og_name)])
## order same_og_freq
same_og_freq <- same_og_freq[order(same_og_freq, decreasing = TRUE)]
most_same_og <- names(same_og_freq)[which(same_og_freq == 25)]
## 在物种间广泛被保守元件调控的og有哪些，有24个，2个抑癌基因，2个hk
ortho_all[ortho_all$Orthogroup %in% most_same_og, "Mus_musculus"]
og_mark[og_mark$Orthogroup %in% most_same_og, ]
## 准备做富集分析

## 在每个物种中收到保守元件调控的基因有哪些，表达的基因有哪些，所占百分比
sp_same_og_num <- tapply(same_og_cre_num$same_og_name, same_og_cre_num$sp, length)
sp_same_og_num <- data.frame(sp = names(sp_same_og_num), num = sp_same_og_num)
high_sp_same_og_num <- tapply(high_same_og_cre_num$same_og_name, high_same_og_cre_num$sp, length)
high_sp_same_og_num <- data.frame(sp = names(high_sp_same_og_num), num = high_sp_same_og_num)

expre_gene_num <- adply(unique(same_og_cre_num$sp), 1, function(sp) {
    print(sp)
    tmp <- raw_FPKM_me_0_sep[raw_FPKM_me_0_sep$species == sp, ]
    tmp <- tmp[!is.na(tmp$Brain), ]
    tmp <- tmp[!(tmp$Brain < 1 & tmp$Kidney < 1 & tmp$Liver < 1), ]
    if (sp %in% species_nots) {
        tmp <- ortho_all[, c("Orthogroup", sp)]
        tmp <- tmp[tmp[[sp]] != "NULL", ]
        tmp <- tmp[tmp[[sp]] != "", ]
    }
    return(data.frame(sp = sp, gene_num = nrow(tmp)))
})
expre_gene_num <- expre_gene_num[, -1]
sp_same_og_num <- merge(sp_same_og_num, expre_gene_num, by = "sp")
sp_same_og_num$rc_ratio <- sp_same_og_num$num / sp_same_og_num$gene_num

colnames(sp_same_og_num) <- c("sp", "num", "gene_num", "rc_ratio", "num_high_ec")
sp_same_og_num$high_ec_ratio <- sp_same_og_num$num_high_ec / sp_same_og_num$num
sp_same_og_num <- merge(sp_same_og_num, high_sp_same_og_num, by = "sp", all = TRUE)
# melt
sp_same_og_num_melt <- melt(sp_same_og_num[, c("sp", "rc_ratio", "high_ec_ratio")], id.vars = "sp", measure.vars = c("rc_ratio", "high_ec_ratio"))
## box and violin plot of sp_same_og_num$rc_ratio

sp_same_og_num_box <- ggplot(sp_same_og_num_melt, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(alpha = 0.5) +
    # geom_violin() +
    theme_classic() +
    theme(
        axis.text.x = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.position = "none"
    ) +
    labs(title = "Ratio of genes regulated \nby conserved elements", y = "Ratio")
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/sp_same_og_num_ratio_box.png", sp_same_og_num_box, width = 8, height = 10, dpi = 300)

egi_nearest_df <- adply(species_ts, 1, function(x) read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/egi_con/", x, "_nearest_egi_con.csv"), header = TRUE))
count(egi_abc_df$element)
tapply(egi_abc_df$element, egi_abc_df$sp_s, count)
egi_nearest_df$element <- factor(egi_nearest_df$element, levels = c("enhancer", "promoter"))
count(egi_nearest_df$element)
tapply(egi_abc_df$element, list(egi_abc_df$element, egi_abc_df$sp_s), length)
tapply(egi_nearest_df$element, list(egi_nearest_df$element, egi_nearest_df$sp_s), length)
peak_num <- tapply(egi_abc_df$peak_s, list(egi_abc_df$sp_s, egi_abc_df$tissue, egi_abc_df$element), function(x) unique(length(x))) %>% as.data.frame()
enh_num <- peak_num[, 1:3, 1] %>% as.data.frame()
enh_num$sum <- rowSums(enh_num, na.rm = TRUE)
pro_num <- peak_num[, 1:3, 2] %>% as.data.frame()
pro_num$sum <- rowSums(pro_num, na.rm = TRUE)
all_peak_num <- tapply(egi_abc_df$peak_s, list(egi_abc_df$sp_s, egi_abc_df$tissue), function(x) unique(length(x))) %>% as.data.frame()
all_peak_num_nearest <- tapply(egi_nearest_df$peak_s, list(egi_nearest_df$sp_s, egi_nearest_df$tissue), function(x) unique(length(x))) %>% as.data.frame()

all_peak_num$sum <- rowSums(all_peak_num, na.rm = TRUE)
all_peak_num_nearest$sum <- rowSums(all_peak_num_nearest, na.rm = TRUE)
egi_sp <- expand.grid(species_ts, species_ts)
colnames(egi_sp) <- c("sp_s", "sp_t")
get_average_ratio <- function(egi_sp, egi_abc_df) {
    sp_s <- egi_sp$sp_s[1]
    sp_t <- egi_sp$sp_t[1]
    tmp <- egi_abc_df[egi_abc_df$sp_s == sp_s & egi_abc_df$sp_t == sp_t, ]
    average_ratio <- ifelse(nrow(tmp) == 0, 0, mean(tmp$same_og_ratio, na.rm = TRUE))
    average_og_s <- ifelse(nrow(tmp) == 0, 0, mean(tmp$og_s, na.rm = TRUE))
    average_og_t <- ifelse(nrow(tmp) == 0, 0, mean(tmp$og_t, na.rm = TRUE))
    average_same_og <- ifelse(nrow(tmp) == 0, 0, mean(tmp$same_og, na.rm = TRUE))
    tmp_enh <- tmp[tmp$element == "enhancer", ]
    tmp_pro <- tmp[tmp$element == "promoter", ]
    average_ratio_enh <- ifelse(nrow(tmp_enh) == 0, 0, mean(tmp_enh$same_og_ratio, na.rm = TRUE))
    average_ratio_pro <- ifelse(nrow(tmp_pro) == 0, 0, mean(tmp_pro$same_og_ratio, na.rm = TRUE))
    average_og_s_enh <- ifelse(nrow(tmp_enh) == 0, 0, mean(tmp_enh$og_s, na.rm = TRUE))
    average_og_t_enh <- ifelse(nrow(tmp_enh) == 0, 0, mean(tmp_enh$og_t, na.rm = TRUE))
    average_og_s_pro <- ifelse(nrow(tmp_pro) == 0, 0, mean(tmp_pro$og_s, na.rm = TRUE))
    average_og_t_pro <- ifelse(nrow(tmp_pro) == 0, 0, mean(tmp_pro$og_t, na.rm = TRUE))
    average_same_og_enh <- ifelse(nrow(tmp_enh) == 0, 0, mean(tmp_enh$same_og, na.rm = TRUE))
    average_same_og_pro <- ifelse(nrow(tmp_pro) == 0, 0, mean(tmp_pro$same_og, na.rm = TRUE))
    res <- data.frame(
        sp_s = sp_s, sp_t = sp_t, average_ratio = average_ratio, average_og_s = average_og_s, average_og_t = average_og_t,
        average_ratio_enh = average_ratio_enh, average_ratio_pro = average_ratio_pro, average_og_s_enh = average_og_s_enh,
        average_og_t_enh = average_og_t_enh, average_og_s_pro = average_og_s_pro, average_og_t_pro = average_og_t_pro,
        average_same_og = average_same_og, average_same_og_enh = average_same_og_enh, average_same_og_pro = average_same_og_pro,
        n = nrow(tmp), n_enh = nrow(tmp_enh), n_pro = nrow(tmp_pro)
    )
    return(res)
}
egi_sp_df <- adply(egi_sp, 1, function(x) get_average_ratio(x, egi_abc_df = egi_abc_df), .parallel = TRUE)
egi_sp_nearest_df <- adply(egi_sp, 1, function(x) get_average_ratio(x, egi_abc_df = egi_nearest_df), .parallel = TRUE)

write.csv(egi_sp_df, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/egi_sp_df.csv", row.names = FALSE)
write.csv(egi_sp_nearest_df, "/media/Data/zhangz/chip/analysis/summary2/sum_all/egi_con/egi_sp_nearest_df.csv", row.names = FALSE)
combinations <- combn(species_ts, 2)

get_half_mean <- function(com, egi_sp_df) {
    sp_s <- com[1]
    sp_t <- com[2]
    tmp1 <- egi_sp_df[egi_sp_df$sp_s == sp_s & egi_sp_df$sp_t == sp_t, ]
    tmp2 <- egi_sp_df[egi_sp_df$sp_s == sp_t & egi_sp_df$sp_t == sp_s, ]
    average_ratio <- mean(c(tmp1$average_ratio, tmp2$average_ratio), na.rm = TRUE)
    average_og_s <- mean(c(tmp1$average_og_s, tmp2$average_og_t), na.rm = TRUE)
    average_og_t <- mean(c(tmp1$average_og_t, tmp2$average_og_s), na.rm = TRUE)
    average_same_og <- mean(c(tmp1$average_same_og, tmp2$average_same_og), na.rm = TRUE)
    average_ratio_enh <- mean(c(tmp1$average_ratio_enh, tmp2$average_ratio_enh), na.rm = TRUE)
    average_ratio_pro <- mean(c(tmp1$average_ratio_pro, tmp2$average_ratio_pro), na.rm = TRUE)
    average_og_s_enh <- mean(c(tmp1$average_og_s_enh, tmp2$average_og_t_enh), na.rm = TRUE)
    average_og_t_enh <- mean(c(tmp1$average_og_t_enh, tmp2$average_og_s_enh), na.rm = TRUE)
    average_og_s_pro <- mean(c(tmp1$average_og_s_pro, tmp2$average_og_t_pro), na.rm = TRUE)
    average_og_t_pro <- mean(c(tmp1$average_og_t_pro, tmp2$average_og_s_pro), na.rm = TRUE)
    average_same_og_enh <- mean(c(tmp1$average_same_og_enh, tmp2$average_same_og_enh), na.rm = TRUE)
    average_same_og_pro <- mean(c(tmp1$average_same_og_pro, tmp2$average_same_og_pro), na.rm = TRUE)
    n <- mean(c(tmp1$n, tmp2$n), na.rm = TRUE)
    n_enh <- mean(c(tmp1$n_enh, tmp2$n_enh), na.rm = TRUE)
    n_pro <- mean(c(tmp1$n_pro, tmp2$n_pro), na.rm = TRUE)
    res <- data.frame(
        sp_s = sp_s, sp_t = sp_t, average_ratio = average_ratio, average_og_s = average_og_s, average_og_t = average_og_t,
        average_ratio_enh = average_ratio_enh, average_ratio_pro = average_ratio_pro, average_og_s_enh = average_og_s_enh,
        average_og_t_enh = average_og_t_enh, average_og_s_pro = average_og_s_pro, average_og_t_pro = average_og_t_pro,
        average_same_og = average_same_og, average_same_og_enh = average_same_og_enh, average_same_og_pro = average_same_og_pro,
        n = n, n_enh = n_enh, n_pro = n_pro
    )
    return(res)
}
egi_sp_df_half <- adply(combinations, 2, function(x) get_half_mean(x, egi_sp_df = egi_sp_df), .parallel = TRUE)
write.csv(egi_sp_df_half, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/egi_sp_df_half.csv", row.names = FALSE)
egi_sp_nearest_df_half <- adply(combinations, 2, function(x) get_half_mean(x, egi_sp_df = egi_sp_nearest_df), .parallel = TRUE)
exchange_sp_name <- function(com, egi_sp_df_half) {
    sp_1 <- com$sp1[1]
    sp_2 <- com$sp2[1]
    tmp <- which(egi_sp_df_half$sp_s == sp_1 & egi_sp_df_half$sp_t == sp_2)
    res <- egi_sp_df_half
    res[tmp, "sp_t"] <- sp_1
    res[tmp, "sp_s"] <- sp_2
    return(res)
}
name_com <- data.frame(
    sp1 = c(
        "Macaca_mulatta", "Lama_glama", "Sus_scrofa", "Bos_taurus", "Ovis_aries",
        "Lama_glama", "Sus_scrofa", "Bos_taurus", "Ovis_aries",
        "Lama_glama", "Sus_scrofa", "Bos_taurus", "Ovis_aries",
        "Lama_glama", "Sus_scrofa", "Bos_taurus", "Ovis_aries",
        "Lama_glama", "Sus_scrofa", "Bos_taurus", "Ovis_aries"
    ),
    sp2 = c(
        "Tupaia_belangeri", "Equus_caballus", "Equus_caballus", "Equus_caballus", "Equus_caballus",
        "Equus_asinus", "Equus_asinus", "Equus_asinus", "Equus_asinus",
        "Felis_catus", "Felis_catus", "Felis_catus", "Felis_catus",
        "Canis_lupus", "Canis_lupus", "Canis_lupus", "Canis_lupus",
        "Mustela_putorius", "Mustela_putorius", "Mustela_putorius", "Mustela_putorius"
    )
)
for (i in seq_len(nrow(name_com))) {
    egi_sp_df_half <- exchange_sp_name(name_com[i, ], egi_sp_df_half)
}
for (i in seq_len(nrow(name_com))) {
    egi_sp_nearest_df_half <- exchange_sp_name(name_com[i, ], egi_sp_nearest_df_half)
}
egi_sp_df_half <- egi_sp_df_half[, -1]
egi_sp_nearest_df_half <- egi_sp_nearest_df_half[, -1]
egi_sp_df_half <- rbind(egi_sp_df_half, data.frame(
    sp_s = "Petaurus_breviceps", sp_t = "Petaurus_breviceps",
    average_ratio = 0, average_og_s = 0, average_og_t = 0,
    average_ratio_enh = 0, average_ratio_pro = 0, average_og_s_enh = 0,
    average_og_t_enh = 0, average_og_s_pro = 0, average_og_t_pro = 0,
    average_same_og = 0, average_same_og_enh = 0, average_same_og_pro = 0,
    n = 0, n_enh = 0, n_pro = 0
))
egi_sp_nearest_df_half <- rbind(egi_sp_nearest_df_half, data.frame(
    sp_s = "Petaurus_breviceps", sp_t = "Petaurus_breviceps",
    average_ratio = 0, average_og_s = 0, average_og_t = 0,
    average_ratio_enh = 0, average_ratio_pro = 0, average_og_s_enh = 0,
    average_og_t_enh = 0, average_og_s_pro = 0, average_og_t_pro = 0,
    average_same_og = 0, average_same_og_enh = 0, average_same_og_pro = 0,
    n = 0, n_enh = 0, n_pro = 0
))
egi_nearest_ratio_heat <- ggplot(egi_sp_nearest_df_half, aes(sp_s, sp_t, fill = average_ratio)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "#02263e") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Average ratio of nearest OGs", x = "Source species", y = "Target species")
egi_ratio_heat <- ggplot(egi_sp_df_half, aes(sp_s, sp_t, fill = average_ratio)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "#02263e") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Average ratio of same OGs", x = "EGI predicted by ABC and 1Mts", y = "EGI predicted by 'nearest' method") +
    theme(legend.position = c(0.8, -0.2))
sp_tre <- read.newick("/media/Data/zhangz/chip/scripts2/info/sps.nwk", node.label = "label")
sp_ts_tre <- keep.tip(sp_tre, species_ts)
tre_y <- ggtree(sp_ts_tre, layout = "rectangular")
tre_x <- ggtree(sp_ts_tre, layout = "rectangular") + coord_flip()
p <- egi_ratio_heat %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)
p_nearest <- egi_nearest_ratio_heat %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)
test_heat <- egi_ratio_heat + geom_tile(data = egi_sp_nearest_df_half, aes(sp_t, sp_s, fill = average_ratio))
test_p <- test_heat %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/egi_ratio_heat.pdf", test_p, width = 10, height = 10)
egi_enh_ratio_heat <- ggplot(egi_sp_df_half, aes(sp_s, sp_t, fill = average_ratio_enh)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "#c01d2e") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Average ratio of same OGs in enhancers", x = "Source species", y = "Target species") +
    geom_text(aes(label = round(average_ratio_enh, 2)))
p_enh <- egi_enh_ratio_heat %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)
egi_pro_ratio_heat <- ggplot(egi_sp_df_half, aes(sp_s, sp_t, fill = average_ratio_pro)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "#02263e") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Average ratio of same OGs in promoters", x = "Source species", y = "Target species")
p_pro <- egi_pro_ratio_heat %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)
enh_heat_test <- egi_enh_ratio_heat + geom_tile(data = egi_sp_nearest_df_half, aes(sp_t, sp_s, fill = average_ratio_enh)) +
    scale_fill_gradient2(low = "white", high = "#c01d2e") +
    geom_text(data = egi_sp_nearest_df_half, aes(sp_t, sp_s, label = round(average_ratio_enh, 2)))
enh_test_p <- enh_heat_test %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/egi_enh_ratio_heat.pdf", enh_test_p, width = 10, height = 10)
pro_heat_test <- egi_pro_ratio_heat + geom_tile(data = egi_sp_nearest_df_half, aes(sp_t, sp_s, fill = average_ratio_pro)) +
    scale_fill_gradient2(low = "white", high = "#02263e")
pro_test_p <- pro_heat_test %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/egi_pro_ratio_heat.pdf", pro_test_p, width = 10, height = 10)

same_og <- unique(str_split(paste(egi_abc_df$same_og_name, collapse = ";"), ";")[[1]])
same_og <- same_og[same_og != ""]
same_og_ep <- tapply(egi_abc_df$same_og_name, egi_abc_df$element, function(x) unique(str_split(paste(x, collapse = ";"), ";")[[1]]))
same_og_ep <- lapply(same_og_ep, function(x) x[x != ""])
# $enhancer
# [1] 45

# $promoter
# [1] 5753
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/fpkm_ortho.RData")
colnames(ortho_me_adj) <- gsub("Macaca_fascicularis", "Macaca_mulatta", colnames(ortho_me_adj))
raw_FPKM_me_0_sep$species <- gsub("Macaca_fascicularis", "Macaca_mulatta", raw_FPKM_me_0_sep$species)
length(which(same_og_ep$promoter %in% ortho_me_20_adj$id))
length(which(same_og_ep$enhancer %in% ortho_me_20_adj$id))
hk_og_all <- read.csv("/media/Data/zhangz/chip/analysis/hk/hk_og_all.csv", header = TRUE)
length(which(same_og_ep$promoter %in% hk_og_all$OGID))
length(which(same_og_ep$enhancer %in% hk_og_all$OGID))
ortho_all <- read.csv("/data1/Genome/orthofinder/Results_Oct22/Orthogroups/Orthogroups.tsv", sep = "\t", header = TRUE)
colnames(ortho_all) <- gsub(".uniq.pr", "", colnames(ortho_all))

# mark og with cancer genes, housekeeping genes, and essential genes
onco_g <- read.csv("/media/Data/zhangz/chip/analysis/onco/cancerGeneList.tsv", sep = "\t", header = TRUE)
og_mark <- data.frame(
    Orthogroup = ortho_all$Orthogroup, cancer = FALSE, onco = FALSE, TSG = FALSE,
    hk = FALSE, essential = FALSE, o2o20 = FALSE
)
find_onco_og <- function(onco_row) {
    sym <- onco_row$Hugo.Symbol
    if (TRUE %in% grepl(paste0("(^|, )", sym, "[|]"), ortho_all$Homo_sapiens)) {
        res <- ortho_all[grep(paste0("(^|, )", sym, "[|]"), ortho_all$Homo_sapiens), "Orthogroup"]
    } else {
        for (g in str_split(onco_row$Gene.Aliases, ", ")[[1]]) {
            if (TRUE %in% grepl(paste0("(^|, )", g, "[|]"), ortho_all$Homo_sapiens)) {
                res <- ortho_all[grep(paste0("(^|, )", g, "[|]"), ortho_all$Homo_sapiens), "Orthogroup"]
                break
            }
        }
    }
    return(res)
}

for (i in seq_len(nrow(onco_g))) {
    og <- find_onco_og(onco_g[i, ])
    og_mark[og_mark$Orthogroup %in% og, "cancer"] <- TRUE
    og_mark[og_mark$Orthogroup %in% og, "symbol"] <- onco_g[i, "Hugo.Symbol"]
    if (onco_g[i, "Is.Oncogene"] == "Yes") {
        og_mark[og_mark$Orthogroup %in% og, "onco"] <- TRUE
    }
    if (onco_g[i, "Is.Tumor.Suppressor.Gene"] == "Yes") {
        og_mark[og_mark$Orthogroup %in% og, "TSG"] <- TRUE
    }
}
# mark hk gene
load("/media/Data/zhangz/chip/analysis/hk/Housekeeping_Genes_Mouse.RData")
for (i in seq_len(nrow(Mouse_HK_genes))) {
    hk_ogid <- ortho_all[grep(paste0("(^|, )", Mouse_HK_genes$Gene[i], "[|]"), ortho_all$Mus_musculus), "Orthogroup"]
    og_mark[og_mark$Orthogroup %in% hk_ogid, "hk"] <- TRUE
}
# mark ogs have one2one orthologs in over 20 species
og_mark$o2o20[og_mark$Orthogroup %in% ortho_me_20_adj$id] <- TRUE
write.csv(og_mark, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_mark.csv", row.names = FALSE)
tapply(og_mark$cancer, og_mark$o2o20, count)
tapply(og_mark$onco, og_mark$o2o20, count)
tapply(og_mark$TSG, og_mark$o2o20, count)
tapply(og_mark$hk, og_mark$o2o20, count)

library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "asia")
og_to_mgene <- function(og) {
    # get mouse gene from ortho
    mgene <- ortho_all[ortho_all$Orthogroup == og, "Mus_musculus"] %>% str_split(., ",")
    # if (mgene[[1]][1] == "NULL" | mgene[[1]][1] == "") {
    #     mgene <- ortho_all[ortho_all$Orthogroup == og, "Rattus_norvegicus"] %>% str_split(., ",")
    # } else {
    #     mgene <- mgene[[1]][1] %>% str_split(., "[|]")
    #     res <- mgene[[1]][1]
    # }
    # if (mgene[[1]][1] == "NULL" | mgene[[1]][1] == "") {
    #     mgene <- ortho_all[ortho_all$Orthogroup == og, "Homo_sapiens"] %>%
    #         str_split(., ",")
    #     res <- getLDS(attributes = "hgnc_symbol", filters = "hgnc_symbol", values = str_split(mgene, "[|]")[[1]][1], mart = human, attributesL = "mgi_symbol", martL = mouse, uniqueRows = TRUE)
    # } else {
    #     mgene <- mgene[[1]][1] %>% str_split(., "[|]")
    #     res <- mgene[[1]][2]
    # }
    # if ()
    # mgene <- mgene[[1]][1] %>% str_split(., "[|]")
    # return(mgene[[1]][1])
    mgene <- mgene[[1]][1] %>% str_split(., "[|]")
    res <- mgene[[1]][1]
    return(res)
}
n <- 0
for (og in nog) {
    if (og_to_mgene(og) == "") {
        print(og_to_mgene(og))
        n <- n + 1
    }
}
same_og_ep_mgene <- lapply(same_og_ep, function(x) sapply(x, og_to_mgene))
same_og_ep_mgene <- lapply(same_og_ep_mgene, function(x) x[x != ""])
# convert SYMBOL to ENTREZID

same_og_ep_mgene_id <- lapply(same_og_ep_mgene, function(x) {
    sapply(x, function(y) {
        if (y == "") {
            return("")
        }
        res <- try(bitr(y, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
        if (length(res) == 0) {
            return("")
        }
        if (class(res) == "try-error") {
            print(y)
            res <- ""
        } else {
            res$ENTREZID
        }
        return(res)
    })
})
same_og_ep_mgene_id <- lapply(same_og_ep_mgene_id, function(x) x[x != ""])
same_og_ep_mgene_id_bak <- same_og_ep_mgene_id
# same_og_ep_mgene_id <- lapply(same_og_ep_mgene_id, function(x) t(x) %>% as.data.frame)
kegg_l <- lapply(same_og_ep_mgene_id, function(x) enrichKEGG(as.numeric(unlist(x)), "mmu", keyType = "ENTREZID", use_internal_data = TRUE))
dotplot(kegg_pro, showCategory = 20)
go_l <- lapply(same_og_ep_mgene_id, function(x) enrichGO(as.numeric(unlist(x)), OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP"))
dotplot(kegg_l$promoter, showCategory = 10)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/kegg_pro.pdf", dotplot(kegg_l$promoter, showCategory = 10), width = 10, height = 10)
# dotplot(kegg_l$enhancer, showCategory = 1)
dotplot(go_l$enhancer, showCategory = 20)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/go_enh.pdf", dotplot(go_l$enhancer, showCategory = 20), width = 10, height = 10)
dotplot(go_l$promoter, showCategory = 20)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/go_pro.pdf", dotplot(go_l$promoter, showCategory = 20), width = 10, height = 10)
ncg_l <- lapply(same_og_ep_mgene_id, function(x) enrichNCG(as.numeric(unlist(x))))
# nog <- names(same_og_ep_mgene$promoter[same_og_ep_mgene$promoter == "NULL"])
# northo <- ortho[ortho$id %in% nog, ]
# non_null_cols <- apply(northo[, -c(1, 2)], 2, function(x) length(which(x != "NULL")))
# which.max(non_null_cols)

library(phytools)
# phylosig(sp_tre, ent1, method = "K", test = TRUE)

td <- 1
egi_abc_df_td1 <- adply(species_ts, 1, function(x) read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", x, "_abc_egi_con_t", td, ".csv"), header = TRUE))
tapply(egi_abc_df_td1$peak_s, list(egi_abc_df_td1$sp_s, egi_abc_df_td1$tissue, egi_abc_df_td1$element), function(x) unique(length(x))) %>% as.data.frame()
egi_sp_df_td1 <- adply(egi_sp, 1, function(x) get_average_ratio(x, egi_abc_df = egi_abc_df_td1), .parallel = TRUE)
egi_sp_df_half_td1 <- adply(combinations, 2, function(x) get_half_mean(x, egi_sp_df = egi_sp_df_td1), .parallel = TRUE)
for (i in seq_len(nrow(name_com))) {
    egi_sp_df_half_td1 <- exchange_sp_name(name_com[i, ], egi_sp_df_half_td1)
}

egi_nearest_df_td1 <- adply(species_ts, 1, function(x) read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/egi_con/", x, "_nearest_egi_con_t", td, ".csv"), header = TRUE))
egi_nearest_df_td1$element <- factor(egi_nearest_df_td1$element, levels = c("enhancer", "promoter"))
count(egi_nearest_df_td1$element)
tapply(egi_nearest_df_td1$peak_s, list(egi_nearest_df_td1$sp_s, egi_nearest_df_td1$tissue, egi_nearest_df_td1$element), function(x) unique(length(x))) %>% as.data.frame()
egi_sp_nearest_df_td1 <- adply(egi_sp, 1, function(x) get_average_ratio(x, egi_abc_df = egi_nearest_df_td1), .parallel = TRUE)
egi_sp_nearest_df_half_td1 <- adply(combinations, 2, function(x) get_half_mean(x, egi_sp_df = egi_sp_nearest_df_td1), .parallel = TRUE)
for (i in seq_len(nrow(name_com))) {
    egi_sp_nearest_df_half_td1 <- exchange_sp_name(name_com[i, ], egi_sp_nearest_df_half_td1)
}

# what gene are in same og
same_og_td1 <- unique(str_split(paste(egi_abc_df_td1$same_og_name, collapse = ";"), ";")[[1]])
same_og_td1 <- same_og_td1[same_og_td1 != ""]
same_og_ep_td1 <- tapply(egi_abc_df_td1$same_og_name, egi_abc_df_td1$element, function(x) unique(str_split(paste(x, collapse = ";"), ";")[[1]]))
same_og_ep_td1 <- lapply(same_og_ep_td1, function(x) x[x != ""])
# r$> lapply(same_og_ep_td1, length)
# $enhancer
# [1] 799

# $promoter
# [1] 7207
same_og_ep_mgene_td1 <- lapply(same_og_ep_td1, function(x) sapply(x, og_to_mgene))
same_og_ep_mgene_td1 <- lapply(same_og_ep_mgene_td1, function(x) x[x != ""])
same_og_ep_mgene_id_td1 <- llply(same_og_ep_mgene_td1, function(x) {
    sapply(x, function(y) {
        if (y == "") {
            return("")
        }
        res <- try(bitr(y, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
        if (length(res) == 0) {
            return("")
        }
        if (class(res) == "try-error") {
            print(y)
            res <- ""
        } else {
            res$ENTREZID
        }
        return(res)
    })
}, .parallel = TRUE)
same_og_ep_mgene_id_td1 <- lapply(same_og_ep_mgene_id_td1, function(x) x[x != ""])
same_og_ep_mgene_id_td1_bak <- same_og_ep_mgene_id_td1
og_mark <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_mark.csv", header = TRUE)

# calc enrichment foldchange and p-value
enh_onco_td10_df <- data.frame(matrix(c(
    length(which(same_og_ep$enhancer %in% og_mark$Orthogroup[og_mark$onco])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep$enhancer) & og_mark$onco)),
    length(which(same_og_ep$enhancer %in% og_mark$Orthogroup[!og_mark$onco])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep$enhancer) & !og_mark$onco))
), 2, 2))
enh_onco_td10_fc <- (enh_onco_td10_df[1, 1] / (enh_onco_td10_df[1, 1] + enh_onco_td10_df[1, 2])) /
    ((enh_onco_td10_df[1, 1] + enh_onco_td10_df[2, 1]) / (enh_onco_td10_df[1, 1] + enh_onco_td10_df[1, 2] + enh_onco_td10_df[2, 1] + enh_onco_td10_df[2, 2]))
enh_onco_td10 <- fisher.test(enh_onco_td10_df, alternative = "greater")
pro_onco_td10_df <- data.frame(matrix(c(
    length(which(same_og_ep$promoter %in% og_mark$Orthogroup[og_mark$onco])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep$promoter) & og_mark$onco)),
    length(which(same_og_ep$promoter %in% og_mark$Orthogroup[!og_mark$onco])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep$promoter) & !og_mark$onco))
), 2, 2))
pro_onco_td10_fc <- (pro_onco_td10_df[1, 1] / (pro_onco_td10_df[1, 1] + pro_onco_td10_df[1, 2])) /
    ((pro_onco_td10_df[1, 1] + pro_onco_td10_df[2, 1]) / (pro_onco_td10_df[1, 1] + pro_onco_td10_df[1, 2] + pro_onco_td10_df[2, 1] + pro_onco_td10_df[2, 2]))
pro_onco_td10 <- fisher.test(pro_onco_td10_df, alternative = "greater")

enh_onco_td1_df <- data.frame(matrix(c(
    length(which(same_og_ep_td1$enhancer %in% og_mark$Orthogroup[og_mark$onco])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep_td1$enhancer) & og_mark$onco)),
    length(which(same_og_ep_td1$enhancer %in% og_mark$Orthogroup[!og_mark$onco])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep_td1$enhancer) & !og_mark$onco))
), 2, 2))
enh_onco_td1_fc <- (enh_onco_td1_df[1, 1] / (enh_onco_td1_df[1, 1] + enh_onco_td1_df[1, 2])) /
    ((enh_onco_td1_df[1, 1] + enh_onco_td1_df[2, 1]) / (enh_onco_td1_df[1, 1] + enh_onco_td1_df[1, 2] + enh_onco_td1_df[2, 1] + enh_onco_td1_df[2, 2]))
enh_onco_td1 <- fisher.test(enh_onco_td1_df, alternative = "greater")

pro_onco_td1_df <- data.frame(matrix(c(
    length(which(same_og_ep_td1$promoter %in% og_mark$Orthogroup[og_mark$onco])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep_td1$promoter) & og_mark$onco)),
    length(which(same_og_ep_td1$promoter %in% og_mark$Orthogroup[!og_mark$onco])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep_td1$promoter) & !og_mark$onco))
), 2, 2))
pro_onco_td1_fc <- (pro_onco_td1_df[1, 1] / (pro_onco_td1_df[1, 1] + pro_onco_td1_df[1, 2])) /
    ((pro_onco_td1_df[1, 1] + pro_onco_td1_df[2, 1]) / (pro_onco_td1_df[1, 1] + pro_onco_td1_df[1, 2] + pro_onco_td1_df[2, 1] + pro_onco_td1_df[2, 2]))
pro_onco_td1 <- fisher.test(pro_onco_td1_df, alternative = "greater")

enh_hk_td10_df <- data.frame(matrix(c(
    length(which(same_og_ep$enhancer %in% og_mark$Orthogroup[og_mark$hk])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep$enhancer) & og_mark$hk)),
    length(which(same_og_ep$enhancer %in% og_mark$Orthogroup[!og_mark$hk])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep$enhancer) & !og_mark$hk))
), 2, 2))
enh_hk_td10_fc <- (enh_hk_td10_df[1, 1] / (enh_hk_td10_df[1, 1] + enh_hk_td10_df[1, 2])) /
    ((enh_hk_td10_df[1, 1] + enh_hk_td10_df[2, 1]) / (enh_hk_td10_df[1, 1] + enh_hk_td10_df[1, 2] + enh_hk_td10_df[2, 1] + enh_hk_td10_df[2, 2]))
enh_hk_td10 <- fisher.test(enh_hk_td10_df, alternative = "greater")

enh_hk_td1_df <- data.frame(matrix(c(
    length(which(same_og_ep_td1$enhancer %in% og_mark$Orthogroup[og_mark$hk])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep_td1$enhancer) & og_mark$hk)),
    length(which(same_og_ep_td1$enhancer %in% og_mark$Orthogroup[!og_mark$hk])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep_td1$enhancer) & !og_mark$hk))
), 2, 2))
enh_hk_td1_fc <- (enh_hk_td1_df[1, 1] / (enh_hk_td1_df[1, 1] + enh_hk_td1_df[1, 2])) /
    ((enh_hk_td1_df[1, 1] + enh_hk_td1_df[2, 1]) / (enh_hk_td1_df[1, 1] + enh_hk_td1_df[1, 2] + enh_hk_td1_df[2, 1] + enh_hk_td1_df[2, 2]))
enh_hk_td1 <- fisher.test(enh_hk_td1_df, alternative = "greater")

pro_hk_td10_df <- data.frame(matrix(c(
    length(which(same_og_ep$promoter %in% og_mark$Orthogroup[og_mark$hk])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep$promoter) & og_mark$hk)),
    length(which(same_og_ep$promoter %in% og_mark$Orthogroup[!og_mark$hk])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep$promoter) & !og_mark$hk))
), 2, 2))
pro_hk_td10_fc <- (pro_hk_td10_df[1, 1] / (pro_hk_td10_df[1, 1] + pro_hk_td10_df[1, 2])) /
    ((pro_hk_td10_df[1, 1] + pro_hk_td10_df[2, 1]) / (pro_hk_td10_df[1, 1] + pro_hk_td10_df[1, 2] + pro_hk_td10_df[2, 1] + pro_hk_td10_df[2, 2]))
pro_hk_td10 <- fisher.test(pro_hk_td10_df, alternative = "greater")

pro_hk_td1_df <- data.frame(matrix(c(
    length(which(same_og_ep_td1$promoter %in% og_mark$Orthogroup[og_mark$hk])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep_td1$promoter) & og_mark$hk)),
    length(which(same_og_ep_td1$promoter %in% og_mark$Orthogroup[!og_mark$hk])),
    length(which(!(og_mark$Orthogroup %in% same_og_ep_td1$promoter) & !og_mark$hk))
), 2, 2))
pro_hk_td1_fc <- (pro_hk_td1_df[1, 1] / (pro_hk_td1_df[1, 1] + pro_hk_td1_df[1, 2])) /
    ((pro_hk_td1_df[1, 1] + pro_hk_td1_df[2, 1]) / (pro_hk_td1_df[1, 1] + pro_hk_td1_df[1, 2] + pro_hk_td1_df[2, 1] + pro_hk_td1_df[2, 2]))
pro_hk_td1 <- fisher.test(pro_hk_td1_df, alternative = "greater")

onco_hk_enrich_df <- data.frame(
    element = c("enhancer", "promoter", "enhancer", "promoter"),
    gene = c("onco", "onco", "hk", "hk"),
    td10_fc = c(enh_onco_td10_fc, pro_onco_td10_fc, enh_hk_td10_fc, pro_hk_td10_fc),
    td1_fc = c(enh_onco_td1_fc, pro_onco_td1_fc, enh_hk_td1_fc, pro_hk_td1_fc),
    td10_p = c(enh_onco_td10$p.value, pro_onco_td10$p.value, enh_hk_td10$p.value, pro_hk_td10$p.value),
    td1_p = c(enh_onco_td1$p.value, pro_onco_td1$p.value, enh_hk_td1$p.value, pro_hk_td1$p.value),
    gene_num_td10 = c(enh_onco_td10_df[1, 1], pro_onco_td10_df[1, 1], enh_hk_td10_df[1, 1], pro_hk_td10_df[1, 1]),
    gene_num_td1 = c(enh_onco_td1_df[1, 1], pro_onco_td1_df[1, 1], enh_hk_td1_df[1, 1], pro_hk_td1_df[1, 1])
)
ggplot(onco_hk_enrich_df, aes(td10_fc, gene, color = td10_p < 0.05, size = log(gene_num_td10))) +
    geom_point() +
    geom_text(aes(label = round(td10_p, 2)), hjust = 1) +
    theme_minimal() +
    labs(title = "Enrichment foldchange and p-value of cancer genes and housekeeping genes", x = "Foldchange", y = "Element and gene") +
    # scale_color_discretec(c()) +
    facet_wrap(element ~ .)
ggplot(onco_hk_enrich_df, aes(td1_fc, gene, color = td1_p < 0.05, size = log(gene_num_td1))) +
    geom_point() +
    geom_text(aes(label = round(td1_p, 2)), hjust = 1) +
    theme_minimal() +
    labs(title = "Enrichment foldchange and p-value of cancer genes and housekeeping genes", x = "Foldchange", y = "Element and gene") +
    # scale_color_gradient(low = "white", high = "#c01d2e") +
    facet_wrap(element ~ .)
egi_enh_ratio_heat_td1 <- ggplot(egi_sp_df_half_td1, aes(sp_s, sp_t, fill = average_ratio_enh)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "#c01d2e") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Average ratio of same OGs in enhancers", x = "Source species", y = "Target species") # +
# geom_text(aes(label=round(average_ratio_enh,2)))
p_enh_td1 <- egi_enh_ratio_heat_td1 %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)
egi_pro_ratio_heat_td1 <- ggplot(egi_sp_df_half_td1, aes(sp_s, sp_t, fill = average_ratio_pro)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "#02263e") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Average ratio of same OGs in promoters", x = "Source species", y = "Target species")
p_pro_td1 <- egi_pro_ratio_heat_td1 %>%
    insert_bottom(tre_x, height = 0.2) %>%
    insert_left(tre_y, width = 0.2)

save.image("/media/Data/zhangz/chip/analysis/summary2/abc_1M/egi_con_assess.RData")


ent1_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ent_ouwie/longevity/ent1_all.csv", header = TRUE)
ent2_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ent_ouwie/longevity/ent2_all.csv", header = TRUE)
ent1_melt <- reshape2::melt(ent1_all, id.vars = "OGID", variable.name = "species", value.name = "ent1")
ent2_melt <- reshape2::melt(ent2_all, id.vars = "OGID", variable.name = "species", value.name = "ent2")
ent1_violin <- ggplot(drop_na(ent1_melt[ent1_melt$ent1 != 0, ]), aes(x = species)) +
    geom_violin(aes(y = ent1)) +
    geom_boxplot(width = 0.1) +
    theme_minimal() +
    labs(title = "Entropy of OGs in species", x = "Species", y = "Entropy") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 30)
    )
ent1_tre_violin <- ent1_violin %>%
    insert_bottom(tre_x, height = 0.2)
library(preprocessCore)
ent1_all_qn <- normalize.quantiles(as.matrix(ent1_all[, -1])) %>%
    as.data.frame() %>%
    cbind(ent1_all$OGID, .)
colnames(ent1_all_qn) <- colnames(ent1_all)


## plot the distribution of egi distance
# load egi of a species
sp <- "Mus_musculus"
sp <- "Bos_taurus"
sp <- "Rhinolophus_pusillus"
egi_abc <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_egi_abc.csv"), header = TRUE)
nrow(egi_abc[egi_abc$distance == 0 & egi_abc$consrv > 0, ])
egi_nearest <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_egi_nearest.csv"), header = TRUE)
summary(egi_nearest$distance)
nrow(egi_nearest[egi_nearest$distance == 0 & egi_nearest$consrv > 0, ])
nrow(egi_nearest[egi_nearest$consrv > 0, ])
nrow(egi_nearest[egi_nearest$distance == 0 & egi_nearest$consrv >= 10, ])
nrow(egi_nearest[egi_nearest$consrv >= 10, ])
nrow(egi_abc[egi_abc$consrv > 0, ])
nrow(egi_abc[egi_abc$consrv >= 10, ])
nrow(egi_abc[egi_abc$distance == 0 & egi_abc$consrv >= 10, ])

# sp <- "Mus_musculus"
anno_count <- adply(species, 1, function(sp) {
    overlap_anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE)
    overlap_anno <- overlap_anno[!duplicated(overlap_anno$unique_id), ]
    print(paste(sp, ":", nrow(overlap_anno)))
    overlap_anno$consrv <- as.numeric(overlap_anno$consrv)
    # plot distribution bar plot of distanceToTSS absolute value (log10) of enhancers and promoters
    sp_tmp <- gsub("_", " ", sp)
    dist_bar_p_high_ec <- ggplot(overlap_anno[overlap_anno$consrv >= 10, ], aes(x = log10(abs(distanceToTSS) + 1), fill = element)) +
        geom_histogram(aes(y = after_stat(density)), position = "dodge", bins = 12, na.rm = TRUE) +
        scale_fill_manual(values = c("#ef4343", "#73b8d5")) +
        theme_classic() +
        labs(title = paste0("Highly EC CREs in <i>", sp_tmp, "</i>"), x = "Distance to nearest TSS(log10)", y = "Density", fill = "") +
        theme(
            axis.text.x = element_text(hjust = 1),
            text = element_text(size = 30, family = "Arial"),
            panel.grid = element_blank(),
            plot.title = element_markdown(hjust = 0.5)
        ) +
        scale_x_continuous(breaks = c(0, 3, 6))
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_dist_bar_high_ec.png"), dist_bar_p_high_ec, width = 12, height = 8)
    dist_bar_p_ec <- ggplot(overlap_anno[overlap_anno$consrv >= 1, ], aes(x = log10(abs(distanceToTSS) + 1), fill = element)) +
        geom_histogram(aes(y = after_stat(density)), position = "dodge", bins = 12, na.rm = TRUE) +
        scale_fill_manual(values = c("#ef4343", "#73b8d5")) +
        theme_classic() +
        labs(title = paste0("EC CREs in <i>", sp_tmp, "</i>"), x = "Distance to nearest TSS(log10)", y = "Density", fill = "") +
        theme(
            axis.text.x = element_text(hjust = 1),
            text = element_text(size = 30, family = "Arial"),
            panel.grid = element_blank(),
            plot.title = element_markdown(hjust = 0.5)
        ) +
        scale_x_continuous(breaks = c(0, 3, 6))
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_dist_bar_ec.png"), dist_bar_p_ec, width = 12, height = 8)
    dist_bar_p <- ggplot(overlap_anno, aes(x = log10(abs(distanceToTSS) + 1), fill = element)) +
        geom_histogram(aes(y = after_stat(density)), position = "dodge", bins = 12, na.rm = TRUE) +
        scale_fill_manual(values = c("#ef4343", "#73b8d5")) +
        theme_classic() +
        labs(title = paste0("All CREs in <i>", sp_tmp, "</i>"), x = "Distance to nearest TSS(log10)", y = "Density", fill = "") +
        theme(
            axis.text.x = element_text(hjust = 1),
            text = element_text(size = 30, family = "Arial"),
            panel.grid = element_blank(),
            plot.title = element_markdown(hjust = 0.5)
        ) +
        scale_x_continuous(breaks = c(0, 3, 6))
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_dist_bar.png"), dist_bar_p, width = 12, height = 8)
    anno_type <- c("1M_w_ne", "1M_wo_ne", "abc_w_ne", "abc_wo_ne", "abc_1M_w_ne", "abc_1M_wo_ne", "nearest")

    # nec_enh <- length(overlap_anno$unique_id[overlap_anno$consrv == 0 & overlap_anno$element == "enhancer"])
    # nec_pro <- length(overlap_anno$unique_id[overlap_anno$consrv == 0 & overlap_anno$element == "promoter"])
    # ec_enh <- length(overlap_anno$unique_id[overlap_anno$consrv > 0 & overlap_anno$consrv < 10 & overlap_anno$element == "enhancer"])
    # ec_pro <- length(overlap_anno$unique_id[overlap_anno$consrv > 0 & overlap_anno$consrv < 10 & overlap_anno$element == "promoter"])
    # ec10_enh <- length(overlap_anno$unique_id[overlap_anno$consrv >= 10 & overlap_anno$element == "enhancer"])
    # ec10_pro <- length(overlap_anno$unique_id[overlap_anno$consrv >= 10 & overlap_anno$element == "promoter"])
    # all_enh <- length(overlap_anno$unique_id[overlap_anno$element == "enhancer"])
    # all_pro <- length(overlap_anno$unique_id[overlap_anno$element == "promoter"])
    # res <- data.frame(species = sp, nec = c(nec_enh, nec_pro), ec = c(ec_enh, ec_pro), ec10 = c(ec10_enh, ec10_pro),
    #     all = c(all_enh, all_pro), nec_ratio = c(nec_enh / all_enh, nec_pro / all_pro),
    #     ec_ratio = c(ec_enh / all_enh, ec_pro / all_pro), ec10_ratio = c(ec10_enh / all_enh, ec10_pro / all_pro),
    #     element = c("enhancer", "promoter"))
    tes <- tapply(overlap_anno$note_abc, list(overlap_anno$element, overlap_anno$note_abc), length) %>% data.frame()
    tes$element <- rownames(tes)
    tes$species <- sp
    colnames(tes) <- gsub("X", "", colnames(tes))
    return(tes)
})
consrv_count <- adply(species, 1, function(sp) {
    overlap_anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE)
    overlap_anno <- overlap_anno[!duplicated(overlap_anno$unique_id), ]
    print(paste(sp, ":", nrow(overlap_anno)))
    overlap_anno$consrv <- as.numeric(overlap_anno$consrv)
    # plot distribution bar plot of distanceToTSS absolute value (log10) of enhancers and promoters
    # sp_tmp <- gsub("_", " ", sp)
    # dist_bar_p_high_ec <- ggplot(overlap_anno[overlap_anno$consrv >= 10,], aes(x = log10(abs(distanceToTSS) + 1), fill = element)) +
    #     geom_histogram(aes(y = after_stat(density)), position = "dodge", bins = 12, na.rm = TRUE) +
    #     scale_fill_manual(values = c("#ef4343", "#73b8d5")) +
    #     theme_classic() +
    #     labs(title = paste0("Highly EC CREs in <i>", sp_tmp, "</i>"), x = "Distance to nearest TSS(log10)", y = "Density", fill = "") +
    #     theme(axis.text.x = element_text(hjust = 1),
    #         text = element_text(size = 30, family = "Arial"),
    #         panel.grid = element_blank(),
    #         plot.title = element_markdown(hjust = 0.5)) +
    #     scale_x_continuous(breaks = c(0, 3, 6))
    # ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_dist_bar_high_ec.png"), dist_bar_p_high_ec, width = 12, height = 8)
    # dist_bar_p_ec <- ggplot(overlap_anno[overlap_anno$consrv >= 1,], aes(x = log10(abs(distanceToTSS) + 1), fill = element)) +
    #     geom_histogram(aes(y = after_stat(density)), position = "dodge", bins = 12, na.rm = TRUE) +
    #     scale_fill_manual(values = c("#ef4343", "#73b8d5")) +
    #     theme_classic() +
    #     labs(title = paste0("EC CREs in <i>", sp_tmp, "</i>"), x = "Distance to nearest TSS(log10)", y = "Density", fill = "") +
    #     theme(axis.text.x = element_text(hjust = 1),
    #         text = element_text(size = 30, family = "Arial"),
    #         panel.grid = element_blank(),
    #         plot.title = element_markdown(hjust = 0.5)) +
    #     scale_x_continuous(breaks = c(0, 3, 6))
    # ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_dist_bar_ec.png"), dist_bar_p_ec, width = 12, height = 8)
    # dist_bar_p <- ggplot(overlap_anno, aes(x = log10(abs(distanceToTSS) + 1), fill = element)) +
    #     geom_histogram(aes(y = after_stat(density)), position = "dodge", bins = 12, na.rm = TRUE) +
    #     scale_fill_manual(values = c("#ef4343", "#73b8d5")) +
    #     theme_classic() +
    #     labs(title = paste0("All CREs in <i>", sp_tmp, "</i>"), x = "Distance to nearest TSS(log10)", y = "Density", fill = "") +
    #     theme(axis.text.x = element_text(hjust = 1),
    #         text = element_text(size = 30, family = "Arial"),
    #         panel.grid = element_blank(),
    #         plot.title = element_markdown(hjust = 0.5)) +
    #     scale_x_continuous(breaks = c(0, 3, 6))
    # ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_dist_bar.png"), dist_bar_p, width = 12, height = 8)

    nec_enh <- length(overlap_anno$unique_id[overlap_anno$consrv == 0 & overlap_anno$element == "enhancer"])
    nec_pro <- length(overlap_anno$unique_id[overlap_anno$consrv == 0 & overlap_anno$element == "promoter"])
    ec_enh <- length(overlap_anno$unique_id[overlap_anno$consrv > 0 & overlap_anno$consrv < 10 & overlap_anno$element == "enhancer"])
    ec_pro <- length(overlap_anno$unique_id[overlap_anno$consrv > 0 & overlap_anno$consrv < 10 & overlap_anno$element == "promoter"])
    ec10_enh <- length(overlap_anno$unique_id[overlap_anno$consrv >= 10 & overlap_anno$element == "enhancer"])
    ec10_pro <- length(overlap_anno$unique_id[overlap_anno$consrv >= 10 & overlap_anno$element == "promoter"])
    all_enh <- length(overlap_anno$unique_id[overlap_anno$element == "enhancer"])
    all_pro <- length(overlap_anno$unique_id[overlap_anno$element == "promoter"])
    res <- data.frame(
        species = sp, nec = c(nec_enh, nec_pro), ec = c(ec_enh, ec_pro), ec10 = c(ec10_enh, ec10_pro),
        all = c(all_enh, all_pro), nec_ratio = c(nec_enh / all_enh, nec_pro / all_pro),
        ec_ratio = c(ec_enh / all_enh, ec_pro / all_pro), ec10_ratio = c(ec10_enh / all_enh, ec10_pro / all_pro),
        element = c("enhancer", "promoter")
    )
    return(res)
})
# venn_list_enh <- list(all = overlap_anno$unique_id,
#     ec = overlap_anno$unique_id[overlap_anno$consrv > 0],
#     ec10 = overlap_anno$unique_id[overlap_anno$consrv >= 10])
# enh_venn_p <- venn.diagram(venn_list_enh, category.names = c("All", "Conserved", "Conserved >= 10"),
#     filename = "/media/Data/zhangz/chip/analysis/summary2/abc_1M/Mus_musculus_venn_enh.png",
#     imagetype = "png", resolution = 600, fill = c("#73b8d5", "#0d4c6d", "#033250"),
#     cat.pos = c(0, 0, 0),
#     scaled = TRUE)
# venn_list_pro <- list(all = overlap_anno$unique_id,
#     ec = overlap_anno$unique_id[overlap_anno$consrv > 0],
#     ec10 = overlap_anno$unique_id[overlap_anno$consrv >= 10])
# pro_venn_p <- venn.diagram(venn_list_pro, category.names = c("All", "Conserved", "Conserved >= 10"),
#     filename = "/media/Data/zhangz/chip/analysis/summary2/abc_1M/Mus_musculus_venn_pro.png",
#     imagetype = "png", resolution = 600, fill = c("#ef4343", "#c01d2e", "#c53340"),
#     cat.pos = c(0, 0, 0),
#     scaled = TRUE)
tree <- read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
base_tre <- ggtree(tree, layout = "rectangular") + geom_tiplab(size = 2) + theme_tree2()
svg_df <- data.frame(species = tree$tip.label, svg = paste0("/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/", tree$tip.label, ".svg"))
svg_df <- svg_df[-c(14, 15, 16, 17), ]
for (dd in svg_df$svg) {
    if (!file.exists(dd)) {
        print(dd)
    }
}
phypic <- sapply(tree$tip.label, function(x) phylopic_uid(str_split(x, "_")[[1]][1]))
phypic$species <- tree$tip.label
phypic <- phypic[, c(3, 2, 1)]
base_tre_horizontal <- ggtree(tree, layout = "rectangular") %<+% svg_df +
    coord_flip() +
    geom_tiplab(aes(image = svg), geom = "image", offset = 2, size = 0.03)
base_tre_horizontal2 <- ggtree(tree, layout = "rectangular") %<+% phypic +
    coord_flip() +
    geom_tiplab(aes(image = uid), geom = "phylopic", offset = 2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/base_tre_horizontal.pdf", base_tre_horizontal, width = 12, height = 4)
consrv_count <- consrv_count[, -1]
write.csv(consrv_count, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/consrv_count.csv", row.names = FALSE)
consrv_count <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/consrv_count.csv", header = TRUE)
# box plot of enhancer and promoter ec10 num in consrv_count
consrv_num_box_high <- ggplot(consrv_count, aes(x = element, y = ec10, fill = element)) +
    geom_boxplot() +
    # facet_grid(element~.) +
    theme_classic() +
    theme(
        axis.text.x = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.position = "none"
    ) +
    labs(title = "EC10 number of CREs", x = "", y = "Number", fill = "Element") +
    scale_fill_manual(values = c("#ef4343", "#73b8d5"))
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/consrv_num_box_high.png", consrv_num_box_high, width = 8, height = 12)
consrv_num_box <- ggplot(consrv_count, aes(x = element, y = ec, fill = element)) +
    geom_boxplot() +
    # facet_grid(element~.) +
    theme_classic() +
    theme(
        axis.text.x = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.position = "none"
    ) +
    labs(title = "EC number of CREs", x = "", y = "Number", fill = "Element") +
    scale_fill_manual(values = c("#ef4343", "#73b8d5"))
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/consrv_num_box.png", consrv_num_box, width = 8, height = 12)


consrv_count_melt <- melt(consrv_count[, c(1, 6:9)], id.vars = c("species", "element"), variable.name = "EC", value.name = "ratio")

write.csv(consrv_count_melt, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/consrv_count_melt.csv", row.names = FALSE)
consrv_count <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/consrv_count.csv", header = TRUE)
consrv_count_melt$EC <- factor(consrv_count_melt$EC, levels = c("nec_ratio", "ec_ratio", "ec10_ratio"))
stack_p <- ggplot(consrv_count_melt, aes(x = species, y = ratio, fill = as.factor(EC))) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(element ~ .) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = "EC distribution of CREs", x = "", y = "Ratio", fill = "Epigenetically conserved ratio") +
    scale_fill_manual(values = c("#f6f8fa", "#959da5", "#27292f")) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) #+
# theme(text = element_text(size = 15))
p <- stack_p %>%
    insert_bottom(base_tre_horizontal, height = 0.2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/consrv_count.pdf", p, width = 13, height = 8)

stack_p_v <- ggplot(consrv_count_melt, aes(y = species, x = ratio, fill = as.factor(EC))) +
    geom_barh(stat = "identity", position = "stack") +
    facet_grid(. ~ element) +
    theme_classic() +
    # theme(axis.text.x = element_text(angle = 90)) +
    labs(title = "EC distribution of CREs", y = "", x = "Ratio", fill = "EC") +
    scale_fill_manual(values = c("#f6f8fa", "#959da5", "#27292f")) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    theme(text = element_text(size = 15))
p_v <- stack_p_v %>%
    insert_left(base_tre, width = 0.2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/consrv_count_v.pdf", p_v, width = 10, height = 8)

anno_count <- anno_count[, -1]
anno_count_melt <- melt(anno_count, id.vars = c("species", "element"), variable.name = "type", value.name = "count")
anno_count_melt[is.na(anno_count_melt)] <- 0
anno_count_ratio <- ddply(anno_count_melt, .(species, element), transform, ratio = count / sum(count) * 100)
# write.csv(anno_count_ratio, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/anno_count_ratio.csv", row.names = FALSE)
anno_count_ratio <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/anno_count_ratio.csv", header = TRUE)
anno_count_ratio$type <- factor(anno_count_ratio$type, levels = c("1M_w_ne", "1M_wo_ne", "nearest", "abc_w_ne", "abc_wo_ne", "abc_1M_w_ne", "abc_1M_wo_ne"))
stack_anno_p <- ggplot(anno_count_ratio, aes(x = species, y = ratio / 100, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(element ~ .) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = "Annotation distribution of CREs", x = "", y = "Proportion", fill = "Anno type") +
    scale_fill_manual(
        values = c("#959da5", "#959da5", "#959da5", "#2f363d", "#2f363d", "#2f363d", "#2f363d"),
        labels = c("TS-1M", "TS-1M", "Nearest", "ABC", "ABC", "ABC & TS-1M", "ABC & TS-1M")
    ) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    scale_y_continuous(labels = scales::percent, n.breaks = 3) +
    theme(axis.text.x = element_text(face = "italic", size = 8)) +
    theme(axis.text.y = element_text(size = 8)) +
    theme(legend.text = element_text(size = 8)) +
    theme(legend.title = element_text(size = 9)) +
    theme(axis.title = element_text(size = 9)) +
    # theme(legend.position = "top") +
    theme(plot.title = element_text(size = 9))
anno_p <- stack_anno_p %>% insert_bottom(base_tre_horizontal, height = 0.2)
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/anno_count_ratio2.pdf", anno_p, width = 18, height = 10, units = "cm")
wo_ne_type <- c("1M_wo_ne", "abc_wo_ne", "abc_1M_wo_ne")
wo_ne_ratio <- adply(unique(anno_count_ratio$species), 1, function(sp) {
    data.frame(
        species = sp,
        ratio = c(
            sum(anno_count_ratio$ratio[anno_count_ratio$species == sp & anno_count_ratio$type %in% wo_ne_type & anno_count_ratio$element == "enhancer"]),
            sum(anno_count_ratio$ratio[anno_count_ratio$species == sp & anno_count_ratio$type %in% wo_ne_type & anno_count_ratio$element == "promoter"])
        ),
        element = c("enhancer", "promoter")
    )
})
wo_ne_ratio <- wo_ne_ratio[, -1]
# r$> tapply(wo_ne_ratio$ratio, wo_ne_ratio$element, mean)
# enhancer promoter
# 49.81956 39.74662
tapply(wo_ne_ratio$ratio[wo_ne_ratio$species %in% species_ts], wo_ne_ratio$element[wo_ne_ratio$species %in% species_ts], mean)

# plot distribution of overlap_og_num in overlap_anno
consrv_count <- adply(species, 1, function(sp) {
    overlap_anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE)
    # overlap_anno <- overlap_anno[!duplicated(overlap_anno$unique_id), ]
    # overlap_anno <- overlap_anno[!duplicated(overlap_anno[, c("peak", "chr", "start", "end", "species", "tissue", "element")]), ]
    # dup_row <- which(duplicated(overlap_anno[, c("peak", "chr", "start", "end", "species", "tissue", "element")], fromLast = TRUE))
    test <- overlap_anno %>%
        group_by(
            unique_id, peak, chr, start, end,
            species, tissue, element, ele_pattern, id, ec_id, geneId, distanceToTSS,
            lib_fc, align, consrv, gc, mean_fc, length, annotation, motif_num, mean_phylop,
            qn_fc, anno_gene_1M, anno_tau_1M, anno_pattern_1M, anno_distance_1M,
            og_1M, gene_num_1M, og_nearest, overlap_gene, overlap_distance, overlap_num, note_abc
        ) %>%
        summarise(
            overlap_tau = paste(overlap_tau, collapse = ";"),
            overlap_pattern = paste(overlap_pattern, collapse = ";"),
            overlap_og = paste(overlap_og, collapse = ";")
        )
    write.csv(test, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), row.names = FALSE)
    # tes2 <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csvtest"), header = TRUE)
})

all_info <- read.csv("/media/Data/zhangz/chip/analysis/summary2/sum_all/all_1M.csv", header = TRUE, nrows = 5)
colnames(all_info)
classes <- sapply(all_info, class)

cols <- c("peak", "species", "tissue", "element", "align", "consrv", "chr", "start", "end")
classes[-(which(colnames(all_info) %in% cols))] <- rep("NULL", length(classes) - 9)
all_info <- all_info[, cols]
all_info <- read.csv("/media/Data/zhangz/chip/analysis/summary2/sum_all/all_1M.csv", header = TRUE, colClasses = classes)
## plot heat figure of alignment and conservation
align_count <- tapply(all_info$align, list(all_info$align, all_info$element), length) %>% data.frame()
align_count <- data.frame()
for (ele in c("enhancer", "promoter")) {
    for (align in 1:24) {
        for (con in 1:align) {
            align_count <- rbind(
                align_count,
                c(align, con, nrow(all_info[all_info$align == align & all_info$consrv == con & all_info$element == ele, ]), ele)
            )
        }
    }
}

colnames(align_count) <- c("align", "consrv", "count", "element")
align_count$count <- as.numeric(align_count$count)
align_count$align <- factor(align_count$align, levels = 1:24)
align_count$consrv <- factor(align_count$consrv, levels = 1:24)
# classInt::classIntervals(align_count$count[align_count$element == "enhancer"], n = 6, style = "jenks", warnLargeN = FALSE)
# style: jenks
#      [0,5180]  (5180,14414] (14414,27061] (27061,42503] (42503,66464] (66464,99099]
#           164            57            41            25            12             1
# classInt::classIntervals(align_count$count[align_count$element == "promoter"], n = 6, style = "jenks", warnLargeN = FALSE)
# style: jenks
#      [0,2151]   (2151,4629]   (4629,7343]  (7343,13344] (13344,33500] (33500,62810]
#            86           105            78            28             2             1

heat_consrv <- ggplot() +
    geom_tile(data = align_count[align_count$element == "enhancer", ], aes(x = align, y = consrv, fill = count)) +
    theme_bw() +
    scale_fill_steps(low = "#F5DEE4", high = "#c53340", breaks = c(1, 5180, 14414, 27061, 42503, 66464, 99099)) +
    scale_x_discrete(breaks = c(5, 10, 15, 20)) +
    scale_y_discrete(breaks = c(5, 10, 15, 20)) +
    theme(legend.position = "right", legend.text = element_text(size = 10)) +
    new_scale_fill() +
    geom_tile(data = align_count[align_count$element == "promoter", ], aes(x = consrv, y = align, fill = count)) +
    scale_fill_steps(low = "#E2EDF5", high = "#033250", breaks = c(1, 2151, 4629, 7343, 13344, 33500, 62810)) +
    theme_bw() +
    labs(title = "Alignment and conservation of CREs", x = "Alignment", y = "Conservation", fill = "Count") +
    scale_x_discrete(breaks = c(5, 10, 15, 20)) +
    scale_y_discrete(breaks = c(5, 10, 15, 20)) +
    theme(legend.position = "left", legend.text = element_text(size = 10)) +
    theme(text = element_text(size = 25, family = "Arial"))
ggsave("/media/Data/zhangz/chip/analysis/summary2/sum_all/align_consrv.png", heat_consrv, width = 12, height = 8)
