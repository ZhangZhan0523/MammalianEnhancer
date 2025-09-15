library(optparse)

option_list <- list(
    make_option(c("-s", "--species"),
        type = "character", default = "Mus_musculus",
        help = "The species of the data"
    )
)
opt_parser <- parse_args(OptionParser(option_list = option_list, usage = "%prog [options]"))
sp <- opt_parser$species
print(sp)

## maybe load some packages
library(ggplot2)
library(nortest)
library(PerformanceAnalytics)
library(ggpubr)
library(forcats)
library(stringr)

## read in the data
species <- c(
    "Ovis_aries", "Bos_taurus",
    "Sus_scrofa", "Lama_glama", "Mustela_putorius",
    "Canis_lupus", "Felis_catus", "Equus_asinus",
    "Equus_caballus", "Rhinolophus_pusillus", "Rhinolophus_ferrumequinum",
    "Hipposideros_larvatus", "Myotis_ricketti", "Myotis_chinensis",
    "Atelerix_albiventris", "Mus_musculus", "Rattus_norvegicus",
    "Cavia_porcellus", "Oryctolagus_cuniculus", "Macaca_mulatta",
    "Tupaia_belangeri",
    "Procavia_capensis", "Petaurus_breviceps"
)
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/ortho_me_adj.RData")
colnames(ortho_me_adj)[which(colnames(ortho_me_adj) == "Macaca_fascicularis")] <- "Macaca_mulatta"
# old version, including multi copy genes and jiangtun
ortho_old <- read.table("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/ortho_old.tsv", sep = "\t", header = TRUE)
colnames(ortho_old)[which(colnames(ortho_old) == "Macaca_fascicularis")] <- "Macaca_mulatta"

## read in anno stats and concat
# res_stats <- data.frame()
# for (s in species) {
#     for (tis in c("Brain", "Kidney", "Liver")) {
#         for (ele in c("enhancer", "promoter")) {
#             tmp <- try(read.table(paste0(
#                 "/media/Data/zhangz/chip/analysis/", s,
#                 "/anno/peaks2/anno2/", s, "_", tis, "_", ele, "_flank_anno_stats.csv"
#             ), sep = ",", header = TRUE))
#             if (class(tmp) == "try-error") {
#                 next
#             } else {
#                 tmp$species <- s
#                 tmp$tissue <- tis
#                 tmp$element <- ele
#                 res_stats <- rbind(res_stats, tmp)
#             }
#         }
#     }
# }
# res_stats2 <- res_stats[res_stats$anno_total != 0, ]
# write.csv(res_stats2, file = "/media/Data/zhangz/chip/analysis/summary2/anno_sum2/res_stats.csv", row.names = FALSE)

# for (sp in species) {
# sp <- "Mus_musculus"
# ele <- "enhancer"
sp_all <- data.frame()
for (tis in c("Brain", "Kidney", "Liver")) {
    for (ele in c("enhancer", "promoter")) {
        # annotation
        tmp <- read.csv(paste0(
            "/media/Data/zhangz/chip/analysis/", sp, "/anno/peaks2/anno2/", sp, "_", tis, "_", ele, "_flank_anno.csv"
        ), header = TRUE)
        # other infomation
        tmp2 <- read.csv(paste0(
            "/media/Data/zhangz/chip/analysis/", sp, "/compare2/",
            sp, "_", tis, "_", ele, "_info.csv"
        ), header = TRUE, row.names = NULL)
        tmp$tissue <- tis
        tmp$element <- ele
        colnames(tmp)[1] <- "peak"
        tmp3 <- merge(
            tmp[, c(
                "peak", "flank_geneIds", "flank_gene_distances", "ele_pattern", "flank_tau", "flank_pattern",
                "anno_gene_l1", "anno_gene_l2", "l1_num", "l2_num", "note", "anno_tau_l1", "anno_tau_l2",
                "anno_pattern_l1", "anno_pattern_l2",
                "anno_distance_l1", "anno_distance_l2"
            )],
            tmp2,
            by = "peak"
        )
        tmp3$anno_gene_l1 <- gsub("NULL", "", tmp3$anno_gene_l1)
        tmp3$anno_gene_l2 <- gsub("NULL", "", tmp3$anno_gene_l2)
        tmp3$anno_tau_l1 <- gsub("NULL", "", tmp3$anno_tau_l1)
        tmp3$anno_tau_l2 <- gsub("NULL", "", tmp3$anno_tau_l2)
        tmp3$anno_pattern_l1 <- gsub("NULL", "", tmp3$anno_pattern_l1)
        tmp3$anno_pattern_l2 <- gsub("NULL", "", tmp3$anno_pattern_l2)
        tmp3$anno_distance_l1 <- gsub("NULL", "", tmp3$anno_distance_l1)
        tmp3$anno_distance_l2 <- gsub("NULL", "", tmp3$anno_distance_l2)
        tmp3$anno_gene <- paste(tmp3$anno_gene_l1, tmp3$anno_gene_l2, sep = ";")
        tmp3$anno_gene <- gsub("^;", "", tmp3$anno_gene)
        tmp3$anno_gene <- gsub(";$", "", tmp3$anno_gene)
        tmp3$anno_tau <- paste(tmp3$anno_tau_l1, tmp3$anno_tau_l2, sep = ";")
        tmp3$anno_tau <- gsub("^;", "", tmp3$anno_tau)
        tmp3$anno_tau <- gsub(";$", "", tmp3$anno_tau)
        tmp3$anno_pattern <- paste(tmp3$anno_pattern_l1, tmp3$anno_pattern_l2, sep = ";")
        tmp3$anno_pattern <- gsub("^;", "", tmp3$anno_pattern)
        tmp3$anno_pattern <- gsub(";$", "", tmp3$anno_pattern)
        tmp3$anno_distance <- paste(tmp3$anno_distance_l1, tmp3$anno_distance_l2, sep = ";")
        tmp3$anno_distance <- gsub("^;", "", tmp3$anno_distance)
        tmp3$anno_distance <- gsub(";$", "", tmp3$anno_distance)
        tmp3$anno_gene_num <- ifelse(tmp3$anno_gene != "", str_split(tmp3$anno_gene, ";") |> sapply(length), 0)
        tmp3$ec_id <- ""
        for (i in 1:nrow(tmp3)) {
            if (tmp3$anno_gene_num[i] == 0) {
                tmp3$anno_gene[i] <- tmp3$geneId[i]
                tmp3$anno_orthogroup[i] <- ifelse(TRUE %in% grepl(paste0("^", tmp3$geneId[i], "[|]"), ortho_me_adj[, sp]),
                    ortho_me_adj[grep(paste0("^", tmp3$geneId[i], "[|]"), ortho_me_adj[, sp]), "id"],
                    ""
                )
                if (sp == "Tupaia_belangeri") {
                    ## orthogroup file does not contain Tupaia gene names like "TSDBGID"
                    ## so use transcript ID "TSDBTID" as proxies
                    gene_tmp <- gsub("TSDBGID.", "TSDBTID.", tmp3$geneId[i])
                    tmp3$orthogroup_old[i] <- ifelse(TRUE %in% grepl(
                        regex(paste0(";", gene_tmp)),
                        ortho_old[, sp]
                    ),
                    ortho_old[grep(regex(paste0(";", gene_tmp, "[.]")),
                        ortho_old[, sp],
                        perl = TRUE
                    ), "Orthogroup"],
                    ""
                    )
                } else {
                    tmp3$orthogroup_old[i] <- ifelse(TRUE %in% grepl(
                        regex(paste0("(^| )", tmp3$geneId[i], ";")),
                        ortho_old[, sp]
                    ),
                    ortho_old[grep(regex(paste0("(^| )", tmp3$geneId[i], ";")),
                        ortho_old[, sp],
                        perl = TRUE
                    ), "Orthogroup"],
                    ""
                    )
                }
                tmp3$anno_distance[i] <- tmp3$distanceToTSS[i]
                tmp3$anno_gene_num[i] <- 1
                tmp3$note2[i] <- "nearest"
                next
            } else {
                gs <- str_split(tmp3$anno_gene[i], ";")[[1]]
                if (sp == "Tupaia_belangeri") {
                    gs_tmp <- sapply(gs, function(x) gsub("TSDBGID.", "TSDBTID.", x))
                    os_old <- sapply(
                        gs_tmp,
                        function(x) ortho_old[grep(regex(paste0(";", x, "[.]")), ortho_old[, sp], perl = TRUE), "Orthogroup"]
                    )
                    tmp3$orthogroup_old[i] <- paste(os_old, collapse = ";")
                } else {
                    os_old <- sapply(
                        gs,
                        function(x) ortho_old[grep(regex(paste0("(^| )", x, ";")), ortho_old[, sp], perl = TRUE), "Orthogroup"]
                    )
                    tmp3$orthogroup_old[i] <- paste(os_old, collapse = ";")
                }
                os <- sapply(gs, function(x) ortho_me_adj[grep(paste0("^", x, "[|]"), ortho_me_adj[, sp]), "id"])
                tmp3$anno_orthogroup[i] <- paste(os, collapse = ";")
                tmp3$note2[i] <- "ts"
            }
        }
        # epigenetically conserved groups
        ec_group <- read.csv(paste0(
            "/media/Data/zhangz/chip/analysis/", sp, "/compare2/",
            sp, "_", tis, "_", ele, "_func_group.tsv"
        ), header = TRUE, sep = "\t")
        ec_group$count <- apply(ec_group[, 2:ncol(ec_group)], 1, function(x) sum(x != ""))
        ec_group$ec_id <- paste("ec", tolower(substr(tis, 1, 1)), 1:nrow(ec_group), sep = "_")
        write.csv(ec_group,
            paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", tis, "_", ele, "_ec_group.csv"),
            row.names = FALSE
        )
        tmp3$ec_id <- sapply(
            tmp3$peak,
            function(x) ifelse(x %in% ec_group[, sp], ec_group[which(ec_group[, sp] == x), "ec_id"], "")
        )
        # tmp3$anno_gene_num2 <- tmp3$l1_num + tmp3$l2_num
        phylop_tmp <- read.csv(paste0(
            "/media/Data/zhangz/chip/analysis/", sp, "/phyloP/",
            sp, "_", tis, "_", ele, "_phylop.tab"
        ), header = FALSE, sep = "\t")
        colnames(phylop_tmp) <- c("peak", "size", "covered", "sum_phylop", "mean0_phylop", "mean_phylop")
        tmp4 <- merge(tmp3, phylop_tmp, by = "peak")
        sp_all <- rbind(sp_all, tmp4)
    }
}
pattern2tau <- function(pattern) {
    if (pattern %in% c("B", "K", "L")) {
        return(1)
    } else if (pattern %in% c("BK", "BL", "KL")) {
        return(0.5)
    } else if (pattern == "BKL") {
        return(0)
    }
}
sp_all$ele_tau <- sapply(sp_all$ele_pattern, pattern2tau)

# for (i in 1:ncol(tmp3)) {
#     print(class(tmp3[, i]))
# }
## calc the number of enhancers and promoters annnotated to a certain gene,
## take the ratio of enhancer/promoter as the regulation complexity index of this gene
## remove the redundancy of elements according to ele_pattern
## the same enhancer or promoter annotated to different genes in different tissues
## need pair information of elements

bk_en_name <- read.csv(paste0(
    "/media/usb1/chip/analysis/", sp,
    "/anno/peaks/common/", sp, "_Brain_Kidney_enhancer_only.bed"
), header = TRUE, sep = "\t")
bl_en_name <- read.csv(paste0(
    "/media/usb1/chip/analysis/", sp,
    "/anno/peaks/common/", sp, "_Brain_Liver_enhancer_only.bed"
), header = TRUE, sep = "\t")
kl_en_name <- read.csv(paste0(
    "/media/usb1/chip/analysis/", sp,
    "/anno/peaks/common/", sp, "_Kidney_Liver_enhancer_only.bed"
), header = TRUE, sep = "\t")
bkl_en_name <- read.csv(paste0(
    "/media/usb1/chip/analysis/", sp,
    "/anno/peaks/common/", sp, "_Brain_Kidney_Liver_enhancer_common.bed"
), header = TRUE, sep = "\t")
bk_pro_name <- read.csv(paste0(
    "/media/usb1/chip/analysis/", sp,
    "/anno/peaks/common/", sp, "_Brain_Kidney_promoter_only.bed"
), header = TRUE, sep = "\t")
bl_pro_name <- read.csv(paste0(
    "/media/usb1/chip/analysis/", sp,
    "/anno/peaks/common/", sp, "_Brain_Liver_promoter_only.bed"
), header = TRUE, sep = "\t")
kl_pro_name <- read.csv(paste0(
    "/media/usb1/chip/analysis/", sp,
    "/anno/peaks/common/", sp, "_Kidney_Liver_promoter_only.bed"
), header = TRUE, sep = "\t")
bkl_pro_name <- read.csv(paste0(
    "/media/usb1/chip/analysis/", sp,
    "/anno/peaks/common/", sp, "_Brain_Kidney_Liver_promoter_common.bed"
), header = TRUE, sep = "\t")
bk_en_name$id <- paste0("enh_bk_", 1:nrow(bk_en_name))
bl_en_name$id <- paste0("enh_bl_", 1:nrow(bl_en_name))
kl_en_name$id <- paste0("enh_kl_", 1:nrow(kl_en_name))
bkl_en_num <- 1
for (i in 1:nrow(bkl_en_name)) {
    if (i == 1) {
        bkl_en_name$id[i] <- paste0("enh_bkl_", bkl_en_num)
        bkl_en_num <- bkl_en_num + 1
    }
    if (bkl_en_name$Brain_peak[i] %in% bkl_en_name$Brain_peak[1:(i - 1)]) {
        bkl_en_name$id[i] <- bkl_en_name$id[which(bkl_en_name$Brain_peak == bkl_en_name$Brain_peak[i])][1]
    } else if (bkl_en_name$Kidney_peak[i] %in% bkl_en_name$Kidney_peak[1:(i - 1)]) {
        bkl_en_name$id[i] <- bkl_en_name$id[which(bkl_en_name$Kidney_peak == bkl_en_name$Kidney_peak[i])][1]
    } else if (bkl_en_name$Liver_peak[i] %in% bkl_en_name$Liver_peak[1:(i - 1)]) {
        bkl_en_name$id[i] <- bkl_en_name$id[which(bkl_en_name$Liver_peak == bkl_en_name$Liver_peak[i])][1]
    } else {
        bkl_en_name$id[i] <- paste0("enh_bkl_", bkl_en_num)
        bkl_en_num <- bkl_en_num + 1
    }
}
if (nrow(bk_pro_name) != 0) {
    bk_pro_name$id <- paste0("pro_bk_", seq_len(nrow(bk_pro_name)))
}
# bkl_en_name$id <- paste0("en_bkl_", 1:nrow(bkl_en_name))
if (nrow(bl_pro_name) != 0) {
    bl_pro_name$id <- paste0("pro_bl_", seq_len(nrow(bl_pro_name)))
}
if (nrow(kl_pro_name) != 0) {
    kl_pro_name$id <- paste0("pro_kl_", seq_len(nrow(kl_pro_name)))
}
bkl_pro_num <- 1
for (i in seq_len(nrow(bkl_pro_name))) {
    if (i == 1) {
        bkl_pro_name$id[i] <- paste0("pro_bkl_", bkl_pro_num)
        bkl_pro_num <- bkl_pro_num + 1
    }
    if (bkl_pro_name$Brain_peak[i] %in% bkl_pro_name$Brain_peak[1:(i - 1)]) {
        bkl_pro_name$id[i] <- bkl_pro_name$id[which(bkl_pro_name$Brain_peak == bkl_pro_name$Brain_peak[i])][1]
    } else if (bkl_pro_name$Kidney_peak[i] %in% bkl_pro_name$Kidney_peak[1:(i - 1)]) {
        bkl_pro_name$id[i] <- bkl_pro_name$id[which(bkl_pro_name$Kidney_peak == bkl_pro_name$Kidney_peak[i])][1]
    } else if (bkl_pro_name$Liver_peak[i] %in% bkl_pro_name$Liver_peak[1:(i - 1)]) {
        bkl_pro_name$id[i] <- bkl_pro_name$id[which(bkl_pro_name$Liver_peak == bkl_pro_name$Liver_peak[i])][1]
    } else {
        bkl_pro_name$id[i] <- paste0("pro_bkl_", bkl_pro_num)
        bkl_pro_num <- bkl_pro_num + 1
    }
}

enh_sp_num <- 1
pro_sp_num <- 1

for (i in 1:nrow(sp_all)) {
    if (sp_all$element[i] == "enhancer") {
        if (sp_all$ele_tau[i] == 1) {
            sp_all$id[i] <- paste0("enh_", tolower(sp_all$ele_pattern[i]), "_", enh_sp_num)
            enh_sp_num <- enh_sp_num + 1
        } else if (sp_all$ele_pattern[i] == "BK") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- bk_en_name$id[which(bk_en_name$Brain_peak == sp_all$peak[i])]
            } else if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- bk_en_name$id[which(bk_en_name$Kidney_peak == sp_all$peak[i])]
            }
        } else if (sp_all$ele_pattern[i] == "BL") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- bl_en_name$id[which(bl_en_name$Brain_peak == sp_all$peak[i])]
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- bl_en_name$id[which(bl_en_name$Liver_peak == sp_all$peak[i])]
            }
        } else if (sp_all$ele_pattern[i] == "KL") {
            if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- kl_en_name$id[which(kl_en_name$Kidney_peak == sp_all$peak[i])]
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- kl_en_name$id[which(kl_en_name$Liver_peak == sp_all$peak[i])]
            }
        } else if (sp_all$ele_pattern[i] == "BKL") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- unique(bkl_en_name$id[which(bkl_en_name$Brain_peak == sp_all$peak[i])])
            } else if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- unique(bkl_en_name$id[which(bkl_en_name$Kidney_peak == sp_all$peak[i])])
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- unique(bkl_en_name$id[which(bkl_en_name$Liver_peak == sp_all$peak[i])])
            }
        }
    } else if (sp_all$element[i] == "promoter") {
        if (sp_all$ele_tau[i] == 1) {
            sp_all$id[i] <- paste0("pro_", tolower(sp_all$ele_pattern[i]), "_", pro_sp_num)
            pro_sp_num <- pro_sp_num + 1
        } else if (sp_all$ele_pattern[i] == "BK") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- bk_pro_name$id[which(bk_pro_name$Brain_peak == sp_all$peak[i])]
            } else if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- bk_pro_name$id[which(bk_pro_name$Kidney_peak == sp_all$peak[i])]
            }
        } else if (sp_all$ele_pattern[i] == "BL") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- bl_pro_name$id[which(bl_pro_name$Brain_peak == sp_all$peak[i])]
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- bl_pro_name$id[which(bl_pro_name$Liver_peak == sp_all$peak[i])]
            }
        } else if (sp_all$ele_pattern[i] == "KL") {
            if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- kl_pro_name$id[which(kl_pro_name$Kidney_peak == sp_all$peak[i])]
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- kl_pro_name$id[which(kl_pro_name$Liver_peak == sp_all$peak[i])]
            }
        } else if (sp_all$ele_pattern[i] == "BKL") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- unique(bkl_pro_name$id[which(bkl_pro_name$Brain_peak == sp_all$peak[i])])
            } else if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- unique(bkl_pro_name$id[which(bkl_pro_name$Kidney_peak == sp_all$peak[i])])
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- unique(bkl_pro_name$id[which(bkl_pro_name$Liver_peak == sp_all$peak[i])])
            }
        }
    }
}

write.csv(sp_all,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_all.csv"),
    row.names = FALSE
)
write.csv(bk_en_name,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_bk_en_name.csv"),
    row.names = FALSE
)
write.csv(bl_en_name,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_bl_en_name.csv"),
    row.names = FALSE
)
write.csv(kl_en_name,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_kl_en_name.csv"),
    row.names = FALSE
)
write.csv(bkl_en_name,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_bkl_en_name.csv"),
    row.names = FALSE
)
write.csv(bk_pro_name,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_bk_pro_name.csv"),
    row.names = FALSE
)
write.csv(bl_pro_name,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_bl_pro_name.csv"),
    row.names = FALSE
)
write.csv(kl_pro_name,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_kl_pro_name.csv"),
    row.names = FALSE
)
write.csv(bkl_pro_name,
    paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_bkl_pro_name.csv"),
    row.names = FALSE
)
