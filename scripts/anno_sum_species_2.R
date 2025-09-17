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
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/ortho_me.RData")
colnames(ortho_me)[which(colnames(ortho_me) == "Macaca_fascicularis")] <- "Macaca_mulatta"

## read in anno stats and concat
res_stats <- data.frame()
for (s in species) {
    for (tis in c("Brain", "Kidney", "Liver")) {
        for (ele in c("enhancer", "promoter")) {
            tmp <- try(read.table(paste0(
                "/media/Data/zhangz/chip/analysis/", s,
                "/anno/peaks2/anno2/", s, "_", tis, "_", ele, "_flank_anno_stats.csv"
            ), sep = ",", header = TRUE))
            if (class(tmp) == "try-error") {
                next
            } else {
                tmp$species <- s
                tmp$tissue <- tis
                tmp$element <- ele
                res_stats <- rbind(res_stats, tmp)
            }
        }
    }
}
res_stats2 <- res_stats[res_stats$anno_total != 0, ]
write.csv(res_stats2, file = "/media/Data/zhangz/chip/analysis/summary2/anno_sum2/res_stats.csv", row.names = FALSE)

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
                tmp3$anno_orthogroup[i] <- ""
                next
            } else {
                gs <- str_split(tmp3$anno_gene[i], ";")[[1]]
                os <- sapply(gs, function(x) ortho_me[grep(paste0("^", x, "[|]"), ortho_me[, sp]), "id"])
                tmp3$anno_orthogroup[i] <- paste(os, collapse = ";")
            }
        }
        # epigenetically conserved groups
        ec_group <- read.csv(paste0(
            "/media/Data/zhangz/chip/analysis/", sp, "/compare2/",
            sp, "_", tis, "_", ele, "_func_group.tsv"
        ), header = TRUE, sep = "\t")
        ec_group$count <- apply(ec_group[, 2:ncol(ec_group)], 1, function(x) sum(x != ""))
        ec_group$ec_id <- paste("ec", tolower(substr(tis, 1, 1)), 1:nrow(ec_group), sep = "_")
        write.csv(ec_group, paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", tis, "_", ele, "_ec_group.csv"), row.names = FALSE)
        tmp3$ec_id <- sapply(tmp3$peak, function(x) ifelse(x %in% ec_group[, sp], ec_group[which(ec_group[, sp] == x), "ec_id"], ""))
        # tmp3$anno_gene_num2 <- tmp3$l1_num + tmp3$l2_num
        sp_all <- rbind(sp_all, tmp3)
    }
}

write.csv(sp_all, paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_all.csv"), row.names = FALSE)

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
                sp_all$id[i] <- bk_en_name$id[which(bk_en_name$Brain_peak == sp_all$V4[i])]
            } else if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- bk_en_name$id[which(bk_en_name$Kidney_peak == sp_all$V4[i])]
            }
        } else if (sp_all$ele_pattern[i] == "BL") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- bl_en_name$id[which(bl_en_name$Brain_peak == sp_all$V4[i])]
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- bl_en_name$id[which(bl_en_name$Liver_peak == sp_all$V4[i])]
            }
        } else if (sp_all$ele_pattern[i] == "KL") {
            if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- kl_en_name$id[which(kl_en_name$Kidney_peak == sp_all$V4[i])]
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- kl_en_name$id[which(kl_en_name$Liver_peak == sp_all$V4[i])]
            }
        } else if (sp_all$ele_pattern[i] == "BKL") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- unique(bkl_en_name$id[which(bkl_en_name$Brain_peak == sp_all$V4[i])])
            } else if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- unique(bkl_en_name$id[which(bkl_en_name$Kidney_peak == sp_all$V4[i])])
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- unique(bkl_en_name$id[which(bkl_en_name$Liver_peak == sp_all$V4[i])])
            }
        }
    } else if (sp_all$element[i] == "promoter") {
        if (sp_all$ele_tau[i] == 1) {
            sp_all$id[i] <- paste0("pro_", tolower(sp_all$ele_pattern[i]), "_", pro_sp_num)
            pro_sp_num <- pro_sp_num + 1
        } else if (sp_all$ele_pattern[i] == "BK") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- bk_pro_name$id[which(bk_pro_name$Brain_peak == sp_all$V4[i])]
            } else if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- bk_pro_name$id[which(bk_pro_name$Kidney_peak == sp_all$V4[i])]
            }
        } else if (sp_all$ele_pattern[i] == "BL") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- bl_pro_name$id[which(bl_pro_name$Brain_peak == sp_all$V4[i])]
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- bl_pro_name$id[which(bl_pro_name$Liver_peak == sp_all$V4[i])]
            }
        } else if (sp_all$ele_pattern[i] == "KL") {
            if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- kl_pro_name$id[which(kl_pro_name$Kidney_peak == sp_all$V4[i])]
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- kl_pro_name$id[which(kl_pro_name$Liver_peak == sp_all$V4[i])]
            }
        } else if (sp_all$ele_pattern[i] == "BKL") {
            if (sp_all$tissue[i] == "Brain") {
                sp_all$id[i] <- unique(bkl_pro_name$id[which(bkl_pro_name$Brain_peak == sp_all$V4[i])])
            } else if (sp_all$tissue[i] == "Kidney") {
                sp_all$id[i] <- unique(bkl_pro_name$id[which(bkl_pro_name$Kidney_peak == sp_all$V4[i])])
            } else if (sp_all$tissue[i] == "Liver") {
                sp_all$id[i] <- unique(bkl_pro_name$id[which(bkl_pro_name$Liver_peak == sp_all$V4[i])])
            }
        }
    }
}


for (ele in c("enhancer", "promoter")) {
    brain_ele_anno <- read.csv(paste0(
        "/media/Data/zhangz/chip/analysis/", sp, "/anno/peaks2/anno2/",
        sp, "_Brain_", ele, "_flank_anno.csv"
    ), header = TRUE)
    kidney_ele_anno <- read.csv(paste0(
        "/media/Data/zhangz/chip/analysis/", sp, "/anno/peaks2/anno2/",
        sp, "_Kidney_", ele, "_flank_anno.csv"
    ), header = TRUE)
    liver_ele_anno <- read.csv(paste0(
        "/media/Data/zhangz/chip/analysis/", sp, "/anno/peaks2/anno2/",
        sp, "_Liver_", ele, "_flank_anno.csv"
    ), header = TRUE)
    bk <- intersect(brain_ele_anno$anno_gene, kidney_ele_anno$anno_gene)
    bk <- bk[-which(bk == "NULL")]
    kl <- intersect(kidney_ele_anno$anno_gene, liver_ele_anno$anno_gene)
    kl <- kl[-which(kl == "NULL")]
    bl <- intersect(brain_ele_anno$anno_gene, liver_ele_anno$anno_gene)
    bl <- bl[-which(bl == "NULL")]
    bkl <- intersect(intersect(brain_ele_anno$anno_gene, kidney_ele_anno$anno_gene), liver_ele_anno$anno_gene)
    bkl <- bkl[-which(bkl == "NULL")]
    # dim(brain_ele_anno[brain_ele_anno$anno_gene == bk[2], ])
    # dim(kidney_ele_anno[kidney_ele_anno$anno_gene == bk[2], ])
    # dim(brain_ele_anno[brain_ele_anno$anno_gene == kl[2], ])
    # dim(brain_ele_anno[brain_ele_anno$anno_gene %in% bk, ])
    # dim(kidney_ele_anno[kidney_ele_anno$anno_gene %in% bk, ])
    # dim(liver_ele_anno[liver_ele_anno$anno_gene %in% bk, ])

    brain_ele_gc <- read.csv(paste0(
        "/media/usb1/chip/analysis/", sp, "/gc/",
        sp, "_Brain_", ele, "_gc.tsv"
    ), header = TRUE, sep = "\t")
    colnames(brain_ele_gc) <- c("seqnames", "start", "end", "V4", "at", "gc", "a", "c", "g", "t", "n", "oth", "seq_len")
    brain_ele <- merge(brain_ele_anno, brain_ele_gc, by = "V4")
    kidney_ele_gc <- read.csv(paste0(
        "/media/usb1/chip/analysis/", sp, "/gc/",
        sp, "_Kidney_", ele, "_gc.tsv"
    ), header = TRUE, sep = "\t")
    colnames(kidney_ele_gc) <- c("seqnames", "start", "end", "V4", "at", "gc", "a", "c", "g", "t", "n", "oth", "seq_len")
    kidney_ele <- merge(kidney_ele_anno, kidney_ele_gc, by = "V4")
    liver_ele_gc <- read.csv(paste0(
        "/media/usb1/chip/analysis/", sp, "/gc/",
        sp, "_Liver_", ele, "_gc.tsv"
    ), header = TRUE, sep = "\t")
    colnames(liver_ele_gc) <- c("seqnames", "start", "end", "V4", "at", "gc", "a", "c", "g", "t", "n", "oth", "seq_len")
    liver_ele <- merge(liver_ele_anno, liver_ele_gc, by = "V4")


    bkl_stats <- data.frame()
    for (i in 1:length(bkl)) {
        b <- brain_ele[brain_ele$anno_gene == bkl[i], ]
        k <- kidney_ele[kidney_ele$anno_gene == bkl[i], ]
        l <- liver_ele[liver_ele$anno_gene == bkl[i], ]
        bkl_stats <- rbind(bkl_stats, data.frame(
            gene = bkl[i], brain_gc = mean(b$gc), kidney_gc = mean(k$gc), liver_gc = mean(l$gc),
            brain_len = mean(b$seq_len), kidney_len = mean(k$seq_len), liver_len = mean(l$seq_len)
        ))
    }
    write.csv(bkl_stats, paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_bkl_stats.csv"), row.names = FALSE)
    # brain_ele_phylop <- read.csv("/media/Data/zhangz/chip/genomes/phyloP/Mus_musculus/musMus_Mus_musculus_Brain_elehancer_optimal.narrowPeak_phyloP.tab",
    #     header = FALSE, sep = "\t"
    # )
    # kidney_ele_phylop <- read.csv("/media/Data/zhangz/chip/genomes/phyloP/Mus_musculus/musMus_Mus_musculus_Kidney_elehancer_optimal.narrowPeak_phyloP.tab",
    #     header = FALSE, sep = "\t"
    # )
    # liver_ele_phylop <- read.csv("/media/Data/zhangz/chip/genomes/phyloP/Mus_musculus/musMus_Mus_musculus_Liver_elehancer_optimal.narrowPeak_phyloP.tab",
    #     header = FALSE, sep = "\t"
    # )
    # colnames(brain_ele_phylop) <- c("V4", "len", "covered", "sum", "mean0", "mean")
    # colnames(kidney_ele_phylop) <- c("V4", "len", "covered", "sum", "mean0", "mean")
    # colnames(liver_ele_phylop) <- c("V4", "len", "covered", "sum", "mean0", "mean")
    # brain_ele <- merge(brain_ele, brain_ele_phylop, by = "V4")
    # kidney_ele <- merge(kidney_ele, kidney_ele_phylop, by = "V4")
    # liver_ele <- merge(liver_ele, liver_ele_phylop, by = "V4")


    brain_ele$tissue <- "Brain"
    kidney_ele$tissue <- "Kidney"
    liver_ele$tissue <- "Liver"
    all_ele <- rbind(brain_ele, kidney_ele, liver_ele)
    all_ele$tissue <- factor(all_ele$tissue, levels = c("Brain", "Kidney", "Liver"))
    aov.gc.all <- aov(gc ~ tissue, data = all_ele)
    summary(aov.gc.all)
    # aov.phylop.all <- aov(mean ~ tissue, data = all_ele)
    # summary(aov.phylop.all)
    # aov.phylop0.all <- aov(mean0 ~ tissue, data = all_ele)
    # summary(aov.phylop0.all)
    aov.len.all <- aov(seq_len ~ tissue, data = all_ele)
    summary(aov.len.all)
    ggplot(all_ele, aes(x = tissue, y = gc)) +
        geom_boxplot()
    # ggplot(all_ele, aes(x = tissue, y = mean)) +
    #     geom_boxplot()
    # ggplot(all_ele, aes(x = tissue, y = mean0)) +
    #     geom_boxplot()
    ggplot(all_ele, aes(x = tissue, y = seq_len)) +
        geom_boxplot()
    bartlett.test(gc ~ tissue, data = all_ele)
    # bartlett.test(mean ~ tissue, data = all_ele)
    # bartlett.test(mean0 ~ tissue, data = all_ele)
    bartlett.test(seq_len ~ tissue, data = all_ele)
    # shapiro.test(sample(all_ele$gc, 5000))
    # shapiro.test(sample(all_ele$mean, 5000))
    # shapiro.test(sample(all_ele$mean0, 5000))
    # shapiro.test(sample(all_ele$seq_len, 5000))
    library(multcomp)
    glht.gc.all <- glht(aov.gc.all, linfct = mcp(tissue = "Tukey")) |> summary()
    # glht(aov.phylop.all, linfct = mcp(tissue = "Tukey")) |> summary()
    # glht(aov.phylop0.all, linfct = mcp(tissue = "Tukey")) |> summary()
    glht(aov.len.all, linfct = mcp(tissue = "Tukey")) |> summary()
    oneway.test(gc ~ tissue, data = all_ele)
    # oneway.test(mean ~ tissue, data = all_ele)
    # oneway.test(mean0 ~ tissue, data = all_ele)
    oneway.test(seq_len ~ tissue, data = all_ele)
    kt.gc.all <- kruskal.test(gc ~ tissue, data = all_ele)
    # kt.phylop.all <- kruskal.test(mean ~ tissue, data = all_ele)
    # kt.phylop0.all <- kruskal.test(mean0 ~ tissue, data = all_ele)
    kt.len.all <- kruskal.test(seq_len ~ tissue, data = all_ele)
    all_wilcox_gc <- with(
        all_ele,
        pairwise.wilcox.test(gc, tissue, pool.sd = FALSE, p.adjust.method = "holm")
    )
    # all_wilcox_phylop <- with(
    #     all_ele,
    #     pairwise.wilcox.test(mean, tissue, pool.sd = FALSE, p.adjust.method = "holm")
    # )
    # all_wilcox_phylop0 <- with(
    #     all_ele,
    #     pairwise.wilcox.test(mean0, tissue, pool.sd = FALSE, p.adjust.method = "holm")
    # )
    all_wilcox_len <- with(
        all_ele,
        pairwise.wilcox.test(seq_len, tissue, pool.sd = FALSE, p.adjust.method = "holm")
    )
    with(
        all_ele,
        pairwise.wilcox.test(gc, tissue, pool.sd = FALSE, p.adjust.method = "bonferroni")
    )
    common_gene_ele <- all_ele[all_ele$anno_gene %in% bkl, ]
    common_gene_aov <- data.frame()
    for (gene in bkl) {
        sel_ele <- all_ele[all_ele$anno_gene == gene, ]
        # p<0.05, not the same variation
        # var_gc <- bartlett.test(gc ~ tissue, data = sel_ele)
        # var_phylop <- bartlett.test(mean ~ tissue, data = sel_ele)
        # var_phylop0 <- bartlett.test(mean0 ~ tissue, data = sel_ele)
        # var_len <- bartlett.test(seq_len ~ tissue, data = sel_ele)
        # # p<0.05, not normal distribution
        # norm_gc <- lillie.test(sel_ele$gc)$p.value
        # norm_phylop <- lillie.test(sel_ele$mean)$p.value
        # norm_phylop0 <- lillie.test(sel_ele$mean0)$p.value
        # norm_len <- lillie.test(sel_ele$seq_len)$p.value
        # aov_gc <- aov(gc ~ tissue, data = sel_ele)
        # aov_phylop <- aov(mean ~ tissue, data = sel_ele)
        # aov_phylop0 <- aov(mean0 ~ tissue, data = sel_ele)
        # aov_len <- aov(seq_len ~ tissue, data = sel_ele)
        # glht_gc <- glht(aov.gc, linfct = mcp(tissue = "Tukey")) |> summary()
        # glht_phylop <- glht(aov.phylop, linfct = mcp(tissue = "Tukey")) |> summary()
        # glht_phylop0 <- glht(aov.phylop0, linfct = mcp(tissue = "Tukey")) |> summary()
        # glht_len <- glht(aov.len, linfct = mcp(tissue = "Tukey")) |> summary()
        # owt_gc <- oneway.test(gc ~ tissue, data = sel_ele)$p.value
        # owt_phylop <- oneway.test(mean ~ tissue, data = sel_ele)$p.value
        # owt_phylop0 <- oneway.test(mean0 ~ tissue, data = sel_ele)$p.value
        # owt_len <- oneway.test(seq_len ~ tissue, data = sel_ele)$p.value
        kt_gc <- kruskal.test(gc ~ tissue, data = sel_ele)$p.value
        # kt_phylop <- kruskal.test(mean ~ tissue, data = sel_ele)$p.value
        # kt_phylop0 <- kruskal.test(mean0 ~ tissue, data = sel_ele)$p.value
        kt_len <- kruskal.test(seq_len ~ tissue, data = sel_ele)$p.value
        wilcox_gc <- with(
            sel_ele,
            pairwise.wilcox.test(gc, tissue, pool.sd = FALSE, p.adjust.method = "holm")
        )
        # wilcox_phylop <- with(
        #     sel_ele,
        #     pairwise.wilcox.test(mean, tissue, pool.sd = FALSE, p.adjust.method = "holm")
        # )
        # wilcox_phylop0 <- with(
        #     sel_ele,
        #     pairwise.wilcox.test(mean0, tissue, pool.sd = FALSE, p.adjust.method = "holm")
        # )
        wilcox_len <- with(
            sel_ele,
            pairwise.wilcox.test(seq_len, tissue, pool.sd = FALSE, p.adjust.method = "holm")
        )
        common_gene_aov <- rbind(common_gene_aov, data.frame(
            gene = gene, total_num = dim(sel_ele)[1], brain_num = dim(sel_ele[sel_ele$tissue == "Brain", ])[1],
            kidney_num = dim(sel_ele[sel_ele$tissue == "Kidney", ])[1], liver_num = dim(sel_ele[sel_ele$tissue == "Liver", ])[1],
            kt_gc = kt_gc, kt_len = kt_len,
            # kt_phylop = kt_phylop, kt_phylop0 = kt_phylop0,
            wilcox_gc_bk = wilcox_gc$p.value["Kidney", "Brain"], wilcox_gc_bl = wilcox_gc$p.value["Liver", "Brain"],
            wilcox_gc_kl = wilcox_gc$p.value["Liver", "Kidney"],
            # wilcox_phylop_bk = wilcox_phylop$p.value["Kidney", "Brain"], wilcox_phylop_bl = wilcox_phylop$p.value["Liver", "Brain"],
            # wilcox_phylop_kl = wilcox_phylop$p.value["Liver", "Kidney"],
            # wilcox_phylop0_bk = wilcox_phylop0$p.value["Kidney", "Brain"], wilcox_phylop0_bl = wilcox_phylop0$p.value["Liver", "Brain"],
            # wilcox_phylop0_kl = wilcox_phylop0$p.value["Liver", "Kidney"],
            wilcox_len_bk = wilcox_len$p.value["Kidney", "Brain"], wilcox_len_bl = wilcox_len$p.value["Liver", "Brain"],
            wilcox_len_kl = wilcox_len$p.value["Liver", "Kidney"]
        ))
    }

    lillie.test(common_gene_ele$gc)
    # lillie.test(common_gene_ele$mean)
    # lillie.test(common_gene_ele$mean0)
    lillie.test(common_gene_ele$seq_len)
    bartlett.test(gc ~ tissue, data = common_gene_ele)
    # bartlett.test(mean ~ tissue, data = common_gene_ele)
    # bartlett.test(mean0 ~ tissue, data = common_gene_ele)
    bartlett.test(seq_len ~ tissue, data = common_gene_ele)
    kruskal.test(gc ~ tissue, data = common_gene_ele)
    # kruskal.test(mean ~ tissue, data = common_gene_ele)
    # kruskal.test(mean0 ~ tissue, data = common_gene_ele)
    kruskal.test(seq_len ~ tissue, data = common_gene_ele)
    oneway.test(gc ~ tissue, data = common_gene_ele)
    # oneway.test(mean ~ tissue, data = common_gene_ele)
    # oneway.test(mean0 ~ tissue, data = common_gene_ele)
    oneway.test(seq_len ~ tissue, data = common_gene_ele)
    common_wilcox_gc <- with(
        common_gene_ele,
        pairwise.wilcox.test(gc, tissue, pool.sd = FALSE, p.adjust.method = "holm")
    )
    # common_wilcox_phylop <- with(
    #     common_gene_ele,
    #     pairwise.wilcox.test(mean, tissue, pool.sd = FALSE, p.adjust.method = "holm")
    # )
    # common_wilcox_phylop0 <- with(
    #     common_gene_ele,
    #     pairwise.wilcox.test(mean0, tissue, pool.sd = FALSE, p.adjust.method = "holm")
    # )
    common_wilcox_len <- with(
        common_gene_ele,
        pairwise.wilcox.test(seq_len, tissue, pool.sd = FALSE, p.adjust.method = "holm")
    )
    with(
        common_gene_ele,
        pairwise.wilcox.test(gc, tissue, pool.sd = FALSE, p.adjust.method = "bonferroni")
    )
    ## plot boxplot of gc, phylop, phylop0, seq_len of selected genes
    png(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_common_gene_aov.png"), width = 12, height = 8, units = "in", res = 300)
    chart.Boxplot(common_gene_aov[, 6:13])
    dev.off()
    my_comarison <- list(c("Brain", "Kidney"), c("Brain", "Liver"), c("Kidney", "Liver"))
    # levels(all_ele$tissue) <- rev(levels(fct_reorder(all_ele$tissue, all_ele$gc, mean)))
    all_gc_boxp <- ggplot(all_ele, aes(x = fct_reorder(tissue, gc, mean), y = gc)) +
        geom_violin(aes(fill = tissue)) +
        geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
        theme_bw() +
        theme() +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times")
        ) +
        xlab("Tissue") +
        ylab("GC content") +
        ggtitle(paste0("GC content of all enhancers in "), sp) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        scale_color_manual(values = c("Brain" = "red", "Kidney" = "blue", "Liver" = "green")) +
        stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_all_gc_boxp.pdf"), all_gc_boxp,
        width = 10, height = 8, units = "in"
    )
    len_ran <- quantile(all_ele$seq_len, c(0.025, 0.975))
    all_ele_no_outlier <- all_ele[which(all_ele$seq_len < len_ran[2] & all_ele$seq_len > len_ran[1]), ]
    all_len_boxp_no_outlier <- ggplot(all_ele_no_outlier, aes(x = fct_reorder(tissue, seq_len, mean), y = seq_len)) +
        geom_violin(aes(fill = tissue), trim = TRUE) +
        geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
        theme_bw() +
        theme() +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times")
        ) +
        # ylim(0, 5000) +
        xlab("Tissue") +
        ylab("Sequence length") +
        ggtitle(paste0("Sequence length of all ", ele, " in "), sp) +
        scale_color_manual(values = c("Brain" = "red", "Kidney" = "blue", "Liver" = "green")) +
        # scale_y_continuous(breaks = c(0, 0.5, 1)) +
        stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    # all_phylop_boxp <- ggplot(all_ele, aes(x = tissue, y = mean)) +
    #     geom_violin(aes(fill = tissue)) +
    #     geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    #     theme_bw() +
    #     theme() +
    #     theme(
    #         axis.text = element_text(size = 20, family = "Times"),
    #         axis.title = element_text(size = 25, family = "Times"),
    #         title = element_text(size = 30, family = "Times"),
    #         legend.title = element_text(size = 20, family = "Times"),
    #         legend.text = element_text(size = 20, family = "Times")
    #     ) +
    #     xlab("Tissue") +
    #     ylab("PhyloP score") +
    #     ggtitle("PhyloP score of all enhancers in three tissues") +
    #     # scale_y_continuous(breaks = c(0, 0.5, 1)) +
    #     stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    # ggsave("/media/Data/zhangz/chip/analysis/summary2/anno_sum/all_phylop_boxp.pdf", all_phylop_boxp,
    #     width = 10, height = 8, units = "in"
    # )
    # all_phylop0_boxp <- ggplot(all_ele, aes(x = tissue, y = mean0)) +
    #     geom_violin(aes(fill = tissue)) +
    #     geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    #     theme_bw() +
    #     theme() +
    #     theme(
    #         axis.text = element_text(size = 20, family = "Times"),
    #         axis.title = element_text(size = 25, family = "Times"),
    #         title = element_text(size = 30, family = "Times"),
    #         legend.title = element_text(size = 20, family = "Times"),
    #         legend.text = element_text(size = 20, family = "Times")
    #     ) +
    #     xlab("Tissue") +
    #     ylab("PhyloP score") +
    #     ggtitle("PhyloP score of all enhancers in three tissues") +
    #     # scale_y_continuous(breaks = c(0, 0.5, 1)) +
    #     stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    # ggsave("/media/Data/zhangz/chip/analysis/summary2/anno_sum/all_phylop0_boxp.pdf", all_phylop0_boxp,
    #     width = 10, height = 8, units = "in"
    # )
    # levels(all_ele$tissue) <- rev(levels(fct_reorder(all_ele$tissue, all_ele$gseq_len, mean)))
    all_len_boxp <- ggplot(all_ele[which(all_ele$seq_len < 5000), ], aes(x = fct_reorder(tissue, seq_len, mean), y = seq_len)) +
        geom_violin(aes(fill = tissue), trim = TRUE) +
        geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
        theme_bw() +
        theme() +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times")
        ) +
        # ylim(0, 5000) +
        xlab("Tissue") +
        ylab("Sequence length") +
        ggtitle(paste0("Sequence length of all enhancers in "), sp) +
        scale_color_manual(values = c("Brain" = "red", "Kidney" = "blue", "Liver" = "green")) +
        # scale_y_continuous(breaks = c(0, 0.5, 1)) +
        stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_all_len_boxp.pdf"), all_len_boxp,
        width = 10, height = 8, units = "in"
    )
    # levels(common_gene_ele$tissue) <- rev(levels(fct_reorder(common_gene_ele$tissue, common_gene_ele$gc, mean)))
    common_gc_boxp <- ggplot(common_gene_ele, aes(x = fct_reorder(tissue, gc, mean), y = gc)) +
        geom_violin(aes(fill = tissue)) +
        geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
        theme_bw() +
        theme() +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times")
        ) +
        xlab("Tissue") +
        ylab("GC content") +
        ggtitle(paste0("GC content of selected genes in "), sp) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        scale_color_manual(values = c("Brain" = "red", "Kidney" = "blue", "Liver" = "green")) +
        stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_common_gc_boxp.pdf"), common_gc_boxp,
        width = 10, height = 8, units = "in"
    )
    # common_phylop_boxp <- ggplot(common_gene_ele, aes(x = tissue, y = mean)) +
    #     geom_violin(aes(fill = tissue)) +
    #     geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    #     theme_bw() +
    #     theme() +
    #     theme(
    #         axis.text = element_text(size = 20, family = "Times"),
    #         axis.title = element_text(size = 25, family = "Times"),
    #         title = element_text(size = 30, family = "Times"),
    #         legend.title = element_text(size = 20, family = "Times"),
    #         legend.text = element_text(size = 20, family = "Times")
    #     ) +
    #     xlab("Tissue") +
    #     ylab("PhyloP score") +
    #     ggtitle("PhyloP score of selected genes in three tissues") +
    #     # scale_y_continuous(breaks = c(0, 0.5, 1)) +
    #     stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    # ggsave("/media/Data/zhangz/chip/analysis/summary2/anno_sum/common_phylop_boxp.pdf", common_phylop_boxp,
    #     width = 10, height = 8, units = "in"
    # )
    # common_phylop0_boxp <- ggplot(common_gene_ele[common_gene_ele$mean0 != 0, ], aes(x = tissue, y = mean0)) +
    #     geom_violin(aes(fill = tissue)) +
    #     geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    #     theme_bw() +
    #     theme() +
    #     theme(
    #         axis.text = element_text(size = 20, family = "Times"),
    #         axis.title = element_text(size = 25, family = "Times"),
    #         title = element_text(size = 30, family = "Times"),
    #         legend.title = element_text(size = 20, family = "Times"),
    #         legend.text = element_text(size = 20, family = "Times")
    #     ) +
    #     xlab("Tissue") +
    #     ylab("PhyloP score") +
    #     ggtitle("PhyloP score of selected genes in three tissues") +
    #     # scale_y_continuous(breaks = c(0, 0.5, 1)) +
    #     stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    # ggsave("/media/Data/zhangz/chip/analysis/summary2/anno_sum/common_phylop0_boxp.pdf", common_phylop0_boxp,
    #     width = 10, height = 8, units = "in"
    # )
    # levels(common_gene_ele$tissue) <- rev(levels(fct_reorder(common_gene_ele$tissue, common_gene_ele$seq_eln, mean)))
    common_len_boxp <- ggplot(common_gene_ele[which(common_gene_ele$seq_len < 5000), ], aes(x = fct_reorder(tissue, seq_len, mean), y = seq_len)) +
        geom_violin(aes(fill = tissue), trim = TRUE) +
        geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
        theme_bw() +
        theme() +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times")
        ) +
        # ylim(0, 5000) +
        xlab("Tissue") +
        ylab("Sequence length") +
        ggtitle(paste0("Sequence length of selected genes in "), sp) +
        scale_color_manual(values = c("Brain" = "red", "Kidney" = "blue", "Liver" = "green")) +
        # scale_y_continuous(breaks = c(0, 0.5, 1)) +
        stat_compare_means(comparisons = my_comarison, method = "wilcox.test", size = 5, label = "p.signif")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_common_len_boxp.pdf"), common_len_boxp,
        width = 10, height = 8, units = "in"
    )
    write.csv(common_gene_aov, paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_common_gene_aov.csv"), row.names = FALSE)
    write.csv(all_ele, paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_all_ele.csv"), row.names = FALSE)
    ## draw a distribution histogram plot of different element pattern of common genes
    common_gene_ele$ele_pattern <- factor(common_gene_ele$ele_pattern,
        levels = c("B", "K", "L", "BK", "BL", "KL", "BKL")
    )
    # levels(common_gene_ele$ele_pattern) <- c("B", "K", "L", "BK", "BL", "KL", "BKL")
    common_ele_pattern_p <- ggplot(common_gene_ele, aes(x = ele_pattern)) +
        geom_histogram(fill = "grey", color = "black", stat = "count") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times")
        ) +
        xlab("Element pattern") +
        ggtitle("Distribution of element pattern annotated \nwith common genes among tissues")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", ele, "_common_ele_pattern_p.pdf"), common_ele_pattern_p,
        width = 10, height = 8, units = "in"
    )
}
# }
