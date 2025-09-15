## for a species, summary all the genes annotated,
## and for each gene, how many enhancers or promoters
## are annotated to it

# library(optparse)

# option_list <- list(
#     make_option(c("-s", "--species"),
#         type = "character", default = "Mus_musculus",
#         help = "The species of the data"
#     )
# )
# opt_parser <- parse_args(OptionParser(option_list = option_list, usage = "%prog [options]"))
# sp <- opt_parser$species
# print(sp)

# load packages
library(stringr, quietly = TRUE)
library(MuMIn, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(plyr, quietly = TRUE)
library(doMC, quietly = TRUE)
library(ggpubr, quietly = TRUE)
library(gghalves, quietly = TRUE)
# read in data
# sp_all <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_all.csv"),
#     header = TRUE
# )
doMC::registerDoMC(cores = 7)
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
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/info_concat.RData")
res <- adply(species, 1, function(sp) {
    ele_info <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno_new.csv"), header = TRUE)
    genes <- unique(as.vector(unlist(sapply(ele_info$overlap_gene, function(x) strsplit(x, ";")[[1]]))))
    genes <- genes[!is.na(genes)]
    sp_res <- adply(genes, 1, function(g) {
        tmp <- ele_info[grep(regex(paste0("(^|;)", g, "(;|$)")), ele_info$overlap_gene), ]
        brain_enhancer <- sum(tmp$element == "enhancer" & tmp$tissue == "Brain")
        brain_promoter <- sum(tmp$element == "promoter" & tmp$tissue == "Brain")
        kidney_enhancer <- sum(tmp$element == "enhancer" & tmp$tissue == "Kidney")
        kidney_promoter <- sum(tmp$element == "promoter" & tmp$tissue == "Kidney")
        liver_enhancer <- sum(tmp$element == "enhancer" & tmp$tissue == "Liver")
        liver_promoter <- sum(tmp$element == "promoter" & tmp$tissue == "Liver")
        all_enhancer <- tmp[tmp$element == "enhancer", "id"] |>
            unique() |>
            length()
        all_promoter <- tmp[tmp$element == "promoter", "id"] |>
            unique() |>
            length()
        data.frame(
            gene = g, Brain_enhancer = brain_enhancer, Brain_promoter = brain_promoter,
            Liver_enhancer = liver_enhancer, Liver_promoter = liver_promoter,
            Kidney_enhancer = kidney_enhancer, Kidney_promoter = kidney_promoter,
            all_enhancer = all_enhancer, all_promoter = all_promoter
        )
    }, .parallel = TRUE)
    sp_res$species <- sp
    write.csv(sp_res, paste0("/media/Data/zhangz/chip/analysis/summary2/each_sp/", sp, "_cre_each_gene_filt.csv"), row.names = FALSE)
    return(sp_res)
}, .parallel = TRUE)
write.csv(res, "/media/Data/zhangz/chip/analysis/summary2/each_sp/cre_per_gene.csv", row.names = FALSE)
## plot distribution density of all enhancers and promoters
# res <- read.csv("/media/Data/zhangz/chip/analysis/summary2/each_sp/cre_per_gene.csv")
res_m <- reshape2::melt(res, id.vars = c("gene", "species"), measure.vars = c("all_enhancer", "all_promoter"), value.name = "count")
cre_p <- ggplot(res_m, aes(x = count, fill = variable)) +
    geom_density(alpha = 0.5, linewidth = 0.5) +
    # facet_wrap(~variable, scales = "free") +
    theme_classic() +
    # theme(legend.position = "none") +
    labs(x = "Number of CREs", y = "Density", fill = "CRE") +
    scale_fill_manual(values = c("#c01d2e", "#03324e"), labels = c("Enhancer", "Promoter")) +
    scale_x_continuous(n.breaks = 3, limits = c(0, 10)) +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.8)) +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    theme(text = element_text(family = "Arial"))
# cre_p_zoom <- cre_p +
#     geom_rect(aes(xmin = 0, xmax = 10, ymin = 0, ymax = 0.1), fill = "grey", alpha = 0.2) +
#     coord_cartesian(xlim = c(0, 10))
# gridExtra::grid.arrange(cre_p, cre_p_zoom, ncol = 2)
# ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/cre_per_gene_density_zoom.png", cre_p_zoom, width = 6, height = 6, units = "cm")
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/cre_per_gene_density.png", cre_p, width = 6, height = 6, units = "cm")
cre_hist_p <- ggplot(res_m, aes(x = count, fill = variable)) +
    # geom_density(alpha = 0.5, linewidth = 0.5) +
    geom_histogram(binwidth = 1, position = "dodge", alpha = 0.5) +
    # facet_wrap(~variable, scales = "free") +
    theme_classic() +
    # theme(legend.position = "none") +
    labs(x = "Number of CREs", y = "Density", fill = "CRE") +
    scale_fill_manual(values = c("#c01d2e", "#03324e"), labels = c("Enhancer", "Promoter")) +
    scale_x_continuous(n.breaks = 3, limits = c(0, 20)) +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.8)) +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/cre_per_gene_hist.png", cre_hist_p, width = 6, height = 6, units = "cm")

res_mid <- tapply(res_m$count, list(res_m$variable, res_m$species), median) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var1") %>%
    gather(key = "Var2", value = "value", -Var1)
# box plot
res_mid_p <- ggplot(melt(res_mid), aes(x = Var1, y = value, fill = Var1)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    labs(x = "", y = "Median CREs number per gene", fill = "CRE") +
    scale_fill_manual(values = c("#ef4343", "#457b9d"), labels = c("Enhancer", "Promoter")) +
    scale_x_discrete(labels = c("Enhancer", "Promoter")) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3, label.x = 1.45) +
    theme(legend.position = "None") +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/cre_per_gene_box.pdf", res_mid_p, width = 6, height = 6, units = "cm")

## i need the median distance and annotated genes number per cre of each species
dis_gene <- adply(species, 1, function(sp) {
    ele_info <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno_new.csv"), header = TRUE)
    anno_set <- unique(ele_info$note_abc)
    abc_set <- anno_set[grep("abc", anno_set)]
    other_set <- anno_set[!grepl("abc", anno_set)]
    dis_median <- paste(ele_info$overlap_distance, collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    abc_dis_median <- paste(ele_info[ele_info$note_abc %in% abc_set, "overlap_distance"], collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    other_dis_median <- paste(ele_info[ele_info$note_abc %in% other_set, "overlap_distance"], collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    ele_info$gene_num <- sapply(ele_info$overlap_gene, function(x) length(strsplit(x, ";")[[1]]))
    gene_num_median <- median(ele_info$gene_num, na.rm = TRUE)
    abc_gene_num_median <- median(ele_info[ele_info$note_abc %in% abc_set, "gene_num"], na.rm = TRUE)
    other_gene_num_median <- median(ele_info[ele_info$note_abc %in% other_set, "gene_num"], na.rm = TRUE)
    res <- data.frame(
        species = sp, dis_median = dis_median, abc_dis_median = abc_dis_median,
        other_dis_median = other_dis_median, gene_num_median = gene_num_median,
        abc_gene_num_median = abc_gene_num_median, other_gene_num_median = other_gene_num_median
    )
    res1_m <- melt(res, id.vars = "species", measure.vars = c("dis_median", "abc_dis_median", "other_dis_median", "gene_num_median", "abc_gene_num_median", "other_gene_num_median"), value.name = "median")
    res1_m$ele <- "all"
    dis_median_enhancer <- paste(ele_info[ele_info$element == "enhancer", "overlap_distance"], collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    abc_dis_median_enhancer <- paste(ele_info[ele_info$note_abc %in% abc_set & ele_info$element == "enhancer", "overlap_distance"], collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    other_dis_median_enhancer <- paste(ele_info[ele_info$note_abc %in% other_set & ele_info$element == "enhancer", "overlap_distance"], collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    gene_num_median_enhancer <- median(ele_info[ele_info$element == "enhancer", "gene_num"], na.rm = TRUE)
    abc_gene_num_median_enhancer <- median(ele_info[ele_info$note_abc %in% abc_set & ele_info$element == "enhancer", "gene_num"], na.rm = TRUE)
    other_gene_num_median_enhancer <- median(ele_info[ele_info$note_abc %in% other_set & ele_info$element == "enhancer", "gene_num"], na.rm = TRUE)
    res2 <- data.frame(
        species = sp, dis_median = dis_median_enhancer, abc_dis_median = abc_dis_median_enhancer,
        other_dis_median = other_dis_median_enhancer, gene_num_median = gene_num_median_enhancer,
        abc_gene_num_median = abc_gene_num_median_enhancer, other_gene_num_median = other_gene_num_median_enhancer
    )
    res2_m <- melt(res2, id.vars = "species", measure.vars = c("dis_median", "abc_dis_median", "other_dis_median", "gene_num_median", "abc_gene_num_median", "other_gene_num_median"), value.name = "median")
    res2_m$ele <- "enhancer"
    dis_median_promoter <- paste(ele_info[ele_info$element == "promoter", "overlap_distance"], collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    abc_dis_median_promoter <- paste(ele_info[ele_info$note_abc %in% abc_set & ele_info$element == "promoter", "overlap_distance"], collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    other_dis_median_promoter <- paste(ele_info[ele_info$note_abc %in% other_set & ele_info$element == "promoter", "overlap_distance"], collapse = ";") |>
        strsplit(";") |>
        unlist() |>
        as.numeric() |>
        abs() |>
        median(., na.rm = TRUE)
    gene_num_median_promoter <- median(ele_info[ele_info$element == "promoter", "gene_num"], na.rm = TRUE)
    abc_gene_num_median_promoter <- median(ele_info[ele_info$note_abc %in% abc_set & ele_info$element == "promoter", "gene_num"], na.rm = TRUE)
    other_gene_num_median_promoter <- median(ele_info[ele_info$note_abc %in% other_set & ele_info$element == "promoter", "gene_num"], na.rm = TRUE)
    res3 <- data.frame(
        species = sp, dis_median = dis_median_promoter, abc_dis_median = abc_dis_median_promoter,
        other_dis_median = other_dis_median_promoter, gene_num_median = gene_num_median_promoter,
        abc_gene_num_median = abc_gene_num_median_promoter, other_gene_num_median = other_gene_num_median_promoter
    )
    res3_m <- melt(res3, id.vars = "species", measure.vars = c("dis_median", "abc_dis_median", "other_dis_median", "gene_num_median", "abc_gene_num_median", "other_gene_num_median"), value.name = "median")
    res3_m$ele <- "promoter"
    res4 <- rbind(res1_m, res2_m, res3_m)
    return(res4)
}, .parallel = TRUE)

write.csv(dis_gene, "/media/Data/zhangz/chip/analysis/summary2/each_sp/dis_gene_num.csv", row.names = FALSE)
# draw boxplot of median distance

dis_gene <- dis_gene[!is.na(dis_gene$median), ]
dis_var <- c("abc_dis_median", "other_dis_median") # "dis_median",
gene_num_var <- c("abc_gene_num_median", "other_gene_num_median") # "gene_num_median",
genenum_p <- ggplot(dis_gene[dis_gene$variable %in% gene_num_var, ], aes(x = variable, y = median, fill = variable)) +
    geom_boxplot(outlier.size = 0.5, fill = "transparent", width = 0.1) +
    geom_half_violin(alpha = 0.5, side = "top") +
    geom_half_point(aes(color = variable), side = "1", alpha = 0.4, size = 0.5) +
    facet_grid(ele ~ .) +
    coord_flip() +
    theme_classic() +
    labs(x = "", y = "Median number of associated genes", fill = "Annotationn", color = "Annotation")

# scale_fill_manual(values = c("#ef4343", "#457b9d", "#f4a261"), labels = c("All", "Enhancer", "Promoter")) + ()
# try histogram plot
genenum_hist_p <- ggplot(dis_gene[which(dis_gene$variable %in% gene_num_var & dis_gene$ele != "all"), ], aes(x = median, fill = variable)) +
    geom_histogram(binwidth = 1, position = "dodge", color = "black", linewidth = 0.2) +
    geom_boxploth(width = 3, aes(y = -4), size = 0.1, outlier.shape = NA) +
    # geom_half_point(aes(x = median, y = -6, color = variable), side = "1", alpha = 0.4, size = 0.5) +
    theme_minimal() +
    facet_grid(ele ~ .) +
    labs(x = "Median number of associated genes", y = "Count", fill = "Annotation", color = "") +
    scale_fill_manual(values = c("#ef4343", "#457b9d"), labels = c("With ABC", "Without ABC")) +
    scale_x_continuous(n.breaks = 3, breaks = c(1, 5, 10)) +
    theme(legend.position = "top") + # , legend.position.inside = c(0.8, 0.9)
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    theme(panel.grid = element_blank()) +
    theme(legend.key.size = unit(0.3, "cm"))
# theme(text = element_text(family = "arial"))
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/gene_num_hist.pdf", genenum_hist_p, width = 6, height = 6, units = "cm")

# test, fill by ele, facet by variable
facet_labels1 <- c("abc_gene_num_median" = "With ABC", "other_gene_num_median" = "Without ABC")
genenum_hist_p2 <- ggplot(dis_gene[which(dis_gene$variable %in% gene_num_var & dis_gene$ele != "all"), ], aes(x = median, fill = ele)) +
    geom_histogram(binwidth = 1, position = "dodge", color = "black", linewidth = 0.2) +
    geom_boxploth(width = 3, aes(y = -4), size = 0.1, outlier.shape = NA) +
    # geom_half_point(aes(x = median, y = -6, color = variable), side = "1", alpha = 0.4, size = 0.5) +
    theme_minimal() +
    facet_grid(variable ~ ., labeller = labeller(category = facet_labels1)) +
    labs(x = "Median number of associated genes", y = "Count", fill = "Annotation", color = "") +
    scale_fill_manual(values = c("#ef4343", "#457b9d"), labels = c("Enhancer", "Promoter")) +
    scale_x_continuous(n.breaks = 3, breaks = c(1, 5, 10)) +
    theme(legend.position = "top") + # , legend.position.inside = c(0.8, 0.9)
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    theme(panel.grid = element_blank()) +
    theme(legend.key.size = unit(0.3, "cm"))
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/gene_num_hist2.pdf", genenum_hist_p2, width = 6, height = 6, units = "cm")

summary(dis_gene[which(dis_gene$variable %in% dis_var & dis_gene$ele != "all" & dis_gene$variable == "abc_dis_median"), ])
dis_hist_p <- ggplot(dis_gene[which(dis_gene$variable %in% dis_var & dis_gene$ele != "all"), ], aes(x = median, fill = variable)) +
    geom_histogram(binwidth = 1, position = "dodge", color = "black", linewidth = 0.2) +
    geom_boxploth(width = 3, aes(x = median, y = -4), size = 0.1, outlier.shape = NA) +
    # geom_half_point(aes(x = median, y = -6, color = variable), side = "1", alpha = 0.4, size = 0.5) +
    theme_minimal() +
    facet_grid(ele ~ .) +
    labs(x = "Median distance to candidate target gene", y = "Count", fill = "Annotation", color = "") +
    scale_fill_manual(values = c("#ef4343", "#457b9d"), labels = c("With ABC", "Without ABC")) +
    scale_x_continuous(breaks = c(1, 10, 100, 10000, 1000000), trans = "log10", labels = scales::comma) +
    theme(legend.position = "top") + # , legend.position.inside = c(0.8, 0.9)
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7)) +
    theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) +
    theme(panel.grid = element_blank()) +
    theme(legend.key.size = unit(0.3, "cm"))
# theme(text = element_text(family = "arial"))
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/dis_hist.pdf", dis_hist_p, width = 6, height = 6, units = "cm")
# fill by ele, facet by variable
facet_labels2 <- c("abc_dis_median" = "With ABC", "other_dis_median" = "Without ABC")
dis_hist_p2 <- ggplot(dis_gene[which(dis_gene$variable %in% dis_var & dis_gene$ele != "all"), ], aes(x = median, fill = ele)) +
    geom_histogram(binwidth = 1, position = "dodge", color = "black", linewidth = 0.2) +
    geom_boxploth(width = 3, aes(x = median, y = -4), size = 0.1, outlier.shape = NA) +
    # geom_half_point(aes(x = median, y = -6, color = variable), side = "1", alpha = 0.4, size = 0.5) +
    theme_minimal() +
    facet_grid(variable ~ ., labeller = labeller(category = facet_labels2)) +
    labs(x = "Median distance to candidate target gene", y = "Count", fill = "Annotation", color = "") +
    scale_fill_manual(values = c("#ef4343", "#457b9d"), labels = c("Enhancer", "Promoter")) +
    scale_x_continuous(breaks = c(1, 10, 100, 10000, 1000000), trans = "log10", labels = scales::comma) +
    theme(legend.position = "top") + # , legend.position.inside = c(0.8, 0.9)
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7)) +
    theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) +
    theme(panel.grid = element_blank()) +
    theme(legend.key.size = unit(0.3, "cm"))
ggsave("/media/Data/zhangz/chip/analysis/summary2/each_sp/dis_hist2.pdf", dis_hist_p2, width = 6, height = 6, units = "cm")


# cre_per_gene <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/each_sp/", sp, "_cre_each_gene.csv"))
# res2 <- data.frame()
# for (i in 2:ncol(res)) {
#     res2 <- rbind(res2, data.frame(
#         min = min(res[, i]),
#         f_th = quantile(res[, i], 0.25),
#         median = median(res[, i]),
#         mean = mean(res[, i]),
#         t_th = quantile(res[, i], 0.75),
#         max = max(res[, i]), row.names = colnames(res)[i]
#     ))
# }
# res2$name <- rownames(res2)
# res3 <- t(melt(res2) |> unite(names, name, variable, sep = "_"))
# colnames(res3) <- res3[1, ]
# res3 <- as.data.frame(res3)
# res3["value", ] <- as.numeric(res3[2, ])
# res3 <- res3[-1, ]
# res3$gene_num <- length(genes)
# res3$species <- sp
# write.csv(res3, paste0("/media/Data/zhangz/chip/analysis/summary2/each_sp/", sp, "_cre_each_gene_summary.csv"), row.names = FALSE)
