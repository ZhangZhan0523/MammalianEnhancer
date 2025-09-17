## summarize the correlation result of each species,
## mainly entropy associated, calculated with entropy of
## egi from abc overlapped with 1Mts

library(ggplot2)
library(reshape2)
library(magrittr)

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

ent_rho_spm <- data.frame()
ent_p_spm <- data.frame()
for (sp in species_ts) {
    # for (tis in c("Brain", "Kidney", "Liver")) {
    tmp_spm <- try(read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent_cor.csv"),
        header = TRUE
    ))
    if (class(tmp_spm) == "try-error") {
        next
    }
    tmp_spm$term <- paste(tmp_spm$Var1, tmp_spm$Var2, sep = "VS")
    tmp <- melt(tmp_spm[c("rho", "p", "term")], id.vars = "term", variable.names = "stat_name", value.name = "stat")
    tmp <- tmp[order(tmp$term), ]
    # tmp$name <- paste(tmp$term, tmp$variable, sep = "_")
    rho_tmp <- as.data.frame(t(tmp[tmp$variable == "rho", c("term", "stat")]))
    p_tmp <- tmp[tmp$variable == "p", c("term", "stat")] %>%
        t() %>%
        as.data.frame()
    colnames(rho_tmp) <- rho_tmp[1, ]
    colnames(p_tmp) <- p_tmp[1, ]
    ent_rho_spm <- rbind(ent_rho_spm, cbind(rho_tmp, sp = sp)[-1, ])
    ent_p_spm <- rbind(ent_p_spm, cbind(p_tmp, sp = sp)[-1, ])
    # }
}
rownames(ent_rho_spm) <- 1:nrow(ent_rho_spm)
rownames(ent_p_spm) <- 1:nrow(ent_p_spm)
ent_rho_spm[, 1:78] <- lapply(ent_rho_spm[, 1:78], as.numeric)
ent_p_spm[, 1:78] <- lapply(ent_p_spm[, 1:78], as.numeric)
ent_sum <- data.frame(term = rep(colnames(ent_rho_spm)[1:78], nrow(ent_rho_spm)), rho = 0, p = 0)
summary(ent_sum[ent_sum$term == "unique_ele_numVSent1", ])
terms <- c()
for (term in colnames(ent_rho_spm)[1:78]) {
    if (median(ent_rho_spm[, term]) > 0.5 & max(ent_p_spm[, term]) < 0.05) {
        print(term)
        terms <- c(terms, term)
    }
    ent_sum[ent_sum$term == term, "rho"] <- ent_rho_spm[, term]
    ent_sum[ent_sum$term == term, "p"] <- ent_p_spm[, term]
}
terms <- c(
    "ent1VSfpkm", "ent2VSfpkm", "ent1VSfpkm_sd", "ent2VSfpkm_sd",
    "ent1VSfpkm_cv", "ent2VSfpkm_cv",
    "ent1VSadj_fpkm", "ent2VSadj_fpkm",
    "ent1VSadj_fpkm_sd", "ent2VSadj_fpkm_sd",
    "ent1VSadj_fpkm_cv", "ent2VSadj_fpkm_cv",
    "ent1VSqn_fpkm_sps", "ent2VSqn_fpkm_sps",
    "ent1VSqn_fpkm_sd", "ent2VSqn_fpkm_sd",
    "ent1VSqn_fpkm_cv", "ent2VSqn_fpkm_cv",
    "ent1VStau", "ent2VStau", "ent2VSent1"
)
ent_rho_p_df <- ent_sum[ent_sum$term %in% terms, ]
ent_rho_p_df$term <- factor(ent_rho_p_df$term, levels = c(
    "ent2VSent1",
    "ent1VSfpkm", "ent2VSfpkm", "ent1VSfpkm_sd", "ent2VSfpkm_sd",
    "ent1VSfpkm_cv", "ent2VSfpkm_cv",
    "ent1VSadj_fpkm", "ent2VSadj_fpkm",
    "ent1VSadj_fpkm_sd", "ent2VSadj_fpkm_sd",
    "ent1VSadj_fpkm_cv", "ent2VSadj_fpkm_cv",
    "ent1VSqn_fpkm_sps", "ent2VSqn_fpkm_sps",
    "ent1VSqn_fpkm_sd", "ent2VSqn_fpkm_sd",
    "ent1VSqn_fpkm_cv", "ent2VSqn_fpkm_cv",
    "ent1VStau", "ent2VStau"
))
ent_rho_p <- ggplot(ent_rho_p_df, aes(x = term, y = rho, fill = term)) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = c(
        "ent2VSent1" = "#c01d2e",
        "ent1VSfpkm" = "#ed4343", "ent2VSfpkm" = "#ed4343",
        "ent1VSfpkm_sd" = "#ed4343", "ent2VSfpkm_sd" = "#ed4343",
        "ent1VSfpkm_cv" = "#ed4343", "ent2VSfpkm_cv" = "#ed4343",
        "ent1VSadj_fpkm" = "#ed4343", "ent2VSadj_fpkm" = "#ed4343",
        "ent1VSadj_fpkm_sd" = "#ed4343", "ent2VSadj_fpkm_sd" = "#ed4343",
        "ent1VSadj_fpkm_cv" = "#ed4343", "ent2VSadj_fpkm_cv" = "#ed4343",
        "ent1VSqn_fpkm_sps" = "#033250", "ent2VSqn_fpkm_sps" = "#033250",
        "ent1VSqn_fpkm_sd" = "#033250", "ent2VSqn_fpkm_sd" = "#033250",
        "ent1VSqn_fpkm_se" = "#033250", "ent2VSqn_fpkm_se" = "#033250",
        "ent1VStau" = "#73b8d5", "ent2VStau" = "#73b8d5"
    )) +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 45)
    ) +
    geom_abline(intercept = c(-0.5, 0, 0.5), slope = 0, linetype = "dashed", color = "black") +
    ylim(-1, 1) +
    ggtitle(paste0("spearman test of CRE entropies and \nassociated gene express (raw FPKM) \nacross genes in a biosample (abc)")) +
    xlab("") +
    ylab("rho")
ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_abc.pdf"), plot = ent_rho_p, width = 12, height = 12, dpi = 600)
ent_rho_p_df2 <- plyr::adply(ent_rho_p_df, 1, function(x) {
    if (TRUE %in% grepl("ent1VS", x$term)) {
        x$ent_type <- "ent1"
    } else {
        x$ent_type <- "ent2"
    }
    if (TRUE %in% grepl("adj", x$term)) {
        fpkm_type <- "rep&adj"
    } else if (TRUE %in% grepl("qn", x$term)) {
        fpkm_type <- "qn"
    } else {
        fpkm_type <- "raw"
    }
    x$fpkm_type <- fpkm_type
    if (TRUE %in% grepl("sd", x$term)) {
        x$fpkm_index <- "sd"
    } else if (TRUE %in% grepl("cv", x$term)) {
        x$fpkm_index <- "cv"
    } else {
        x$fpkm_index <- "mean"
    }
    return(x)
})
ent_rho_p_df2 <- ent_rho_p_df2[!(ent_rho_p_df2$term %in% c("ent2VSent1", "ent1VStau", "ent2VStau")), ]
ent1_rho_p_point <- ggplot(ent_rho_p_df2[ent_rho_p_df2$ent_type == "ent1", ], aes(y = rho, x = -log(p))) +
    geom_point(aes(color = fpkm_type, shape = fpkm_index)) +
    theme_bw() +
    scale_color_manual(values = c("raw" = "#ed4343", "qn" = "#033250", "rep&adj" = "#73b8d5")) +
    scale_shape_manual(values = c("mean" = 16, "sd" = 17, "cv" = 18)) +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank()
    ) +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5)) +
    facet_grid(fpkm_index ~ fpkm_type) +
    ggtitle("spearman test of CRE entropies and \nassociated gene express (raw FPKM) \nacross genes in a biosample (abc)") +
    xlab("rho") +
    ylab("-log(p)")
ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_abc_point.pdf"), plot = ent1_rho_p_point, width = 12, height = 12, dpi = 600)

ent1_rho_box_p <- ggplot(ent_rho_p_df2[ent_rho_p_df2$ent_type == "ent1", ], aes(x = interaction(fpkm_type, fpkm_index), y = rho, fill = interaction(fpkm_type, fpkm_index))) +
    geom_boxplot(linewidth = 0.1, outlier.size = 0.1) +
    theme_bw() +
    scale_fill_manual(values = c(
        "raw.mean" = "#ed4343", "raw.sd" = "#ed4343", "raw.cv" = "#ed4343",
        "qn.mean" = "#033250", "qn.sd" = "#033250", "qn.cv" = "#033250",
        "rep&adj.mean" = "#73b8d5", "rep&adj.sd" = "#73b8d5", "rep&adj.cv" = "#73b8d5"
    )) +
    theme(
        axis.text = element_text(size = 8, family = "Times"),
        axis.title = element_text(size = 8, family = "Times"),
        title = element_text(size = 8, family = "Times"),
        legend.title = element_text(size = 8, family = "Times"),
        legend.text = element_text(size = 8, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 45)
    ) +
    geom_abline(intercept = c(-0.5, 0, 0.5), slope = 0, linetype = "dashed", color = "black") +
    ylim(-1, 1) +
    scale_x_discrete(labels = c(
        "raw.mean" = "raw mean", "raw.sd" = "raw sd", "raw.cv" = "raw cv",
        "qn.mean" = "qn mean", "qn.sd" = "qn sd", "qn.cv" = "qn cv",
        "rep&adj.mean" = "rep&adj mean", "rep&adj.sd" = "rep&adj sd", "rep&adj.cv" = "rep&adj cv"
    ), limits = c("raw.mean", "raw.sd", "raw.cv", "qn.mean", "qn.sd", "qn.cv", "rep&adj.mean", "rep&adj.sd", "rep&adj.cv")) +
    ggtitle("spearman rho of gene entropies and express (raw FPKM) \nacross genes in a biosample (abc)") +
    xlab("") +
    ylab("spearman rho")
ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_abc_box.pdf"), plot = ent1_rho_box_p, width = 10, height = 8, dpi = 600, units = "cm")


summary(ent_rho_p_df2[which(ent_rho_p_df2$fpkm_type == "raw" & ent_rho_p_df2$fpkm_index == "mean" & ent_rho_p_df2$ent_type == "ent1"), ])
summary(ent_rho_p_df2[which(ent_rho_p_df2$fpkm_type == "raw" & ent_rho_p_df2$fpkm_index == "cv" & ent_rho_p_df2$ent_type == "ent1"), ])
summary(ent_rho_p_df2[which(ent_rho_p_df2$fpkm_type == "rep&adj" & ent_rho_p_df2$fpkm_index == "mean" & ent_rho_p_df2$ent_type == "ent1"), ])
summary(ent_rho_p_df2[which(ent_rho_p_df2$fpkm_type == "rep&adj" & ent_rho_p_df2$fpkm_index == "cv" & ent_rho_p_df2$ent_type == "ent1"), ])
## try to re-draw and save a new plot for publication
ent1_rho_box_p <- ggplot(ent_rho_p_df2[which(ent_rho_p_df2$ent_type == "ent1" & ent_rho_p_df2$fpkm_type %in% c("rep&adj", "qn")), ], aes(x = interaction(fpkm_type, fpkm_index), y = rho, fill = fpkm_type)) +
    geom_boxplot(linewidth = 0.2, outlier.size = 0.1) +
    theme_light() +
    scale_fill_manual(values = c(
        # "raw.mean" = "#ed4343", "raw.sd" = "#ed4343", "raw.cv" = "#ed4343",
        # "qn.mean" = "#ed4343", "qn.sd" = "#ed4343", "qn.cv" = "#ed4343",
        # "rep&adj.mean" = "#73b8d5", "rep&adj.sd" = "#73b8d5", "rep&adj.cv" = "#73b8d5"
        "qn" = "#ed4343", "rep&adj" = "#73b8d5"
    ), labels = c("qn" = "1 rep", "rep&adj" = "multi-rep")) +
    theme(
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        title = element_text(size = 8, hjust = 0.5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "top",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust = 0.5)
    ) +
    geom_abline(intercept = c(-0.5, 0, 0.5), slope = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    # ylim(-0.6, 0.6) +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5)) +
    scale_x_discrete(labels = c(
        # "raw.mean" = "raw mean", "raw.sd" = "raw sd", "raw.cv" = "raw cv",
        "qn.mean" = "mean", "qn.sd" = "sd", "qn.cv" = "cv",
        "rep&adj.mean" = "mean", "rep&adj.sd" = "sd", "rep&adj.cv" = "cv"
    ), limits = c("qn.mean", "qn.sd", "qn.cv", "rep&adj.mean", "rep&adj.sd", "rep&adj.cv")) +
    # "raw.mean", "raw.sd", "raw.cv",
    ggtitle("Correlation of entropies and\nexpression across genes") +
    xlab("") +
    ylab("Spearman rho") +
    labs(fill = "Expression") +
    theme(legend.spacing = unit(0.05, "cm")) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.box.margin = margin(0, 0, 0, 0)) +
    theme(legend.margin = margin(0, 0, 0, 0))
ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_abc_box.pdf"), plot = ent1_rho_box_p, width = 6, height = 6, dpi = 600, units = "cm")


ent2_rho_p_point <- ggplot(ent_rho_p_df2[ent_rho_p_df2$ent_type == "ent2", ], aes(y = rho, x = -log(p))) +
    geom_point(aes(color = fpkm_type, shape = fpkm_index)) +
    theme_bw() +
    scale_color_manual(values = c("raw" = "#ed4343", "qn" = "#033250", "rep&adj" = "#73b8d5")) +
    scale_shape_manual(values = c("mean" = 16, "sd" = 17, "cv" = 18)) +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank()
    ) +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5)) +
    facet_grid(fpkm_index ~ fpkm_type) +
    ggtitle("spearman test of CRE entropies 2 and \nassociated gene express (raw FPKM) \nacross genes in a biosample (abc)") +
    xlab("rho") +
    ylab("-log(p)")
ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_abc_point2.pdf"), plot = ent2_rho_p_point, width = 12, height = 12, dpi = 600)

ent2_rho_box_p <- ggplot(ent_rho_p_df2[ent_rho_p_df2$ent_type == "ent2", ], aes(x = interaction(fpkm_type, fpkm_index), y = rho, fill = interaction(fpkm_type, fpkm_index))) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = c(
        "raw.mean" = "#ed4343", "raw.sd" = "#ed4343", "raw.cv" = "#ed4343",
        "qn.mean" = "#033250", "qn.sd" = "#033250", "qn.cv" = "#033250",
        "rep&adj.mean" = "#73b8d5", "rep&adj.sd" = "#73b8d5", "rep&adj.cv" = "#73b8d5"
    )) +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 45)
    ) +
    geom_abline(intercept = c(-0.5, 0, 0.5), slope = 0, linetype = "dashed", color = "black") +
    ylim(-1, 1) +
    scale_x_discrete(labels = c(
        "raw.mean" = "raw mean", "raw.sd" = "raw sd", "raw.cv" = "raw cv",
        "qn.mean" = "qn mean", "qn.sd" = "qn sd", "qn.cv" = "qn cv",
        "rep&adj.mean" = "rep&adj mean", "rep&adj.sd" = "rep&adj sd", "rep&adj.cv" = "rep&adj cv"
    ), limits = c("raw.mean", "raw.sd", "raw.cv", "qn.mean", "qn.sd", "qn.cv", "rep&adj.mean", "rep&adj.sd", "rep&adj.cv")) +
    ggtitle("spearman test of CRE entropies 2 and \nassociated gene express (raw FPKM) \nacross genes in a biosample (abc)") +
    xlab("") +
    ylab("spearman rho")
ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_abc_box2.pdf"), plot = ent2_rho_box_p, width = 12, height = 12, dpi = 600)

ent_rho_spm_no0 <- data.frame()
ent_p_spm_no0 <- data.frame()
for (sp in species_ts) {
    tmp_spm <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/", sp, "_", dis, "_ent_cor_no0_adj.csv"),
        header = TRUE
    )
    tmp_spm$term <- paste(tmp_spm$Var1, tmp_spm$Var2, sep = "VS")
    tmp <- melt(tmp_spm[c("rho", "p", "term")], id.vars = "term", variable.names = "stat_name", value.name = "stat")
    tmp <- tmp[order(tmp$term), ]
    # tmp$name <- paste(tmp$term, tmp$variable, sep = "_")
    rho_tmp <- as.data.frame(t(tmp[tmp$variable == "rho", c("term", "stat")]))
    p_tmp <- tmp[tmp$variable == "p", c("term", "stat")] %>%
        t() %>%
        as.data.frame()
    colnames(rho_tmp) <- rho_tmp[1, ]
    colnames(p_tmp) <- p_tmp[1, ]
    ent_rho_spm_no0 <- rbind(ent_rho_spm_no0, cbind(rho_tmp, sp = sp)[-1, ])
    ent_p_spm_no0 <- rbind(ent_p_spm_no0, cbind(p_tmp, sp = sp)[-1, ])
}
rownames(ent_rho_spm_no0) <- 1:nrow(ent_rho_spm_no0)
rownames(ent_p_spm_no0) <- 1:nrow(ent_p_spm_no0)
ent_rho_spm_no0[, 1:78] <- lapply(ent_rho_spm_no0[, 1:78], as.numeric)
ent_p_spm_no0[, 1:78] <- lapply(ent_p_spm_no0[, 1:78], as.numeric)
ent_sum_no0 <- data.frame(term = rep(colnames(ent_rho_spm_no0)[1:78], 21), rho = 0, p = 0)
ent_no0_terms <- c()
for (term in colnames(ent_rho_spm_no0)[1:78]) {
    if (median(ent_rho_spm_no0[, term]) > 0.5 & max(ent_p_spm_no0[, term]) < 0.05) {
        print(term)
        ent_no0_terms <- c(ent_no0_terms, term)
    }
    ent_sum_no0[ent_sum_no0$term == term, "rho"] <- ent_rho_spm_no0[, term]
    ent_sum_no0[ent_sum_no0$term == term, "p"] <- ent_p_spm_no0[, term]
}
terms <- c(
    "no0_ent1VSno0_fpkm", "no0_ent2VSno0_fpkm", "no0_ent1VSno0_fpkm_sd",
    "no0_ent2VSno0_fpkm_sd", "no0_ent1VSno0_fpkm_cv", "no0_ent2VSno0_fpkm_cv",
    "no0_ent1VSno0_qn_fpkm_sps", "no0_ent2VSno0_qn_fpkm_sps", "no0_ent1VSno0_qn_fpkm_sd",
    "no0_ent2VSno0_qn_fpkm_sd", "no0_ent1VSno0_qn_fpkm_cv", "no0_ent2VSno0_qn_fpkm_cv",
    "no0_adj_fpkmVSno0_ent1", "no0_adj_fpkmVSno0_ent2", "no0_adj_fpkm_sdVSno0_ent1",
    "no0_adj_fpkm_sdVSno0_ent2", "no0_adj_fpkm_cvVSno0_ent1", "no0_adj_fpkm_cvVSno0_ent2",
    "no0_ent1VSno0_tau", "no0_ent2VSno0_tau", "no0_ent2VSno0_ent1"
)
ent_rho_spm_no0_df <- ent_sum_no0[ent_sum_no0$term %in% terms, ]
ent_rho_spm_no0_df$term <- factor(ent_rho_spm_no0_df$term, levels = c(
    "no0_ent1VSno0_fpkm", "no0_ent2VSno0_fpkm",
    "no0_ent1VSno0_fpkm_sd", "no0_ent2VSno0_fpkm_sd",
    "no0_ent1VSno0_fpkm_cv", "no0_ent2VSno0_fpkm_cv",
    "no0_ent1VSno0_qn_fpkm_sps", "no0_ent2VSno0_qn_fpkm_sps",
    "no0_ent1VSno0_qn_fpkm_sd", "no0_ent2VSno0_qn_fpkm_sd",
    "no0_ent1VSno0_qn_fpkm_cv", "no0_ent2VSno0_qn_fpkm_cv",
    "no0_adj_fpkmVSno0_ent1", "no0_adj_fpkmVSno0_ent2",
    "no0_adj_fpkm_sdVSno0_ent1", "no0_adj_fpkm_sdVSno0_ent2",
    "no0_adj_fpkm_cvVSno0_ent1", "no0_adj_fpkm_cvVSno0_ent2",
    "no0_ent1VSno0_tau", "no0_ent2VSno0_tau", "no0_ent2VSno0_ent1"
))
ent_rho_no0_p <- ggplot(ent_rho_spm_no0_df, aes(x = term, y = rho, fill = term)) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = c(
        "no0_ent1VSno0_fpkm" = "#ed4343", "no0_ent2VSno0_fpkm" = "#ed4343",
        "no0_ent1VSno0_fpkm_sd" = "#ed4343", "no0_ent2VSno0_fpkm_sd" = "#ed4343",
        "no0_ent1VSno0_fpkm_cv" = "#ed4343", "no0_ent2VSno0_fpkm_cv" = "#ed4343",
        "no0_ent1VSno0_qn_fpkm_sps" = "#033250", "no0_ent2VSno0_qn_fpkm_sps" = "#033250",
        "no0_ent1VSno0_qn_fpkm_sd" = "#033250", "no0_ent2VSno0_qn_fpkm_sd" = "#033250",
        "no0_ent1VSno0_qn_fpkm_cv" = "#033250", "no0_ent2VSno0_qn_fpkm_cv" = "#033250",
        "no0_adj_fpkmVSno0_ent1" = "#73b8d5", "no0_adj_fpkmVSno0_ent2" = "#73b8d5",
        "no0_adj_fpkm_sdVSno0_ent1" = "#73b8d5", "no0_adj_fpkm_sdVSno0_ent2" = "#73b8d5",
        "no0_adj_fpkm_cvVSno0_ent1" = "#73b8d5", "no0_adj_fpkm_cvVSno0_ent2" = "#73b8d5",
        "no0_ent1VSno0_tau" = "#73b8d5", "no0_ent2VSno0_tau" = "#73b8d5",
        "no0_ent2VSno0_ent1" = "#c01d2e"
    )) +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 45)
    ) +
    geom_abline(intercept = c(-0.5, 0, 0.5), slope = 0, linetype = "dashed", color = "black") +
    ylim(-1, 1) +
    ggtitle("spearman test of CRE entropies and \nassociated gene express (adjusted FPKM) \nacross genes in a biosample") +
    xlab("") +
    ylab("rho")
ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_no0_", dis, ".pdf"), plot = ent_rho_no0_p, width = 12, height = 12, dpi = 600)


sum_rho_spm <- data.frame()
sum_p_spm <- data.frame()
for (sp in species_ts) {
    for (tis in c("Brain", "Kidney", "Liver")) {
        tmp_spm <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/", sp, "_", tis, "_chiq.csv"),
            header = TRUE
        )
        tmp_spm$term <- paste(tmp_spm$Var1, tmp_spm$Var2, sep = "VS")
        tmp <- melt(tmp_spm[c("rho", "p", "term")], id.vars = "term", variable.names = "stat_name", value.name = "stat")
        tmp <- tmp[order(tmp$term), ]
        # tmp$name <- paste(tmp$term, tmp$variable, sep = "_")
        rho_tmp <- as.data.frame(t(tmp[tmp$variable == "rho", c("term", "stat")]))
        p_tmp <- tmp[tmp$variable == "p", c("term", "stat")] %>%
            t() %>%
            as.data.frame()
        if (TRUE %in% (grepl("....b_fpkm_diff", rho_tmp[1, ]))) {
            rho_tmp[1, ] <- gsub("....b_fpkm_diff", "", rho_tmp[1, ])
        }
        if (TRUE %in% (grepl("....b_fpkm_diff", p_tmp[1, ]))) {
            p_tmp[1, ] <- gsub("....b_fpkm_diff", "", p_tmp[1, ])
        }
        if (TRUE %in% (grepl("....k_fpkm_diff", rho_tmp[1, ]))) {
            rho_tmp[1, ] <- gsub("....k_fpkm_diff", "", rho_tmp[1, ])
        }
        if (TRUE %in% (grepl("....k_fpkm_diff", p_tmp[1, ]))) {
            p_tmp[1, ] <- gsub("....k_fpkm_diff", "", p_tmp[1, ])
        }
        if (TRUE %in% (grepl("....l_fpkm_diff", rho_tmp[1, ]))) {
            rho_tmp[1, ] <- gsub("....l_fpkm_diff", "", rho_tmp[1, ])
        }
        if (TRUE %in% (grepl("....l_fpkm_diff", p_tmp[1, ]))) {
            p_tmp[1, ] <- gsub("....l_fpkm_diff", "", p_tmp[1, ])
        }
        colnames(rho_tmp) <- rho_tmp[1, ]
        colnames(p_tmp) <- p_tmp[1, ]
        sum_rho_spm <- rbind(sum_rho_spm, cbind(rho_tmp, sp = sp, tis = tis)[-1, ])
        sum_p_spm <- rbind(sum_p_spm, cbind(p_tmp, sp = sp, tis = tis)[-1, ])
    }
}
rownames(sum_rho_spm) <- 1:nrow(sum_rho_spm)
rownames(sum_p_spm) <- 1:nrow(sum_p_spm)
sum_rho_spm[, 1:632] <- lapply(sum_rho_spm[, 1:632], as.numeric)
sum_p_spm[, 1:632] <- lapply(sum_p_spm[, 1:632], as.numeric)
sum <- data.frame(term = rep(colnames(sum_rho_spm)[1:632], 63), rho = 0, p = 0)
for (term in colnames(sum_rho_spm)[1:632]) {
    if (median(sum_rho_spm[, term]) > 0.5 & max(sum_p_spm[, term]) < 0.05) {
        print(term)
    }
    # p_tmp <- data.frame(p = sum_p_spm[, term], rho = sum_rho_spm[, term])
    # p <- ggplot(p_tmp, aes(x = rho, y = p)) +
    #     geom_point() +
    #     theme_bw()
    sum[sum$term == term, "rho"] <- sum_rho_spm[, term]
    sum[sum$term == term, "p"] <- sum_p_spm[, term]
}
plot_sum <- sum[sum$term %in% c("enhancer_numVSfpkm", "promoter_numVSfpkm", "median_enh_fcVSfpkm", "median_pro_fcVSfpkm"), ]
plot_sum$term <- gsub("VSfpkm", "", plot_sum$term)
plot_sum$term <- factor(plot_sum$term, levels = c("enhancer_num", "promoter_num", "median_enh_fc", "median_pro_fc"))
rho_box_p <- ggplot(plot_sum, aes(x = term, y = rho, fill = term)) +
    # geom_violin() +
    geom_boxplot() +
    scale_fill_manual(values = c("enhancer_num" = "#ed4343", "promoter_num" = "#73b8d5", "median_enh_fc" = "#c01d2e", "median_pro_fc" = "#033250")) +
    theme_bw() +
    theme() +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank()
    ) +
    theme(axis.text.x = element_text(hjust = 0.5)) +
    xlab("") +
    ylab("rho") +
    ggtitle("spearman test of CRE number and \nassociated gene express (raw FPKM) \namong biosamples")
ggsave("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/cre_fpkm.pdf", plot = rho_box_p, width = 10, height = 10, dpi = 600)

num_spm <- sum[sum$term %in% c("promoter_numVSenhancer_num", "median_pro_fcVSmedian_enh_fc", "median_pro"), ]
num_spm$term <- factor(num_spm$term, levels = c("promoter_numVSenhancer_num", "median_pro_fcVSmedian_enh_fc"))
num_spm_p <- ggplot(num_spm, aes(x = term, y = rho, fill = term)) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = c("promoter_numVSenhancer_num" = "#ed4343", "median_pro_fcVSmedian_enh_fc" = "#73b8d5")) +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank()
    ) +
    ggtitle("spearman test of CRE number and \n max foldchange among biosamples") +
    xlab("")
ggsave("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/cre_num_fc.pdf", plot = num_spm_p, width = 10, height = 10, dpi = 600)

ggplot(sum_rho_spm, aes(x = colnames(sum_rho_spm))) +
    geom_boxplot() +
    facet_wrap(~tis) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

example_point_df <- data.frame()
for (tis in c("Brain", "Kidney", "Liver")) {
    tmp <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/", sp, "_", tis, "_per_gene.csv"),
        header = TRUE
    )
    tmp$tissue <- tis
    example_point_df <- rbind(example_point_df, tmp)
}
example_point_plot <- ggplot(example_point_df, aes(x = enhancer_num, y = log(fpkm))) +
    geom_bin2d() +
    theme_bw() +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        panel.grid = element_blank()
    ) +
    facet_wrap(~tissue) +
    xlab("enhancer number") +
    ylab("log(raw FPKM)") +
    ggtitle("enhancer number and gene expression")
example_enh_num_fpkm_p <- ggplot(example_point_df, aes(x = enhancer_num, y = log(fpkm))) +
    geom_density_2d_filled()
ggsave("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/enh_num_fpkm.pdf", plot = example_enh_num_fpkm_p, width = 10, height = 10, dpi = 600)
example_pro_num_fpkm_p <- ggplot(example_point_df, aes(x = promoter_num, y = log(fpkm))) +
    geom_density_2d_filled()
ggsave("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/pro_num_fpkm.pdf", plot = example_pro_num_fpkm_p, width = 10, height = 10, dpi = 600)
example_point_df$log_fpkm <- log2(example_point_df$fpkm + 1)
cor.test(example_point_df$enhancer_num[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Brain"], example_point_df$log_fpkm[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Brain"], method = "spearman")
cor(example_point_df$promoter_num[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Brain"], example_point_df$log_fpkm[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Brain"], method = "spearman")
cor(example_point_df$enhancer_num[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Kidney"], example_point_df$log_fpkm[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Kidney"], method = "spearman")
cor(example_point_df$promoter_num[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Kidney"], example_point_df$log_fpkm[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Kidney"], method = "spearman")
cor(example_point_df$enhancer_num[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Liver"], example_point_df$log_fpkm[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Liver"], method = "spearman")
cor(example_point_df$promoter_num[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Liver"], example_point_df$log_fpkm[example_point_df$log_fpkm > 1 & example_point_df$tissue == "Liver"], method = "spearman")
cor.test(example_point_df$enhancer_num[example_point_df$tissue == "Brain"], example_point_df$fpkm[example_point_df$tissue == "Brain"], method = "spearman")


align_con_spm <- sum[sum$term %in% c("max_enh_conVSmax_enh_align", "max_pro_alignVSmax_enh_align", "mean_pro_conVSmean_pro_align", "median_pro_conVSmedian_pro_align"), ]
unique(align_con_spm$term)
align_con_spm$term <- factor(align_con_spm$term, levels = c("max_enh_conVSmax_enh_align", "max_pro_alignVSmax_enh_align", "mean_pro_conVSmean_pro_align", "median_pro_conVSmedian_pro_align"))
align_con_p <- ggplot(align_con_spm, aes(x = term, y = rho, fill = term)) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = c("max_enh_conVSmax_enh_align" = "#033250", "max_pro_alignVSmax_enh_align" = "#73b8d5", "mean_pro_conVSmean_pro_align" = "#c01d2e", "median_pro_conVSmedian_pro_align" = "#ed4343")) +
    theme(
        axis.text = element_text(size = 20, family = "Times"),
        axis.title = element_text(size = 25, family = "Times"),
        title = element_text(size = 30, family = "Times"),
        legend.title = element_text(size = 20, family = "Times"),
        legend.text = element_text(size = 20, family = "Times"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 18)
    ) +
    ylim(0, 1) +
    scale_x_discrete(labels = c(
        "max_enh_conVSmax_enh_align" = "max EC & GC\n of enhancers", "max_pro_alignVSmax_enh_align" = "max GC of \nenhancers and promoters",
        "mean_pro_conVSmean_pro_align" = "mean EC &GC \nof promoters", "median_pro_conVSmedian_pro_align" = "median EC & GC \n of promoters"
    )) +
    ggtitle("spearman rhos of CRE \nEC (Num of species epigenetically conserved) and \nGC (Num of species genetically conserved) \nof different biosamples") +
    xlab("")
ggsave("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/cre_con_align.pdf", plot = align_con_p, width = 12, height = 12, dpi = 600)

# example_ec_gc_df <- data.frame()
# for (tis in c("Brain", "Kidney", "Liver")) {
#     tmp <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/", sp, "_", tis, "_per_gene.csv"),
#         header = TRUE
#     )
#     tmp$tissue <- tis
#     example_ec_gc_df <- rbind(example_ec_gc_df, tmp)
# }
sp <- "Petaurus_breviceps"
example_ec_density_p <- ggplot(example_point_df, aes(x = max_enh_align, y = max_enh_con)) +
    geom_density_2d_filled()
example_enh_pro_gc_density_p <- ggplot(example_point_df[example_point_df$max_enh_align != 0 & example_point_df$max_pro_align != 0, ], aes(x = max_pro_align, y = max_enh_align)) +
    geom_density_2d_filled()
cor.test(example_point_df$max_enh_align[example_point_df$tissue == "Brain" & example_point_df$max_enh_align != 0 & example_point_df$max_pro_align != 0], example_point_df$max_pro_align[example_point_df$tissue == "Brain" & example_point_df$max_enh_align != 0 & example_point_df$max_pro_align != 0], method = "spearman")

example_mid_pro_con_align_p <- ggplot(
    example_point_df[example_point_df$median_pro_align > 0 & example_point_df$median_pro_con > 0, ],
    aes(x = median_pro_align, y = median_pro_con)
) +
    geom_density_2d_filled()
cor.test(example_point_df[example_point_df$median_pro_align > 0 & example_point_df$median_pro_con > 0, "median_pro_align"], example_point_df[example_point_df$median_pro_align > 0 & example_point_df$median_pro_con > 0, "median_pro_con"], method = "spearman")

for (dis in c("1M", "500k", "200k", "100k", "50k")) {
    ent_rho_spm <- data.frame()
    ent_p_spm <- data.frame()
    for (sp in species_ts) {
        # for (tis in c("Brain", "Kidney", "Liver")) {
        tmp_spm <- try(read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/", sp, "_", dis, "_filtered_ent_cor.csv"),
            header = TRUE
        ))
        if (class(tmp_spm) == "try-error") {
            next
        }
        tmp_spm$term <- paste(tmp_spm$Var1, tmp_spm$Var2, sep = "VS")
        tmp <- melt(tmp_spm[c("rho", "p", "term")], id.vars = "term", variable.names = "stat_name", value.name = "stat")
        tmp <- tmp[order(tmp$term), ]
        # tmp$name <- paste(tmp$term, tmp$variable, sep = "_")
        rho_tmp <- as.data.frame(t(tmp[tmp$variable == "rho", c("term", "stat")]))
        p_tmp <- tmp[tmp$variable == "p", c("term", "stat")] %>%
            t() %>%
            as.data.frame()
        colnames(rho_tmp) <- rho_tmp[1, ]
        colnames(p_tmp) <- p_tmp[1, ]
        ent_rho_spm <- rbind(ent_rho_spm, cbind(rho_tmp, sp = sp)[-1, ])
        ent_p_spm <- rbind(ent_p_spm, cbind(p_tmp, sp = sp)[-1, ])
        # }
    }
    rownames(ent_rho_spm) <- 1:nrow(ent_rho_spm)
    rownames(ent_p_spm) <- 1:nrow(ent_p_spm)
    ent_rho_spm[, 1:45] <- lapply(ent_rho_spm[, 1:45], as.numeric)
    ent_p_spm[, 1:45] <- lapply(ent_p_spm[, 1:45], as.numeric)
    ent_sum <- data.frame(term = rep(colnames(ent_rho_spm)[1:45], nrow(ent_rho_spm)), rho = 0, p = 0)
    terms <- c()
    for (term in colnames(ent_rho_spm)[1:45]) {
        if (median(ent_rho_spm[, term]) > 0.5 & max(ent_p_spm[, term]) < 0.05) {
            print(term)
            terms <- c(terms, term)
        }
        ent_sum[ent_sum$term == term, "rho"] <- ent_rho_spm[, term]
        ent_sum[ent_sum$term == term, "p"] <- ent_p_spm[, term]
    }
    terms <- c(
        "ent1VSfpkm", "ent2VSfpkm", "ent1VSfpkm_sd", "ent2VSfpkm_sd",
        "ent1VSfpkm_se", "ent2VSfpkm_se",
        "ent1VSqn_fpkm_sps", "ent2VSqn_fpkm_sps",
        "ent1VSqn_fpkm_sd", "ent2VSqn_fpkm_sd",
        "ent1VSqn_fpkm_se", "ent2VSqn_fpkm_se",
        "ent1VStau", "ent2VStau", "ent2VSent1"
    )
    ent_rho_p_df <- ent_sum[ent_sum$term %in% terms, ]
    ent_rho_p_df$term <- factor(ent_rho_p_df$term, levels = c(
        "ent2VSent1",
        "ent1VSfpkm", "ent2VSfpkm",
        "ent1VSfpkm_sd", "ent2VSfpkm_sd",
        "ent1VSfpkm_se", "ent2VSfpkm_se",
        "ent1VSqn_fpkm_sps", "ent2VSqn_fpkm_sps",
        "ent1VSqn_fpkm_sd", "ent2VSqn_fpkm_sd",
        "ent1VSqn_fpkm_se", "ent2VSqn_fpkm_se",
        "ent1VStau", "ent2VStau"
    ))
    ent_rho_p <- ggplot(ent_rho_p_df, aes(x = term, y = rho, fill = term)) +
        geom_boxplot() +
        theme_bw() +
        scale_fill_manual(values = c(
            "ent2VSent1" = "#c01d2e",
            "ent1VSfpkm" = "#ed4343", "ent2VSfpkm" = "#ed4343",
            "ent1VSfpkm_sd" = "#ed4343", "ent2VSfpkm_sd" = "#ed4343",
            "ent1VSfpkm_se" = "#ed4343", "ent2VSfpkm_se" = "#ed4343",
            "ent1VSqn_fpkm_sps" = "#033250", "ent2VSqn_fpkm_sps" = "#033250",
            "ent1VSqn_fpkm_sd" = "#033250", "ent2VSqn_fpkm_sd" = "#033250",
            "ent1VSqn_fpkm_se" = "#033250", "ent2VSqn_fpkm_se" = "#033250",
            "ent1VStau" = "#73b8d5", "ent2VStau" = "#73b8d5"
        )) +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times"),
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(hjust = 1, angle = 45)
        ) +
        geom_abline(intercept = c(-0.5, 0, 0.5), slope = 0, linetype = "dashed", color = "black") +
        ylim(-1, 1) +
        ggtitle(paste0("spearman test of CRE entropies and \nassociated gene express (raw FPKM) \nacross genes in a biosample (", dis, ")")) +
        xlab("") +
        ylab("rho")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_", dis, ".pdf"), plot = ent_rho_p, width = 12, height = 12, dpi = 600)
}
for (dis in c("1M", "500k", "200k", "100k", "50k")) {
    ent_rho_spm_no0 <- data.frame()
    ent_p_spm_no0 <- data.frame()
    for (sp in species_ts) {
        tmp_spm <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/", sp, "_", dis, "_ent_cor_no0_adj.csv"),
            header = TRUE
        )
        tmp_spm$term <- paste(tmp_spm$Var1, tmp_spm$Var2, sep = "VS")
        tmp <- melt(tmp_spm[c("rho", "p", "term")], id.vars = "term", variable.names = "stat_name", value.name = "stat")
        tmp <- tmp[order(tmp$term), ]
        # tmp$name <- paste(tmp$term, tmp$variable, sep = "_")
        rho_tmp <- as.data.frame(t(tmp[tmp$variable == "rho", c("term", "stat")]))
        p_tmp <- tmp[tmp$variable == "p", c("term", "stat")] %>%
            t() %>%
            as.data.frame()
        colnames(rho_tmp) <- rho_tmp[1, ]
        colnames(p_tmp) <- p_tmp[1, ]
        ent_rho_spm_no0 <- rbind(ent_rho_spm_no0, cbind(rho_tmp, sp = sp)[-1, ])
        ent_p_spm_no0 <- rbind(ent_p_spm_no0, cbind(p_tmp, sp = sp)[-1, ])
    }
    rownames(ent_rho_spm_no0) <- 1:nrow(ent_rho_spm_no0)
    rownames(ent_p_spm_no0) <- 1:nrow(ent_p_spm_no0)
    ent_rho_spm_no0[, 1:78] <- lapply(ent_rho_spm_no0[, 1:78], as.numeric)
    ent_p_spm_no0[, 1:78] <- lapply(ent_p_spm_no0[, 1:78], as.numeric)
    ent_sum_no0 <- data.frame(term = rep(colnames(ent_rho_spm_no0)[1:78], 21), rho = 0, p = 0)
    ent_no0_terms <- c()
    for (term in colnames(ent_rho_spm_no0)[1:78]) {
        if (median(ent_rho_spm_no0[, term]) > 0.5 & max(ent_p_spm_no0[, term]) < 0.05) {
            print(term)
            ent_no0_terms <- c(ent_no0_terms, term)
        }
        ent_sum_no0[ent_sum_no0$term == term, "rho"] <- ent_rho_spm_no0[, term]
        ent_sum_no0[ent_sum_no0$term == term, "p"] <- ent_p_spm_no0[, term]
    }
    terms <- c(
        "no0_ent1VSno0_fpkm", "no0_ent2VSno0_fpkm", "no0_ent1VSno0_fpkm_sd",
        "no0_ent2VSno0_fpkm_sd", "no0_ent1VSno0_fpkm_cv", "no0_ent2VSno0_fpkm_cv",
        "no0_ent1VSno0_qn_fpkm_sps", "no0_ent2VSno0_qn_fpkm_sps", "no0_ent1VSno0_qn_fpkm_sd",
        "no0_ent2VSno0_qn_fpkm_sd", "no0_ent1VSno0_qn_fpkm_cv", "no0_ent2VSno0_qn_fpkm_cv",
        "no0_adj_fpkmVSno0_ent1", "no0_adj_fpkmVSno0_ent2", "no0_adj_fpkm_sdVSno0_ent1",
        "no0_adj_fpkm_sdVSno0_ent2", "no0_adj_fpkm_cvVSno0_ent1", "no0_adj_fpkm_cvVSno0_ent2",
        "no0_ent1VSno0_tau", "no0_ent2VSno0_tau", "no0_ent2VSno0_ent1"
    )
    ent_rho_spm_no0_df <- ent_sum_no0[ent_sum_no0$term %in% terms, ]
    ent_rho_spm_no0_df$term <- factor(ent_rho_spm_no0_df$term, levels = c(
        "no0_ent1VSno0_fpkm", "no0_ent2VSno0_fpkm",
        "no0_ent1VSno0_fpkm_sd", "no0_ent2VSno0_fpkm_sd",
        "no0_ent1VSno0_fpkm_cv", "no0_ent2VSno0_fpkm_cv",
        "no0_ent1VSno0_qn_fpkm_sps", "no0_ent2VSno0_qn_fpkm_sps",
        "no0_ent1VSno0_qn_fpkm_sd", "no0_ent2VSno0_qn_fpkm_sd",
        "no0_ent1VSno0_qn_fpkm_cv", "no0_ent2VSno0_qn_fpkm_cv",
        "no0_adj_fpkmVSno0_ent1", "no0_adj_fpkmVSno0_ent2",
        "no0_adj_fpkm_sdVSno0_ent1", "no0_adj_fpkm_sdVSno0_ent2",
        "no0_adj_fpkm_cvVSno0_ent1", "no0_adj_fpkm_cvVSno0_ent2",
        "no0_ent1VSno0_tau", "no0_ent2VSno0_tau", "no0_ent2VSno0_ent1"
    ))
    ent_rho_no0_p <- ggplot(ent_rho_spm_no0_df, aes(x = term, y = rho, fill = term)) +
        geom_boxplot() +
        theme_bw() +
        scale_fill_manual(values = c(
            "no0_ent1VSno0_fpkm" = "#ed4343", "no0_ent2VSno0_fpkm" = "#ed4343",
            "no0_ent1VSno0_fpkm_sd" = "#ed4343", "no0_ent2VSno0_fpkm_sd" = "#ed4343",
            "no0_ent1VSno0_fpkm_cv" = "#ed4343", "no0_ent2VSno0_fpkm_cv" = "#ed4343",
            "no0_ent1VSno0_qn_fpkm_sps" = "#033250", "no0_ent2VSno0_qn_fpkm_sps" = "#033250",
            "no0_ent1VSno0_qn_fpkm_sd" = "#033250", "no0_ent2VSno0_qn_fpkm_sd" = "#033250",
            "no0_ent1VSno0_qn_fpkm_cv" = "#033250", "no0_ent2VSno0_qn_fpkm_cv" = "#033250",
            "no0_adj_fpkmVSno0_ent1" = "#73b8d5", "no0_adj_fpkmVSno0_ent2" = "#73b8d5",
            "no0_adj_fpkm_sdVSno0_ent1" = "#73b8d5", "no0_adj_fpkm_sdVSno0_ent2" = "#73b8d5",
            "no0_adj_fpkm_cvVSno0_ent1" = "#73b8d5", "no0_adj_fpkm_cvVSno0_ent2" = "#73b8d5",
            "no0_ent1VSno0_tau" = "#73b8d5", "no0_ent2VSno0_tau" = "#73b8d5",
            "no0_ent2VSno0_ent1" = "#c01d2e"
        )) +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times"),
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(hjust = 1, angle = 45)
        ) +
        geom_abline(intercept = c(-0.5, 0, 0.5), slope = 0, linetype = "dashed", color = "black") +
        ylim(-1, 1) +
        ggtitle("spearman test of CRE entropies and \nassociated gene express (adjusted FPKM) \nacross genes in a biosample") +
        xlab("") +
        ylab("rho")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_fpkm_no0_", dis, ".pdf"), plot = ent_rho_no0_p, width = 12, height = 12, dpi = 600)
}

for (dis in c("1M", "500k", "200k", "100k", "50k")) {
    ent_rho_spm_adj <- data.frame()
    ent_p_spm_adj <- data.frame()
    for (sp in species_ts) {
        # for (tis in c("Brain", "Kidney", "Liver")) {
        tmp_spm <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/", sp, "_", dis, "_ent_cor_adj.csv"),
            header = TRUE
        )
        tmp_spm$term <- paste(tmp_spm$Var1, tmp_spm$Var2, sep = "VS")
        tmp <- melt(tmp_spm[c("rho", "p", "term")], id.vars = "term", variable.names = "stat_name", value.name = "stat")
        tmp <- tmp[order(tmp$term), ]
        # tmp$name <- paste(tmp$term, tmp$variable, sep = "_")
        rho_tmp <- as.data.frame(t(tmp[tmp$variable == "rho", c("term", "stat")]))
        p_tmp <- tmp[tmp$variable == "p", c("term", "stat")] %>%
            t() %>%
            as.data.frame()
        colnames(rho_tmp) <- rho_tmp[1, ]
        colnames(p_tmp) <- p_tmp[1, ]
        ent_rho_spm_adj <- rbind(ent_rho_spm_adj, cbind(rho_tmp, sp = sp)[-1, ])
        ent_p_spm_adj <- rbind(ent_p_spm_adj, cbind(p_tmp, sp = sp)[-1, ])
        # }
    }
    rownames(ent_rho_spm_adj) <- 1:nrow(ent_rho_spm_adj)
    rownames(ent_p_spm_adj) <- 1:nrow(ent_p_spm_adj)
    ent_rho_spm_adj[, 1:78] <- lapply(ent_rho_spm_adj[, 1:78], as.numeric)
    ent_p_spm_adj[, 1:78] <- lapply(ent_p_spm_adj[, 1:78], as.numeric)
    ent_sum_adj <- data.frame(term = rep(colnames(ent_rho_spm_adj)[1:78], 21), rho = 0, p = 0)
    ent_adj_terms <- c()
    for (term in colnames(ent_rho_spm_adj)[1:78]) {
        if (median(ent_rho_spm_adj[, term]) > 0.5 & max(ent_p_spm_adj[, term]) < 0.05) {
            print(term)
            ent_adj_terms <- c(ent_adj_terms, term)
        }
        ent_sum_adj[ent_sum_adj$term == term, "rho"] <- ent_rho_spm_adj[, term]
        ent_sum_adj[ent_sum_adj$term == term, "p"] <- ent_p_spm_adj[, term]
    }
    terms <- c("adj_fpkmVSent1", "adj_fpkmVSent2", "adj_fpkm_sdVSent1", "adj_fpkm_sdVSent2", "adj_fpkm_cvVSent1", "adj_fpkm_cvVSent2")
    ent_rho_spm_adj_df <- ent_sum_adj[ent_sum_adj$term %in% terms, ]
    ent_rho_spm_adj_df$term <- factor(ent_rho_spm_adj_df$term, levels = c("adj_fpkmVSent1", "adj_fpkmVSent2", "adj_fpkm_sdVSent1", "adj_fpkm_sdVSent2", "adj_fpkm_cvVSent1", "adj_fpkm_cvVSent2"))
    ent_rho_adj_p <- ggplot(ent_rho_spm_adj_df, aes(x = term, y = rho, fill = term)) +
        geom_boxplot() +
        theme_bw() +
        scale_fill_manual(values = c("adj_fpkmVSent1" = "#ed4343", "adj_fpkmVSent2" = "#ed4343", "adj_fpkm_sdVSent1" = "#ed4343", "adj_fpkm_sdVSent2" = "#ed4343", "adj_fpkm_cvVSent1" = "#ed4343", "adj_fpkm_cvVSent2" = "#ed4343")) +
        theme(
            axis.text = element_text(size = 20, family = "Times"),
            axis.title = element_text(size = 25, family = "Times"),
            title = element_text(size = 30, family = "Times"),
            legend.title = element_text(size = 20, family = "Times"),
            legend.text = element_text(size = 20, family = "Times"),
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(hjust = 1, angle = 45)
        ) +
        geom_abline(intercept = c(-0.5, 0, 0.5), slope = 0, linetype = "dashed", color = "black") +
        ylim(-1, 1) +
        ggtitle(paste0("spearman test of CRE entropies and \nassociated gene express (adjusted FPKM) \nacross genes in a biosample (", dis, ")")) +
        xlab("") +
        ylab("rho")
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/cre_ent_adj_fpkm_", dis, ".pdf"), plot = ent_rho_adj_p, width = 12, height = 12, dpi = 600)
}
in_sp_rho_spm <- data.frame()
in_sp_p_spm <- data.frame()
for (sp in species_ts) {
    for (tis in c("Brain", "Kidney", "Liver")) {
        tmp_spm <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/in_sp/", sp, "_", tis, "_in_cor.csv"),
            header = TRUE
        )
        tmp_spm$term <- paste(tmp_spm$Var1, tmp_spm$Var2, sep = "VS")
        tmp <- melt(tmp_spm[c("rho", "p", "term")], id.vars = "term", variable.names = "stat_name", value.name = "stat")
        tmp <- tmp[order(tmp$term), ]
        # tmp$name <- paste(tmp$term, tmp$variable, sep = "_")
        rho_tmp <- as.data.frame(t(tmp[tmp$variable == "rho", c("term", "stat")]))
        p_tmp <- tmp[tmp$variable == "p", c("term", "stat")] %>%
            t() %>%
            as.data.frame()
        colnames(rho_tmp) <- rho_tmp[1, ]
        colnames(p_tmp) <- p_tmp[1, ]
        in_sp_rho_spm <- rbind(in_sp_rho_spm, cbind(rho_tmp, sp = sp, tis = tis)[-1, ])
        in_sp_p_spm <- rbind(in_sp_p_spm, cbind(p_tmp, sp = sp, tis = tis)[-1, ])
    }
}

save.image("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/cor_each_sp_sum.RData")
