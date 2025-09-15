# assess alignment coverage of all species, to see if the assembly of each species is good
library(ggtree)
library(ggplot2)
library(ape)
library(nlme)
library(rr2)
library(tidyr)
library(plyr)
library(dplyr)

library(doMC)
doMC::registerDoMC(cores = 4)

tre <- read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
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
species2 <- c(
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
order_colors <- c(
    "Artiodactyla" = "#3682be", "Carnivora" = "#45a776", "Perissodactyla" = "#f05330",
    "Chiroptera" = "#eed777", "Eulipotyphla" = "#38cb7d", "Rodentia" = "#334f65",
    "Lagomorpha" = "#ddae33", "Primates" = "#b3974e", "Scandentia" = "#844bb3",
    "Hyracoidea" = "#93c555", "Diprotodontia" = "#5f6694"
)


tre$tip.label <- gsub("Macaca_mulatta", "Macaca_fascicularis", tre$tip.label)
tre <- ape::keep.tip(tre, species)
setwd("/media/Data/zhangz/chip/genomes/halcoverage")
align_cov <- read.csv("/media/Data/zhangz/chip/genomes/halcoverage/coverage.log", header = TRUE, sep = ",")

# sum row
align_cov$sum <- rowSums(align_cov[, -1])
align_cov$coverage <- align_cov$sum / max(align_cov$sum)
align_cov <- align_cov[align_cov$Genome %in% species, ]
div_time <- cophenetic.phylo(tre)
for (sp in species) {
    align_cov[align_cov$Genome == sp, "div_time"] <- div_time[sp, "Mus_musculus"] / 2
}
sp_info <- read.csv("/media/Data/zhangz/chip/scripts/info/species.csv", header = TRUE)
sp_info$species <- gsub("Macaca_mulatta", "Macaca_fascicularis", sp_info$species)

div_coverage_gls <- gls(coverage ~ div_time, data = align_cov, method = "ML")
div_coverage_gls
div_coverage_r2 <- R2_lik(mod = div_coverage_gls)
div_coverage_point <- ggplot(align_cov[align_cov$Genome != "Mus_musculus", ], aes(x = div_time, y = coverage)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() +
    xlab("Species divergence time") +
    ylab("Genome alignment coverage") +
    theme(text = element_text(size = 8), se = TRUE) +
    annotate("text", x = 80, y = 0.7, label = paste(
        "R^2 = ", round(div_coverage_r2, 2),
        ", p value =", round(summary(div_coverage_gls)$tTable[2, 4], 2)
    ), size = 3)

ggsave("/media/Data/zhangz/chip/analysis/summary2/align_coverage_div.pdf", div_coverage_point, width = 7, height = 7, units = "cm")

assembly <- read.csv("/media/Data/zhangz/chip/genomes/quality/report.tsv", header = TRUE, sep = "\t")
for (sp in species) {
    align_cov[align_cov$Genome == sp, "N50"] <- assembly[assembly$Assembly == "N50", sp]
    align_cov[align_cov$Genome == sp, "N_per_100_kbp"] <- assembly[assembly$Assembly == "# N's per 100 kbp", sp]
}
coverage_n50_gls_ou <- gls(coverage ~ N50, data = align_cov, method = "ML", correlation = corMartins(1, phy = tre))
coverage_n50_gls_ou_r2 <- R2_lik(mod = coverage_n50_gls_ou)
coverage_n50_gls_bm <- gls(coverage ~ N50, data = align_cov, method = "ML", correlation = corBrownian(1, phy = tre))
coverage_n50_gls_bm_r2 <- R2_lik(mod = coverage_n50_gls_bm)
cor.test(pic(align_cov$coverage, phy = tre), pic(align_cov$N50, phy = tre), method = "spearman")
coverage_n_gls_ou <- gls(coverage ~ N_per_100_kbp, data = align_cov, method = "ML", correlation = corMartins(1, phy = tre))
coverage_n_gls_ou_r2 <- R2_lik(mod = coverage_n_gls_ou)
coverage_n_gls_bm <- gls(coverage ~ N_per_100_kbp, data = align_cov, method = "ML", correlation = corBrownian(1, phy = tre))
coverage_n_gls_bm_r2 <- R2_lik(mod = coverage_n_gls_bm)

cor.test(pic(align_cov$coverage, phy = tre), pic(align_cov$N_per_100_kbp, phy = tre), method = "spearman")
cor.test(pic(align_cov$coverage, phy = tre), pic(align_cov$N_per_100_kbp, phy = tre), method = "pearson")
cor.test(align_cov$coverage, align_cov$N_per_100_kbp, method = "spearman")

# correlation heat map of coverage, N50, N_per_100_kbp
pic_matrix <- data.frame(coverage = pic(align_cov$coverage, phy = tre), N50 = pic(align_cov$N50, phy = tre), N_per_100_kbp = pic(align_cov$N_per_100_kbp, phy = tre))
cor_matrix <- cor(pic_matrix, use = "pairwise.complete.obs", method = "spearman")
library(Hmisc)
corr_p <- rcorr(as.matrix(pic_matrix), type = "spearman")
library(corrplot)
pdf("/media/Data/zhangz/chip/analysis/summary2/align_coverage_cor_heatmap.pdf")
corrplot(cor_matrix, method = "ellipse", type = "upper", tl.col = "black", tl.cex = 0.8, tl.srt = 45, tl.pos = "lt", p.mat = corr_p$P, sig.level = 0.05)
corrplot(cor_matrix, method = "number", type = "lower", add = TRUE, tl.col = "n", tl.cex = 0.8, tl.pos = "n", p.mat = corr_p$P, sig.level = 0.05)
dev.off()


# functions
se <- function(x, model) {
    ## calculate standard error of the slope, x is predictor variable
    vi <- vcov(model)[1, 1] + x * vcov(model)[1, 2] * 1 + (1 * vcov(model)[2, 1] + x * vcov(model)[2, 2]) * x
    se <- sqrt(vi)
    return(se)
}

pgls_pmax_model <- function(var1, var2, data_df, tree, model, lg = FALSE) {
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

PGLS <- function(var1, var2, data_df, lg1 = FALSE, lg2 = FALSE) {
    # fit pgls model to data_df on var1 and var2
    data_df <- data_df[, c("species", var1, var2, "order", "svg_path")]
    data_df <- data_df[!is.na(data_df[[var1]]) & !is.na(data_df[[var2]]), ]
    if (lg1 == TRUE) {
        data_df[[var1]] <- log(data_df[[var1]])
    }
    if (lg2 == TRUE) {
        data_df[[var2]] <- log(data_df[[var2]])
    }
    spp <- data_df$species
    tree <- ape::keep.tip(tre, spp)
    formula_str <- paste(var2, " ~ ", var1, sep = "")
    bm_model <- gls(as.formula(formula_str), correlation = corBrownian(phy = tree), data = data_df, method = "ML", na.action = na.omit)
    ou_model <- try(gls(as.formula(formula_str), correlation = corMartins(1, phy = tree), data = data_df, method = "ML", na.action = na.omit))
    if (class(ou_model) == "try-error") {
        ou_model <- gls(as.formula(formula_str), correlation = corMartins(1, phy = tree, fixed = TRUE), data = data_df, method = "ML", na.action = na.omit)
    }
    models <- list(BM = bm_model, OU = ou_model)
    lrt <- anova(bm_model, ou_model)
    if (lrt$`p-value`[2] < 0.05) {
        model <- names(models)[which.max(lrt$logLik)]
    } else {
        model <- "OU"
    }
    pmax <- pgls_pmax_model(var1, var2, data_df, tree, model)
    res_model <- models[[model]]
    res_summary <- summary(res_model)
    res_r2 <- R2_lik(mod = res_model)
    res_pred <- data.frame(x = data_df[[var1]], pred = res_model$fitted)
    res_se <- se(data_df[[var1]], res_model)
    res_pred <- cbind(res_pred, res_se)
    colnames(res_pred) <- c("x", "pred", "se")
    res_pred$lwr <- res_pred$pred - 1.96 * res_pred$se
    res_pred$upr <- res_pred$pred + 1.96 * res_pred$se
    res_slope <- res_summary$tTable[2, 1]
    res_intercept <- res_summary$tTable[1, 1]
    res_p <- res_summary$tTable[2, 4]
    res_plot <- ggplot(data_df, aes_string(x = var1, y = var2)) +
        geom_point(aes(fill = order, color = order), size = 1, shape = 24) +
        geom_abline(intercept = res_intercept, slope = res_slope, color = "red") +
        geom_ribbon(data = res_pred, aes(x = x, y = pred, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) +
        scale_fill_manual(values = order_colors) +
        scale_color_manual(values = order_colors) +
        # annotate("text", label = paste(
        #     "R^2:", round(res_r2, 4),
        #     "p:", round(res_p, 4),
        #     "pmax:", round(as.numeric(pmax[1]), 4),
        #     "model:", model
        # ), size = 3, hjust = 0, vjust = 1) +
        theme_classic() +
        theme(
            axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12),
            axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.line = element_line(colour = "black"), legend.position = "none"
        ) +
        scale_x_continuous(n.breaks = 3, labels = scales::comma, expand = c(0.15, 0)) +
        scale_y_continuous(n.breaks = 3) +
        geom_image(aes(image = svg_path), size = 0.05, height = 0.03, alpha = 0.7) +
        labs(title = paste(formula_str, ", Model:", model, "\nP-value:", round(res_p, 4), ", max P-value:", round(as.numeric(pmax[1]), 4), "\nSpecies:", pmax[2], ", R2:", round(res_r2, 4)))
    return(list(model = model, pmax = pmax, res = res_model, plot = res_plot))
}

# div_coverage_pred <- data.frame(x = align_cov$div_time, pred = )
coverage_n_point <- ggplot(align_cov, aes(x = N_per_100_kbp, y = coverage)) +
    geom_point() +
    theme_classic() +
    xlab("Number of N's per 100 kbp") +
    ylab("Genome alignment coverage") +
    theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70")) +
    geom_smooth(method = "lm", se = FALSE, color = "red")
ggsave("/media/Data/zhangz/chip/analysis/summary2/align_coverage.pdf", coverage_n_point, width = 7, height = 7, units = "cm")

# align data of Mus_musculus
sp <- "Mus_musculus"


a_ply(species, 1, function(sp) {
    if (sp == "Macaca_fascicularis") {
        ele_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", "Macaca_mulatta", "_overlap_anno.csv"), header = TRUE, sep = ",")
    } else {
        ele_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    }
    # m_ply(expand.grid(tis = unique(ele_df$tissue), ele = c("enhancer", "promoter")), function(tis, ele) {
    conds <- expand.grid(tis = unique(ele_df$tissue), ele = c("enhancer", "promoter"))
    for (i in 1:nrow(conds)) {
        tis <- conds$tis[i]
        ele <- conds$ele[i]
        # print(paste(tis, ele))
        if (sp == "Macaca_fascicularis") {
            align_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/", "Macaca_mulatta", "/compare2/", tis, "/", ele, "/", "Macaca_mulatta", "_", tis, "_", ele, "_stats.tsv"), header = TRUE, sep = "\t")
        } else {
            align_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/", sp, "/compare2/", tis, "/", ele, "/", sp, "_", tis, "_", ele, "_stats.tsv"), header = TRUE, sep = "\t")
        }
        align_df$align22 <- rowSums(align_df[, species2[species != sp]])
        ele_df[which(ele_df$tissue == tis & ele_df$element == ele), "align22"] <- align_df$align22[match(ele_df[which(ele_df$tissue == tis & ele_df$element == ele), "peak"], align_df$peak)]
    }
    # }, .parallel = TRUE)
    ele_df$align22[is.na(ele_df$align22)] <- 0
    if (sp == "Macaca_fascicularis") {
        write.csv(ele_df, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", "Macaca_mulatta", "_overlap_align22.csv"), row.names = FALSE)
    } else {
        write.csv(ele_df, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_align22.csv"), row.names = FALSE)
    }
}, .parallel = TRUE)

a_ply(species, 1, function(sp) {
    if (sp == "Macaca_fascicularis") {
        ele_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", "Macaca_mulatta", "_overlap_align22.csv"), header = TRUE, sep = ",")
    } else {
        ele_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_align22.csv"), header = TRUE, sep = ",")
    }
    # draw distribution of alignable species number
    align_dis_enhancer_plot <- ggplot(ele_df[which(ele_df$align22 > 0 & ele_df$element == "enhancer"), ], aes(x = align22)) +
        geom_histogram(binwidth = 1, fill = "#ef4343", color = "black", linewidth = 0.05) +
        theme_classic() +
        scale_x_continuous(breaks = c(1, 10, 21)) +
        scale_y_continuous(n.breaks = 3) +
        xlab(paste0("Number of species that CREs sequence in\n", gsub("_", " ", sp), " can be aligned to")) +
        theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70"))
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_enhancer_align_dis.pdf"), align_dis_enhancer_plot, width = 7, height = 7, units = "cm")
    align_dis_promoter_plot <- ggplot(ele_df[which(ele_df$align22 > 0 & ele_df$element == "promoter"), ], aes(x = align22)) +
        geom_histogram(binwidth = 1, fill = "#73b8d5", color = "black", linewidth = 0.05) +
        theme_classic() +
        scale_x_continuous(breaks = c(1, 10, 21)) +
        scale_y_continuous(n.breaks = 3) +
        xlab(paste0("Number of species that CREs sequence in\n", gsub("_", " ", sp), " can be aligned to")) +
        theme(text = element_text(size = 8), axis.line = element_line(linewidth = 0.5, color = "grey70"))
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_promoter_align_dis.pdf"), align_dis_promoter_plot, width = 7, height = 7, units = "cm")
}, .parallel = TRUE)


# pgls plot of genome alignment coverage and N50
colnames(align_cov)[1] <- "species"
align_cov$order <- sp_info$order[match(align_cov$species, sp_info$species)]

align_cov$svg_path <- paste0("/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/", align_cov$species, ".svg")
align_cov$svg_path[which(align_cov$species %in% c("Myotis_chinensis", "Myotis_ricketti", "Rhinolophus_ferrumequinum", "Rhinolophus_pusillus"))] <- "/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/Hipposideros_larvatus.svg"
align_cov$svg_path[which(align_cov$species == "Macaca_fascicularis")] <- "/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/Macaca_mulatta.svg"
for (svg in align_cov$svg_path) {
    if (!file.exists(svg)) {
        print(svg)
    }
}

coverage_n50_pgls <- PGLS("N50", "coverage", align_cov)
ggsave("/media/Data/zhangz/chip/genomes/halcoverage/coverage_n50_point.pdf", coverage_n50_pgls$plot, height = 6, width = 6)

coverage_n_pgls <- PGLS("N_per_100_kbp", "coverage", align_cov)
ggsave("/media/Data/zhangz/chip/genomes/halcoverage/coverage_n_point.pdf", coverage_n_pgls$plot, height = 6, width = 6)
