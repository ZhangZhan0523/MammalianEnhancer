## use PGLS to find correlation between gene regulatory entropy and life history traits,
## for example, body mass, longevity, BMR, etc.

## load packages
library(ape)
library(caper)
library(phytools)
library(plyr)
library(doMC)
library(nlme)
library(rr2)
library(geiger)
library(ggrepel)
library(ggimage)
library(ggplot2)
library(magrittr)
library(stringr)
library(UpSetR)
library(data.table)
library(tidyr)
library(ggplotify)
library(rtracklayer)

doMC::registerDoMC(50)

## functions
se <- function(x, model) {
    ## calculate standard error of the slope, x is predictor variable
    vi <- vcov(model)[1, 1] + x * vcov(model)[1, 2] * 1 + (1 * vcov(model)[2, 1] + x * vcov(model)[2, 2]) * x
    se <- sqrt(vi)
    return(se)
}

pgls_pmax_model <- function(var1, var2, data_df, tree, model, lg1 = TRUE, lg2 = TRUE) {
    ## fit PGLS model and return the maximum p-value and the species with the maximum p-value
    data_df <- data_df[!is.na(data_df[[var1]]) & !is.na(data_df[[var2]]), ]
    spp <- data_df$species
    tree <- keep.tip(tree, spp)
    p_list <- adply(spp, 1, function(sp) {
        data_sp <- data_df[data_df$species != sp, ]
        tree_sp <- drop.tip(tree, sp)
        var1 <- ifelse(lg1 == TRUE, paste("log(", var1, ")", sep = ""), var1)
        var2 <- ifelse(lg2 == TRUE, paste("log(", var2, ")", sep = ""), var2)
        formula_str <- paste(var2, " ~ ", var1, sep = "")
        # if (lg == TRUE) {
        #     formula_str <- paste("log(", var2, ") ~ log(", var1, ")", sep = "")
        # } else {
        #     formula_str <- paste(var2, " ~ ", var1, sep = "")
        # }
        if (model == "BM") {
            res <- gls(as.formula(formula_str), correlation = corBrownian(phy = tree_sp), data = data_sp, method = "ML", na.action = na.omit)
        } else if (model == "OU") {
            # glsControl(opt = "nlminb", returnObject = TRUE, singular.ok = TRUE)  # set optimization method
            res <- gls(as.formula(formula_str), correlation = corMartins(1, phy = tree_sp, fixed = TRUE), data = data_sp, method = "ML", na.action = na.omit)
        }
        # print(sp)
        res_summary <- summary(res)
        return(data.frame(species = sp, p_value = res_summary$tTable[2, 4]))
    })
    return(c(max(p_list$p_value), p_list$species[which.max(p_list$p_value)]))
}
# function for pipeline of PGLS, including model fitting, model selection, p-value calculation, and visualization
PGLS <- function(ogid, trait, lg1 = FALSE, lg2 = TRUE) {
    ## fit PGLS model and return the maximum p-value and the species with the maximum p-value
    data_df <- ent_15[, c("species", ogid)]
    traitd <- traits[, c("species", "Order", trait)]
    traitd <- traitd[traitd$species %in% data_df$species, ]
    data_df <- merge(data_df, traitd, by = "species")

    data_df <- data_df[!is.na(data_df[[ogid]]) & !is.na(data_df[[trait]]), ]
    data_df <- data_df[data_df[[ogid]] != 0, ]
    if (lg2 == TRUE) {
        data_df[[trait]] <- log(data_df[[trait]])
    }
    spp <- data_df$species
    svg <- traits[traits$species %in% spp, c("species", "svg_path")]
    data_df <- merge(data_df, svg, by = "species")
    tree <- keep.tip(tre, spp)
    # var1 <- ifelse(lg1 == TRUE, paste("log(", ogid, ")", sep = ""), ogid)
    # var2 <- ifelse(lg2 == TRUE, paste("log(", trait, ")", sep = ""), trait)
    # formula_str <- paste(var2, " ~ ", var1, sep = "")
    formula_str <- paste(trait, " ~ ", ogid, sep = "")
    bm_model <- gls(as.formula(formula_str), correlation = corBrownian(phy = tree), data = data_df, method = "ML", na.action = na.omit)
    ou_model <- try(gls(as.formula(formula_str), correlation = corMartins(1, phy = tree, fixed = FALSE), data = data_df, method = "ML", na.action = na.omit))
    if (inherits(ou_model, "try-error")) {
        ou_model <- gls(as.formula(formula_str), correlation = corMartins(1, phy = tree, fixed = TRUE), data = data_df, method = "ML", na.action = na.omit)
    }
    models <- list(BM = bm_model, OU = ou_model)
    lrt <- anova(bm_model, ou_model)
    if (lrt$`p-value`[2] < 0.05) {
        model <- names(models)[which.max(lrt$logLik)]
    } else {
        model <- "OU"
    }
    pmax <- pgls_pmax_model(ogid, trait, data_df, tree, model, lg1 = FALSE, lg2 = FALSE)
    res_model <- models[[model]]
    res_summary <- summary(res_model)
    # res_r2 <- R2_lik(res_model, phy = tree)
    # res_r2 <- res_r2[1]
    res_r2 <- R2_lik(res_model)
    res_pred <- data.frame(x = data_df[[ogid]], pred = res_model$fitted)
    res_se <- se(data_df[[ogid]], res_model)
    res_pred <- cbind(res_pred, res_se)
    colnames(res_pred) <- c("x", "pred", "se")
    res_pred$lwr <- res_pred$pred - 1.96 * res_pred$se
    res_pred$upr <- res_pred$pred + 1.96 * res_pred$se
    res_slope <- res_summary$tTable[2, 1]
    res_intercept <- res_summary$tTable[1, 1]
    res_p <- res_summary$tTable[2, 4]
    res_plot <- ggplot(data_df, aes_string(x = ogid, y = trait)) +
        geom_point(aes(fill = Order, color = Order), size = 1, shape = 24) +
        geom_abline(intercept = res_intercept, slope = res_slope, color = "red") +
        geom_ribbon(data = res_pred, aes(x = x, y = pred, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) +
        geom_image(image = data_df$svg_path, aes(x = .data[[ogid]], y = .data[[trait]]), size = 0.05, height = 0.03, alpha = 0.5) +
        scale_fill_manual(values = order_colors) +
        scale_color_manual(values = order_colors) +
        annotate("text", label = paste(
            "R^2:", round(res_r2, 4),
            "p:", round(res_p, 4),
            "pmax:", round(as.numeric(pmax[1]), 4),
            "model:", model,
            "log trait:", lg2
        ), size = 3, hjust = 0, vjust = 1) +
        theme_classic() +
        theme(
            axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8),
            plot.title = element_text(size = 10, hjust = 0.5),
            axis.line = element_line(colour = "black"), legend.position = "none"
        ) +
        labs(title = paste("OGID:", ogid, ", Trait:", trait, ", Model:", model, "\nP-value:", round(res_p, 4), ", max P-value:", round(as.numeric(pmax[1]), 4), "\nSpecies:", pmax[2], ", R2:", round(res_r2, 4), ", log trait:", lg2))
    ggsave(paste0("/media/Data/zhangz/chip/analysis/summary2/pgls/res/", trait, "/", trait, "_", ogid, "_", ".pdf"), res_plot, width = 10, height = 10, units = "cm")
    res <- data.frame(names = c("OGID", "Trait", "Model", "P-value", "Max_P-value", "Species", "R2", "Slope", "log_trait"), values = c(ogid, trait, model, res_p, as.numeric(pmax[1]), pmax[2], res_r2, res_slope, lg2))
    write.csv(res, file = paste0("/media/Data/zhangz/chip/analysis/summary2/pgls/res/", trait, "/", trait, "_", ogid, "_", ".csv"), row.names = FALSE)
    res <- data.frame(t(res))
    colnames(res) <- res[1, ]
    res <- res[-1, ]
    return(res)
}

## functions
find_mus_ortho <- function(og) {
    # print(og)
    mus_gene <- ortho_all[ortho_all$Orthogroup == og, "Mus_musculus"]
    mus_gene <- str_split(mus_gene, ", ")[[1]]
    return(paste(str_split(mus_gene[1], "[|]")[[1]], collapse = ";"))
    # return(str_split(mus_gene[1], "[|]")[[1]][1])
}

## load data and preprocess
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

order_colors <- c("Artiodactyla" = "#3682be", "Carnivora" = "#45a776", "Perissodactyla" = "#f05330", "Chiroptera" = "#eed777", "Eulipotyphla" = "#38cb7d", "Rodentia" = "#334f65", "Lagomorpha" = "#ddae33", "Primates" = "#b3974e", "Scandentia" = "#844bb3", "Hyracoidea" = "#93c555", "Diprotodontia" = "#5f6694")

# ent1_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ent_ouwie/longevity/ent1_all.csv", header = TRUE)
# ent1_all[is.na(ent1_all)] <- 0
# # calculate how many species have data for each OGID, how many cols in each row are not 0
# ent1_all$non0 <- apply(ent1_all[, 2:ncol(ent1_all)], 1, function(x) sum(x != 0))
# nrow(ent1_all[ent1_all$non0 >= 15, ])
# 667
# ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)

ent_all <- adply(species, 1, function(sp) {
    # read ent1 from file and concate together, no matter gene copy number
    ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
    res <- data.frame(t(ent[, c("OGID", "ent1")]))
    colnames(res) <- res[1, ]
    res <- res[-1, ]
    rownames(res) <- sp
    new_col <- data.frame(species = sp)
    res <- cbind(new_col, res)
    return(res)
}, .parallel = TRUE)
ent_all[is.na(ent_all)] <- 0
ent_all <- ent_all[, -1]
ent_all[, 2:ncol(ent_all)] <- apply(ent_all[, 2:ncol(ent_all)], 2, function(x) as.numeric(x))
head(ent_all[, 1:5])
ent_non0 <- data.frame(OGID = colnames(ent_all)[2:ncol(ent_all)], non0 = apply(ent_all[, 2:ncol(ent_all)], 2, function(x) sum(x != 0)))
nrow(ent_non0[ent_non0$non0 >= 15, ])
# 1425
ortho_15 <- ent_non0$OGID[ent_non0$non0 >= 15]
ent_15 <- ent_all[, c("species", ent_non0$OGID[ent_non0$non0 >= 15])]

ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)
og_mark <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/og_mark.csv", header = TRUE)


# read tree
tre <- read.newick("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
tre <- keep.tip(tre, species)

# read life history traits
traits <- read.csv("/media/Data/zhangz/chip/analysis/summary2/tissue_ent/traits.csv", header = TRUE)
traits <- traits[traits$species %in% species, c(1, 2, 12:42)]
traits$svg_path <- paste0("/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/", traits$species, ".svg")
traits$svg_path[which(traits$species %in% c("Myotis_chinensis", "Myotis_ricketti", "Rhinolophus_ferrumequinum", "Rhinolophus_pusillus"))] <- "/media/Data/zhangz/chip/analysis/summary2/fig/silhouette/Hipposideros_larvatus.svg"
for (svg in traits$svg_path) {
    if (!file.exists(svg)) {
        print(svg)
    }
}
lg_traits <- c("body_colume", "FTM.d.y", "AW.g", "Gestation.days.", "BMR.O2ml.h.", "BMR_mass.g.", "Weaning.days.", "Birth_weight.g.")

# conduct PGLS
pgls_res <- mdply(expand.grid(ogid = ortho_15, trait = colnames(traits)[3:(ncol(traits) - 1)]), function(ogid, trait) {
    # create trait dir if not exists
    # print(paste(ogid, trait))
    ogid <- as.character(ogid)
    trait <- as.character(trait)
    if (!dir.exists(paste0("/media/Data/zhangz/chip/analysis/summary2/pgls/res/", trait))) {
        dir.create(paste0("/media/Data/zhangz/chip/analysis/summary2/pgls/res/", trait))
    }
    if (trait %in% c("Echolocation", "Hibernation", "Habit", "Diet", "Sociality", "Behavior", "Feed", "svg_path")) {
        return(NULL)
    }
    if (trait %in% lg_traits) {
        res <- try(PGLS(ogid, trait, lg1 = FALSE, lg2 = TRUE))
    } else {
        res <- try(PGLS(ogid, trait, lg1 = FALSE, lg2 = FALSE))
    }
    if (inherits(res, "try-error")) {
        print(paste(trait, ogid))
        return(data.frame(x = ogid, error = as.character(res)))
    }
    return(res)
}, .parallel = TRUE)
write.csv(pgls_res, file = "/media/Data/zhangz/chip/analysis/summary2/pgls/pgls_res.csv", row.names = FALSE)

pgls_res <- read.csv("/media/Data/zhangz/chip/analysis/summary2/pgls/pgls_res.csv", header = TRUE)
pgls_res[["P-value"]] <- as.numeric(pgls_res[["P-value"]])
pgls_res[["Max_P-value"]] <- as.numeric(pgls_res[["Max_P-value"]])
pgls_res[["R2"]] <- as.numeric(pgls_res[["R2"]])
head(pgls_res)
summary(pgls_res)
head(pgls_res[order(pgls_res[["Max_P-value"]]), ])
nrow(pgls_res[pgls_res[["Max_P-value"]] < 0.05, ])
sig_pgls_res <- pgls_res[pgls_res[["Max_P.value"]] < 0.05, ]
summary(sig_pgls_res)
tapply(sig_pgls_res$R2, sig_pgls_res$trait, summary)
tapply(sig_pgls_res$R2, sig_pgls_res$trait, length)
sig_pgls_res <- sig_pgls_res[!is.na(sig_pgls_res$OGID), ]
pgls_res[!is.na(pgls_res$error), ]

## there are four ogid fail to fit PGLS model in Temperature.K.
## OG0009116 OG0010155 OG0011506 OG0003532
## only OG0010155 has a significant p-value in Temperature.K., with a R2 = 0.4054453
## the other three have a p-value > 0.05
## so we just ignore them

## plot the significant PGLS results
sig_pgls_res <- sig_pgls_res %>% arrange(trait, `P-value`)
sig_pgls_res <- sig_pgls_res %>% arrange(trait, P.value)
tapply(sig_pgls_res$trait, sig_pgls_res$Model, length)
#  BM  OU
#   8 329
sig_pgls_res <- sig_pgls_res[, 3:11]
sig_pgls_res$gene <- sapply(sig_pgls_res$OGID, find_mus_ortho)
sig_pgls_res <- sig_pgls_res[!is.na(sig_pgls_res$Trait), ]
sig_pgls_res[sig_pgls_res$Trait == "MLres.y", ]
sig_pgls_res[sig_pgls_res$Trait %in% c("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y"), ]
sig_pgls_res[which(sig_pgls_res$Trait %in% c("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y") & sig_pgls_res$Slope > 0), ]
intersect(sig_pgls_res[sig_pgls_res$Trait == "ML.yrs.y", "OGID"], sig_pgls_res[sig_pgls_res$Trait == "MLres.y", "OGID"])
upset_matrix <- table(sig_pgls_res$OGID, sig_pgls_res$Trait)
upset_df <- as.data.frame(upset_matrix)
colnames(upset_df) <- c("OGID", "Trait", "Frequency")

# long df upset_df to wide df
wide_df <- spread(upset_df, Trait, Frequency)
# pdf("/media/Data/zhangz/chip/analysis/summary2/pgls/upset.pdf", width = 10, height = 10)
require(UpSetR)
p1 <- upset(wide_df,
    sets = c("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y", "AW.g"),
    nset = 23, nintersects = 60,
    queries = list(
        list(
            query = intersects, params = list("ML.yrs.y", "MLres.y"),
            color = "red", active = TRUE
        ),
        list(
            query = intersects, params = list("MLres.y", "FTMres.y"),
            color = "blue", active = TRUE
        ),
        list(
            query = intersects, params = list("ML.yrs.y", "FTM.d.y"),
            active = TRUE
        ),
        list(
            query = intersects, params = list("FTMres.y", "FTM.d.y"),
            active = TRUE
        )
    ),
)
# require(ggplotify)
# g1 <- as.ggplot(p1)
# ggsave("/media/Data/zhangz/chip/analysis/summary2/pgls/upset.pdf", g1, width = 10, height = 10, units = "cm")
pdf("/media/Data/zhangz/chip/analysis/summary2/pgls/upset_longevity.pdf", width = 10, height = 10)
print(p1)
dev.off()
p2 <- upset(wide_df,
    nset = 23, nintersects = 60,
    queries = list(
        # list(
        #     query = intersects, params = list("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y"),
        #     color = "red", active = TRUE
        # ),
        list(
            query = intersects, params = list("ML.yrs.y", "MLres.y"),
            active = TRUE
        ),
        list(
            query = intersects, params = list("MLres.y", "FTMres.y"),
            active = TRUE
        )
    ),
)
pdf("/media/Data/zhangz/chip/analysis/summary2/pgls/upset_all.pdf", width = 16, height = 9)
print(p2)
dev.off()

write.csv(sig_pgls_res, file = "/media/Data/zhangz/chip/analysis/summary2/pgls/sig_pgls_res.csv", row.names = FALSE)

## gene ontology analysis of high entropy genes in Mus_musculus
library(clusterProfiler)
library(org.Mm.eg.db)
detach("package:KEGG.db", unload = TRUE)
library(KEGG.db)

## prepare kegg msa db
remotes::install_github("YuLab-SMU/createKEGGdb")
library(createKEGGdb)
setwd("/media/Data/zhangz/chip/analysis/summary2/ent_gene")
createKEGGdb::create_kegg_db("mmu")
install.packages("./KEGG.db_1.0.tar.gz", repos = NULL, type = "source")

## load data
sp <- "Mus_musculus"
mus_ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
mus_ent <- mus_ent[mus_ent$ent1 > 0, ]
summary(mus_ent)
plot(density(mus_ent$ent1))
mus_ent <- mus_ent[order(mus_ent$unique_ele_num, decreasing = TRUE), ]

## get the top 100 high entropy genes
top100 <- mus_ent$gene[1:100]

## do GO and KEGG analysis
top100 <- unique(top100)
top100_GO <- enrichGO(top100, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
top100_KEGG <- enrichKEGG(top100,
    organism = "mmu",
    keyType = "SYMBOL", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
    use_internal_data = TRUE
)
# plot top 10 GO terms
go_dot <- dotplot(top100_GO, showCategory = 15)
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/go_top100_dot.pdf", go_dot, width = 8, height = 9)
go_cnet_p <- cnetplot(top100_GO, categorySize = "pvalue", foldChange = "qvalue", showCategory = 15)
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/go_top100_cnet_p.pdf", go_cnet_p, width = 12, height = 9)
head(og_mark)
og_mark$gene <- sapply(og_mark$Orthogroup, find_mus_ortho)
og_mark[og_mark$Orthogroup %in% c("OG0009639", "OG0006213", "OG0007610"), ]
top100_gene <- og_mark[og_mark$Orthogroup %in% mus_ent$OGID[1:100], ]
top100_gene$gene <- sapply(top100_gene$Orthogroup, find_mus_ortho)

## find a well studied example enhancer of mouse genome
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/fpkm_ortho.RData")
mus_anno <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE)
head(mus_anno[which(mus_anno$chr == "NC_000078.7" & mus_anno$start > 112500000 & mus_anno$start < 112600000), c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "act_time")])
head(mus_anno[grepl("Nenf", mus_anno$overlap_gene), c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "act_time")])
mus_anno[582, ]
head(mus_anno[grepl("Shh", mus_anno$overlap_gene), c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "act_time")])
head(mus_anno[grepl("Abcf3", mus_anno$overlap_gene), c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "act_time")])
head(mus_anno[grepl("Inf2", mus_anno$overlap_gene), c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "note_abc", "overlap_distance", "act_time", "s0_0", "s0_1", "s1_0", "s1_1")], 20)
head(mus_anno[grepl("Lactb2", mus_anno$overlap_gene), c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "note_abc", "overlap_distance", "act_time", "s0_0", "s0_1", "s1_0", "s1_1")], 20)
grep("abc", mus_anno[grepl("Lactb2", mus_anno$overlap_gene), "note_abc"])
mus_anno[grepl("Lactb2", mus_anno$overlap_gene), c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "note_abc", "overlap_distance", "act_time", "s0_0", "s0_1", "s1_0", "s1_1")][c(3, 10, 18, 23, 33, 34), ]
mus_anno[879, ]
mus_anno[30657, ]

# OG0014307 = Shh
## There are some enhancers around Shh gene, but none of them receive support from abc model, and average expression is not higher than 1
## OG0007610 = Abcf3
## OG0009639 = Lactb2
# OG0008870 = Inf2
mus_ent[mus_ent$OGID == "OG0008870", ]
#          OGID gene     fpkm  fpkm_sd  fpkm_cv adj_fpkm adj_fpkm_sd adj_fpkm_cv qn_fpkm_sps qn_fpkm_sd qn_fpkm_cv      tau      ent1 ent2 unique_ele_num
# 277 OG0008870 Inf2 15.52847 19.69557 1.268353 1.960296   0.8643462    44.09264    17.65747    23.1586   1.311547 0.886453 0.1798256    1              1
# maybe Inf2 is good choice

library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
mus_anno$specificity <- "N"
mus_anno$specificity[which(mus_anno$act_time > 0 &
    (mus_anno$s1_1 > 0 | mus_anno$s0_1 > 0 | !mus_anno$ele_pattern %in% c("B", "K", "L")) &
    mus_anno$element == "enhancer")] <- "evol_pleio_e"
mus_anno$specificity[which(mus_anno$act_time > 0 &
    (mus_anno$s1_1 > 0 | mus_anno$s0_1 > 0 | !mus_anno$ele_pattern %in% c("B", "K", "L")) &
    mus_anno$element == "promoter")] <- "evol_pleio_p"
mus_anno$specificity[which(mus_anno$act_time > 0 &
    (mus_anno$s1_1 == 0 & mus_anno$s0_1 == 0 & mus_anno$ele_pattern %in% c("B", "K", "L")) &
    mus_anno$element == "enhancer")] <- "evol_ts_e"
mus_anno$specificity[which(mus_anno$act_time > 0 &
    (mus_anno$s1_1 == 0 & mus_anno$s0_1 == 0 & mus_anno$ele_pattern %in% c("B", "K", "L")) &
    mus_anno$element == "promoter")] <- "evol_ts_p"
mus_anno$specificity[which(mus_anno$act_time == 0 &
    mus_anno$ele_pattern %in% c("B", "K", "L") &
    mus_anno$element == "enhancer")] <- "sp_ts_e"
mus_anno$specificity[which(mus_anno$act_time == 0 &
    mus_anno$ele_pattern %in% c("B", "K", "L") &
    mus_anno$element == "promoter")] <- "sp_ts_p"
mus_anno$specificity[which(mus_anno$act_time == 0 &
    !mus_anno$ele_pattern %in% c("B", "K", "L") &
    mus_anno$element == "enhancer")] <- "sp_pleio_e"
mus_anno$specificity[which(mus_anno$act_time == 0 &
    !mus_anno$ele_pattern %in% c("B", "K", "L") &
    mus_anno$element == "promoter")] <- "sp_pleio_p"

spec_color <- c(
    "evol_pleio_e" = "#b23f3f", "evol_ts_e" = "#ce6868", "sp_ts_e" = "#ffbebe", "sp_pleio_e" = "#ffd7d7",
    "evol_pleio_p" = "#3a6c82", "evol_ts_p" = "#548195", "sp_ts_p" = "#99b7c5", "sp_pleio_p" = "#c0d3dc"
)

# select CREs in NC_000078.7，`mus_anno$start > 112500000 & mus_anno$start < 112600000`
cre_inf <- mus_anno[
    which(mus_anno$chr == "NC_000078.7" & mus_anno$start > 112500000 & mus_anno$start < 112700000),
    c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "act_time", "specificity")
]
cre_inf$strand <- "."
cre_inf <- cre_inf[order(cre_inf$start), c("chr", "start", "end", "strand", "overlap_og", "overlap_gene", "tissue", "element", "align", "act_time", "specificity")]
tapply(cre_inf$specificity, cre_inf$element, table)
# $enhancer
# evol_pleio_e    evol_ts_e   sp_pleio_e      sp_ts_e
#           11            2            1           11

# $promoter
# evol_pleio_p    evol_ts_p      sp_ts_p
#            2            1            1
unique(cre_inf$tissue)
# mark cre receive support from abc model
cre_inf$abc <- FALSE
cre_inf$abc[which(cre_inf$overlap_og == "Inf2")] <- TRUE
colnames(cre_inf) <- c("seqnames", "start", "end", "strand", "OGID", "gene", "tissue", "element_type", "align", "act_time", "specificity", "abc")
cre_inf$fill <- spec_color[match(cre_inf$specificity, names(spec_color))]
# Inf2 is on chr12
cre_inf$seqnames <- "chr12"
cre_inf_brain <- cre_inf[which(cre_inf$tissue == "Brain"), ]
cre_inf_liver <- cre_inf[which(cre_inf$tissue == "Liver"), ]
cre_inf_kidney <- cre_inf[which(cre_inf$tissue == "Kidney"), ]
cre_gr_brain <- makeGRangesFromDataFrame(cre_inf_brain, keep.extra.columns = TRUE)
cre_gr_liver <- makeGRangesFromDataFrame(cre_inf_liver, keep.extra.columns = TRUE)
cre_gr_kidney <- makeGRangesFromDataFrame(cre_inf_kidney, keep.extra.columns = TRUE)
cre_gr <- makeGRangesFromDataFrame(cre_inf, keep.extra.columns = TRUE)


# plot using Gviz
cre_atr_brain <- AnnotationTrack(cre_gr_brain, name = "Brain CRE", feature = as.vector(cre_inf_brain$fill), stacking = "dense")
cre_atr_liver <- AnnotationTrack(cre_gr_liver, name = "Liver CRE", feature = as.vector(cre_inf_liver$fill), stacking = "dense")
cre_atr_kidney <- AnnotationTrack(cre_gr_kidney, name = "Kidney CRE", feature = as.vector(cre_inf_kidney$fill), stacking = "dense")
gr_plot <- plotTracks(list(cre_atr_brain, cre_atr_liver, cre_atr_kidney))
genome(cre_atr_brain) <- "mm39"
gen <- "mm39"
chr <- "chr12"
itrack <- IdeogramTrack(genome = gen, chromosome = chr, centromereShape = "Circle")
gtrack <- GenomeAxisTrack()
cre_gr$group <- cre_inf$tissue
cre_track <- AnnotationTrack(cre_gr,
    name = "CRE", feature = as.vector(cre_inf$fill),
    showId = TRUE, stacking = "squish"
)
genome(cre_track) <- "mm39"
displayPars(cre_atr_brain) <- list(fill = cre_inf_brain$fill)
displayPars(cre_atr_liver) <- list(fill = cre_inf_liver$fill)
displayPars(cre_atr_kidney) <- list(fill = cre_inf_kidney$fill)
displayPars(cre_track) <- list(fill = cre_inf$fill)
plotTracks(list(itrack, gtrack, cre_track), from = 112400000, to = 112700000)


# add gene track
library(org.Mm.eg.db)
gene_track <- GeneRegionTrack(
    "/media/Data/zhangz/chip/genomes/Mus_musculus/GCF_000001635.27_GRCm39_genomic_chr12.gff3",
    genome = gen, chromosome = chr, name = "Gene",
    transcriptAnnotation = "symbol", showId = TRUE,
    showFeatureId = TRUE, showFeatureType = TRUE,
    groupAnnotation = "gene"
)
# Error : subscript contains out-of-bounds indices
# Warning message:
# In .import.gff3(file) :
#   File '/media/Data/zhangz/chip/genomes/Mus_musculus/GCF_000001635.27_GRCm39_genomic_chr12.gff3' is not valid according to the GFF3 standard and can not be properly parsed.
# Results may not be what you expected!

plotTracks(gene_track,
    from = 112400000, to = 112700000, stacking = "full",
    transcriptAnnotation = "symbol", showId = TRUE, groupAnnotation = "gene",
    just.group = "above"
)
# whole piece, no gene
# biom_track <- BiomartGeneRegionTrack(
#     genome = gen, chromosome = chr, start = 112400000, end = 112700000,
#     name = "Gene",
# )
# txdb <- loadDb("/media/Data/zhangz/chip/genomes/Mus_musculus/Mumu_AH84139_refGene.sqlite")
gene2_track <- GeneRegionTrack(
    txdb,
    genome = gen, chromosome = chr, name = "Gene",
    transcriptAnnotation = "SYMBOL"
)
plotTracks(gene2_track, from = 112400000, to = 112700000, stacking = "squish", transcriptAnnotation = "SYMBOL")
# gene names are NM_177039 style, cannot read gene names
knownGenes <- UcscTrack(
    genome = "mm39", chromosome = "chr12",
    track = "knownGene", table = "knownGene",
    from = 112400000, to = 112700000,
    trackType = "GeneRegionTrack",
    rstarts = "exonStarts", rends = "exonEnds",
    gene = "name", symbol = "name",
    transcript = "name", strand = "strand",
    fill = "#8282d2", name = "UCSC Genes"
)
plotTracks(refGenes)

# read gff file as GRanges
library(rtracklayer)
mus_gtf <- import("/media/Data/zhangz/chip/genomes/Mus_musculus/GCF_000001635.27_GRCm39_genomic.gff3", format = "GFF")
# change the seqnames from NCBI style to chr style
seqlevels(mus_gtf)
chr_names <- read.csv("/media/Data/zhangz/chip/genomes/Mus_musculus/mm39.chromAlias.txt", header = TRUE, sep = "\t")
colnames(chr_names)[1] <- c("ucsc")
chr_map <- setNames(chr_names$ucsc, chr_names$refseq)
seqlevels(mus_gtf) <- chr_map[seqlevels(mus_gtf)]
mus_chr12 <- mus_gtf[seqnames(mus_gtf) == "chr12"]
# selecet genes and exons, utr, cds
mus_chr12 <- mus_chr12[which(mus_chr12$type %in% c("gene")), ]
mus_chr12 <- mus_chr12[!grepl("Gm", mus_chr12$gene), ]
mus_chr12 <- mus_chr12[!grepl("Mir", mus_chr12$gene), ]
chr12_track <- GeneRegionTrack(mus_chr12,
    genome = gen, chromosome = chr,
    name = "Gene", transcriptAnnotation = "gene",
    showId = TRUE, showFeatureId = TRUE, showFeatureType = TRUE,
    groupAnnotation = "gene", shape = "arrow", stackHeight = 0.3
)
plotTracks(chr12_track,
    chr = "chr12", from = 112400000, to = 112700000,
    stacking = "squish", transcriptAnnotation = "gene"
)

# tis_p <- plotTracks(list(itrack, gtrack, gene_track, cre_atr_brain, cre_atr_liver, cre_atr_kidney),
#     from = 112400000, to = 112700000, sizes = c(3, 3, 5, 10, 10, 10)
# )
group(gene_track)
head(symbol(gene_track))

ht <- HighlightTrack(
    trackList = list(chr12_track, cre_atr_brain, cre_atr_liver, cre_atr_kidney),
    start = 112554000, width = 1000, fill = "white", alpha = 0.3,
    chromosome = "chr12", sizes = c(1, 1, 1, 1)
)
tis_p <- plotTracks(list(itrack, gtrack, ht),
    from = 112500000, to = 112700000
)
tis_gg <- as.ggplot(grid2grob(plotTracks(list(itrack, gtrack, ht),
    from = 112500000, to = 112700000
)))
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/Inf2_cre.pdf", tis_gg, width = 8, height = 6)


## read bigwig file and plot histone peak
## mouse brain H3K27ac bigwig file: /DELL_EMC/zhangz/chip/analysis/Mus_musculus/chip2/Brain/H3K27ac/chip/3c581f2b-85df-4df3-a55c-11786f452e68/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## mouse brain H3K4me3 bigwig file: /DELL_EMC/zhangz/chip/analysis/Mus_musculus/chip2/Brain/H3K4me3/chip/21d7a3e9-3267-4e0d-b9ba-d3446acdc52c/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## rat brain H3K27ac bigwig file: /DELL_EMC/zhangz/chip/analysis/Rattus_norvegicus/chip2/Brain/H3K27ac/chip/fd70c7c8-f3b5-45fe-944f-a0287fdfbab9/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## rat brain H3K4me3 bigwig file: /DELL_EMC/zhangz/chip/analysis/Rattus_norvegicus/chip2/Brain/H3K4me3/chip/67bd2ad7-3e63-45e0-a072-8816d81d6257/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## rat liver H3K27ac bigwig file: /DELL_EMC/zhangz/chip/analysis/Rattus_norvegicus/chip2/Liver/H3K27ac/chip/ff442024-254a-4382-a1b2-c3e54c7a3874/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## rat liver H3K4me3 bigwig file: /DELL_EMC/zhangz/chip/analysis/Rattus_norvegicus/chip2/Liver/H3K4me3/chip/0e6825f7-de68-423c-9ff1-a2237155e12e/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Tupaia_belangeri brain H3K27ac bigwig file: /DELL_EMC/zhangz/chip/analysis/Tupaia_belangeri/chip2/Brain/H3K27ac/chip/f7ab1ff2-5b9c-4852-918a-7ff72dfa9993/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Tupaia_belangeri brain H3K4me3 bigwig file: /DELL_EMC/zhangz/chip/analysis/Tupaia_belangeri/chip2/Brain/H3K4me3/chip/3f7b06ca-bbae-4538-a8f3-9d66be0d16bf/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Tupaia_belangeri kidney H3K27ac bigwig file: /DELL_EMC/zhangz/chip/analysis/Tupaia_belangeri/chip2/Kidney/H3K27ac/chip/d6b600a5-3f3c-4d77-b455-ad5523d1d72a/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Tupaia_belangeri kidney H3K4me3 bigwig file: /DELL_EMC/zhangz/chip/analysis/Tupaia_belangeri/chip2/Kidney/H3K4me3/chip/82d8f37f-e4d0-4f62-8468-f12c79263c2a/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
mouse_brain_H3K27ac_bw <- "/DELL_EMC/zhangz/chip/analysis/Mus_musculus/chip2/Brain/H3K27ac/chip/3c581f2b-85df-4df3-a55c-11786f452e68/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
mouse_brain_H3K4me3_bw <- "/DELL_EMC/zhangz/chip/analysis/Mus_musculus/chip2/Brain/H3K4me3/chip/21d7a3e9-3267-4e0d-b9ba-d3446acdc52c/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
rat_brain_H3K27ac_bw <- "/DELL_EMC/zhangz/chip/analysis/Rattus_norvegicus/chip2/Brain/H3K27ac/chip/fd70c7c8-f3b5-45fe-944f-a0287fdfbab9/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
rat_brain_H3K4me3_bw <- "/DELL_EMC/zhangz/chip/analysis/Rattus_norvegicus/chip2/Brain/H3K4me3/chip/67bd2ad7-3e63-45e0-a072-8816d81d6257/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
rat_liver_H3K27ac_bw <- "/DELL_EMC/zhangz/chip/analysis/Rattus_norvegicus/chip2/Liver/H3K27ac/chip/ff442024-254a-4382-a1b2-c3e54c7a3874/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
rat_liver_H3K4me3_bw <- "/DELL_EMC/zhangz/chip/analysis/Rattus_norvegicus/chip2/Liver/H3K4me3/chip/0e6825f7-de68-423c-9ff1-a2237155e12e/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
tupaia_brain_H3K27ac_bw <- "/DELL_EMC/zhangz/chip/analysis/Tupaia_belangeri/chip2/Brain/H3K27ac/chip/f7ab1ff2-5b9c-4852-918a-7ff72dfa9993/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
tupaia_brain_H3K4me3_bw <- "/DELL_EMC/zhangz/chip/analysis/Tupaia_belangeri/chip2/Brain/H3K4me3/chip/3f7b06ca-bbae-4538-a8f3-9d66be0d16bf/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
tupaia_kidney_H3K27ac_bw <- "/DELL_EMC/zhangz/chip/analysis/Tupaia_belangeri/chip2/Kidney/H3K27ac/chip/d6b600a5-3f3c-4d77-b455-ad5523d1d72a/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
tupaia_kidney_H3K4me3_bw <- "/DELL_EMC/zhangz/chip/analysis/Tupaia_belangeri/chip2/Kidney/H3K4me3/chip/82d8f37f-e4d0-4f62-8468-f12c79263c2a/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"

options(ucscChromosomeNames = FALSE)
mouse_brain_H3K27ac_dr <- DataTrack(
    range = mouse_brain_H3K27ac_bw, genome = "mm39", chr = "NC_000078.7",
    type = "mountain", name = "Mouse Brain H3K27ac", alpha = 0.5
)
mouse_brain_H3K4me3_dr <- DataTrack(
    range = mouse_brain_H3K4me3_bw, genome = "mm39", chr = "NC_000078.7",
    type = "mountain", name = "Mouse Brain H3K4me3", alpha = 0.5
)
# NC_000078.7 112554156 112555002
plotTracks(c(mouse_brain_H3K27ac_dr, mouse_brain_H3K4me3_dr), from = 112554000, to = 112555000, ylim = c(2, 6))

ot <- OverlayTrack(trackList = list(mouse_brain_H3K27ac_dr, mouse_brain_H3K4me3_dr))
plotTracks(ot, from = 112554000, to = 112555000, ylim = c(2, 6))


# another gene: OG0009639 = Lactb2, OG0013300 = Xkr9
#      OGID               species    Brain   Kidney    Liver
# OG0009639          Mus_musculus   5.3644 349.5253 222.6958
# OG0013300          Mus_musculus   0.0000   0.1924  10.7914
mus_anno[30657, ]
mus_anno[30655:30660, ]
# enh_k_130303_K_Peak_15875 Peak_15875 Kidney enhancer NC_000067.7=chr1 13740465 13741263

# specificity and fill color
mus_anno$specificity <- "N"
mus_anno$specificity[which(mus_anno$act_time > 0 &
    (mus_anno$s1_1 > 0 | mus_anno$s0_1 > 0 | !mus_anno$ele_pattern %in% c("B", "K", "L")) &
    mus_anno$element == "enhancer")] <- "evol_pleio_e"
mus_anno$specificity[which(mus_anno$act_time > 0 &
    (mus_anno$s1_1 > 0 | mus_anno$s0_1 > 0 | !mus_anno$ele_pattern %in% c("B", "K", "L")) &
    mus_anno$element == "promoter")] <- "evol_pleio_p"
mus_anno$specificity[which(mus_anno$act_time > 0 &
    (mus_anno$s1_1 == 0 & mus_anno$s0_1 == 0 & mus_anno$ele_pattern %in% c("B", "K", "L")) &
    mus_anno$element == "enhancer")] <- "evol_ts_e"
mus_anno$specificity[which(mus_anno$act_time > 0 &
    (mus_anno$s1_1 == 0 & mus_anno$s0_1 == 0 & mus_anno$ele_pattern %in% c("B", "K", "L")) &
    mus_anno$element == "promoter")] <- "evol_ts_p"
mus_anno$specificity[which(mus_anno$act_time == 0 &
    mus_anno$ele_pattern %in% c("B", "K", "L") &
    mus_anno$element == "enhancer")] <- "sp_ts_e"
mus_anno$specificity[which(mus_anno$act_time == 0 &
    mus_anno$ele_pattern %in% c("B", "K", "L") &
    mus_anno$element == "promoter")] <- "sp_ts_p"
mus_anno$specificity[which(mus_anno$act_time == 0 &
    !mus_anno$ele_pattern %in% c("B", "K", "L") &
    mus_anno$element == "enhancer")] <- "sp_pleio_e"
mus_anno$specificity[which(mus_anno$act_time == 0 &
    !mus_anno$ele_pattern %in% c("B", "K", "L") &
    mus_anno$element == "promoter")] <- "sp_pleio_p"

spec_color <- c(
    "evol_pleio_e" = "#b23f3f", "evol_ts_e" = "#ce6868", "sp_ts_e" = "#ffbebe", "sp_pleio_e" = "#ffd7d7",
    "evol_pleio_p" = "#3a6c82", "evol_ts_p" = "#548195", "sp_ts_p" = "#99b7c5", "sp_pleio_p" = "#c0d3dc"
)

# select region around Lactb2
cre_inf2 <- mus_anno[
    which(mus_anno$chr == "NC_000067.7" & mus_anno$start > 13690000 & mus_anno$start < 13750000),
    c("chr", "start", "end", "overlap_og", "overlap_gene", "tissue", "element", "align", "act_time", "specificity")
]
# 这里的起止位置必须和最终显示是一样的，因为Gviz查找fill参数时会按照最终绘图的from决定的第一个位置和fill的第一个参数匹配
# 而不是在指定fill参数后就固定键值对，所以每次更改绘图范围都需要重新设置fill参数， 否则颜色会错乱
cre_inf2$strand <- "."
cre_inf2 <- cre_inf2[order(cre_inf2$start), c("chr", "start", "end", "strand", "overlap_og", "overlap_gene", "tissue", "element", "align", "act_time", "specificity")]
tapply(cre_inf2$specificity, cre_inf2$element, table)
# $enhancer

# evol_pleio_e    evol_ts_e   sp_pleio_e      sp_ts_e
#           26            2            1           17

# $promoter

# evol_pleio_p    evol_ts_p   sp_pleio_p      sp_ts_p
#            5            2            4            3

unique(cre_inf$tissue)
# mark cre receive support from abc model
cre_inf2$abc <- FALSE
cre_inf2$abc[which(cre_inf2$overlap_og == "Lactb2")] <- TRUE
colnames(cre_inf2) <- c("seqnames", "start", "end", "strand", "OGID", "gene", "tissue", "element_type", "align", "act_time", "specificity", "abc")
cre_inf2$fill <- spec_color[match(cre_inf2$specificity, names(spec_color))]
# Inf2 is on chr12
cre_inf2$seqnames <- "chr1"
cre_inf_brain2 <- cre_inf2[which(cre_inf2$tissue == "Brain"), ]
cre_inf_liver2 <- cre_inf2[which(cre_inf2$tissue == "Liver"), ]
cre_inf_kidney2 <- cre_inf2[which(cre_inf2$tissue == "Kidney"), ]
cre_gr_brain2 <- makeGRangesFromDataFrame(cre_inf_brain2, keep.extra.columns = TRUE)
cre_gr_liver2 <- makeGRangesFromDataFrame(cre_inf_liver2, keep.extra.columns = TRUE)
cre_gr_kidney2 <- makeGRangesFromDataFrame(cre_inf_kidney2, keep.extra.columns = TRUE)
cre_gr2 <- makeGRangesFromDataFrame(cre_inf2, keep.extra.columns = TRUE)

# plot using Gviz
cre_atr_brain2 <- AnnotationTrack(cre_gr_brain2, name = "Brain CRE", feature = as.vector(cre_gr_brain2$fill), stacking = "squish")
cre_atr_liver2 <- AnnotationTrack(cre_gr_liver2, name = "Liver CRE", feature = as.vector(cre_gr_liver2$fill), stacking = "squish")
cre_atr_kidney2 <- AnnotationTrack(cre_gr_kidney2, name = "Kidney CRE", feature = as.vector(cre_gr_kidney2$fill), stacking = "squish")
# cre_atr <- AnnotationTrack(cre_gr2, name = "CRE", feature = as.vector(cre_inf2$fill), stacking = "squish")
# gr_plot2 <- plotTracks(list(cre_atr_brain2, cre_atr_liver2, cre_atr_kidney2))
genome(cre_atr_brain2) <- "mm39"
gen <- "mm39"
chr2 <- "chr1"
itrack2 <- IdeogramTrack(genome = gen, chromosome = chr2, centromereShape = "Circle")
gtrack <- GenomeAxisTrack()
cre_gr2$group <- cre_inf2$tissue
cre_track2 <- AnnotationTrack(cre_gr2,
    name = "CRE", feature = as.vector(cre_gr2$fill),
    showId = TRUE, stacking = "squish"
)
# 这里需要用cre_gr2$fill，而不是cre_inf2$fill, 否则着色不对
genome(cre_track2) <- "mm39"
displayPars(cre_atr_brain2) <- list(fill = cre_gr_brain2$fill)
displayPars(cre_atr_liver2) <- list(fill = cre_gr_liver2$fill)
displayPars(cre_atr_kidney2) <- list(fill = cre_gr_kidney2$fill)
displayPars(cre_track2) <- list(fill = cre_gr2$fill)
plotTracks(list(itrack2, gtrack, cre_track2), chr = chr2, from = 13740000, to = 13750000)
plotTracks(list(itrack2, gtrack, cre_atr_brain2, cre_atr_kidney2, cre_atr_liver2), chr = chr2, from = 13000000, to = 14000000)

mus_chr1 <- mus_gtf[seqnames(mus_gtf) == "chr1"]
# selecet genes and exons, utr, cds
mus_chr1 <- mus_chr1[which(mus_chr1$type %in% c("gene")), ]
mus_chr1 <- mus_chr1[!grepl("Gm", mus_chr1$gene), ]
mus_chr1 <- mus_chr1[!grepl("Mir", mus_chr1$gene), ]
chr1_track <- GeneRegionTrack(mus_chr1,
    genome = gen, chromosome = chr2,
    name = "Gene", transcriptAnnotation = "gene",
    showId = TRUE, showFeatureId = TRUE, showFeatureType = TRUE,
    groupAnnotation = "gene", shape = "arrow", stackHeight = 0.3
)
plotTracks(chr1_track,
    chr = "chr1", from = 13600000, to = 13900000,
    stacking = "squish", transcriptAnnotation = "gene"
)

ht2 <- HighlightTrack(
    trackList = list(chr1_track, cre_atr_brain2, cre_atr_kidney2, cre_atr_liver2),
    start = 13740465, width = 800, fill = "white", alpha = 0.3,
    chromosome = "chr1", sizes = c(1, 1, 1, 1), inBackground = TRUE
)
tis_p <- plotTracks(list(itrack2, gtrack, ht2),
    from = 13690000, to = 13750000
)
tis_gg2 <- as.ggplot(grid2grob(plotTracks(list(itrack2, gtrack, ht2),
    from = 13690000, to = 13750000
)))
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/Lactb2_cre.pdf", tis_gg2, width = 8, height = 6)

mouse_Kidney_H3K27ac_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Mus_musculus/chip2/Kidney/H3K27ac/chip/a84087e3-364e-4bed-94d9-5f8f48faba33/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    genome = "mm39", chr = "NC_000067.7",
    type = "polygon", name = "Mouse Kidney H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE
)
mouse_Kidney_H3K4me3_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Mus_musculus/chip2/Kidney/H3K4me3/chip/06ecdbe7-7413-4367-a219-2ab8600b4522/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    genome = "mm39", chr = "NC_000067.7",
    type = "polygon", name = "Mouse Kidney H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE
)
ot_mouse_kidney <- OverlayTrack(trackList = list(mouse_Kidney_H3K27ac_dr, mouse_Kidney_H3K4me3_dr), baseline = 2)
mouse_brain_H3K27ac_dr <- DataTrack(
    range = mouse_brain_H3K27ac_bw, genome = "mm39", chr = "NC_000067.7",
    type = "polygon", name = "Mouse Brain H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE
)
mouse_brain_H3K4me3_dr <- DataTrack(
    range = mouse_brain_H3K4me3_bw, genome = "mm39", chr = "NC_000067.7",
    type = "polygon", name = "Mouse Brain H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE
)
ot_mouse_brain <- OverlayTrack(trackList = list(mouse_brain_H3K27ac_dr, mouse_brain_H3K4me3_dr), baseline = 2)
mouse_liver_H3K27ac_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Mus_musculus/chip2/Liver/H3K27ac/chip/8c0d8db8-36d4-4513-ae18-148027a64b59/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    genome = "mm39", chr = "NC_000067.7",
    type = "polygon", name = "Mouse Liver H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE
)
mouse_liver_H3K4me3_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Mus_musculus/chip2/Liver/H3K4me3/chip/400549e4-6803-4407-b33e-2ce5c27f0b44/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    genome = "mm39", chr = "NC_000067.7",
    type = "polygon", name = "Mouse Liver H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE
)
ot_mouse_liver <- OverlayTrack(trackList = list(mouse_liver_H3K27ac_dr, mouse_liver_H3K4me3_dr), baseline = 2)

plotTracks(c(gtrack, ot_mouse_brain, ot_mouse_kidney, ot_mouse_liver),
    chr = "NC_000067.7",
    from = 13740000, to = 13742000,
    ylim = c(0, 10),
    sizes = c(1, 1, 1, 1), baseline = 2, legend = TRUE
)
mouse_peak_gg <- as.ggplot(grid2grob(plotTracks(c(gtrack, ot_mouse_brain, ot_mouse_kidney, ot_mouse_liver),
    chr = "NC_000067.7",
    from = 13740000, to = 13742000,
    ylim = c(0, 10),
    sizes = c(1, 1, 1, 1), baseline = 2, legend = TRUE
)))
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/Lactb2_peak.pdf", mouse_peak_gg, width = 8, height = 6)

load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/fpkm_ortho.RData")
# raw_FPKM_me_0_sep raw expression
# adj_FPKM_me_0_sep adjusted expression
# ortho_me_adj orthogroup information
ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)
# OG0009639 = Lactb2, OG0013300 = Xkr9
fpkm_df <- adj_FPKM_me_0_sep[which(adj_FPKM_me_0_sep$OGID %in% c("OG0009639", "OG0013300") & adj_FPKM_me_0_sep$species %in% c("Ovis_aries", "Equus_caballus", "Rhinolophus_pusillus", "Rhinopithecus_roxellana", "Sus_scrofa", "Felis_catus", "Mus_musculus")), ]
# turn into wide format
fpkm_df <- reshape2::melt(fpkm_df, id.vars = c("OGID", "species"), variable.name = "tissue", value.name = "fpkm")
# barplot()
fpkm_df$species <- factor(fpkm_df$species, levels = c("Ovis_aries", "Sus_scrofa", "Equus_caballus", "Felis_catus", "Rhinolophus_pusillus", "Mus_musculus"))
fpkm_bar <- ggplot(fpkm_df, aes(x = OGID, y = fpkm, fill = tissue)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    theme_bw() +
    scale_fill_manual(values = c("Brain" = "#FFD2A8", "Kidney" = "#7E8BB4", "Liver" = "#84C990")) +
    scale_x_discrete(labels = c("OG0009639" = "Lactb2", "OG0013300" = "Xkr9")) +
    facet_grid(species ~ ., scales = "fixed") +
    scale_y_continuous(breaks = c(0, 3, 6)) +
    theme(legend.position = "none") +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Gene orthologs", y = "log2(FPKM + 1)")
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/Lactb2_Xkr9_fpkm.pdf", fpkm_bar, width = 3, height = 6)
# Felis_catus_Kidney Felis_catus_Liver Equus_caballus_Liver Ovis_aries_Liver Sus_scrofa_Liver Rhinolophus_pusillus_Liver
# Rhinopithecus_roxellana_Liver
# Felis_catus, Equus_caballus Ovis_aries, Sus_scrofa, Rhinolophus_pusillus, Rhinopithecus_roxellana
# ec_k_12562 enh_k_130303_K_Peak_15875
mus_anno[30657, ]
mouse_func_group <- read.csv("/media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Mus_musculus_Kidney_enhancer_func_group.tsv", header = TRUE, sep = "\t")
mouse_func_group[mouse_func_group$Mus_musculus == "Peak_15875", ]

## Felis_catus orthologous region of Peak_15875
## grep Peak_15875 /media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Kidney/enhancer/Mus_musculus2Felis_catus_Kidney_enhancer_halper.bed4
## felCat.NC_058385.1      20835969        20836512        Peak_15875
## Felis_catus Kidney H3K27ac /DELL_EMC/zhangz/chip/analysis/Felis_catus/chip2/Kidney/H3K27ac/chip/3bc65576-6171-490d-a3f5-bc21254014d8/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Felis_catus Kidney H3K4me3 /DELL_EMC/zhangz/chip/analysis/Felis_catus/chip2/Kidney/H3K4me3/chip/4feb7e57-19c9-49da-b7fa-ac4cf7bdbaed/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Felis_catus Liver H3K27ac /DELL_EMC/zhangz/chip/analysis/Felis_catus/chip2/Liver/H3K27ac/chip/bc424a13-8c88-470c-b708-7650754a40dd/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Felis_catus Liver H3K4me3 /DELL_EMC/zhangz/chip/analysis/Felis_catus/chip2/Liver/H3K4me3/chip/fa3bb32d-0163-48b2-913a-0c94966d773f/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
Felis_Kidney_H3K27ac_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Felis_catus/chip2/Kidney/H3K27ac/chip/3bc65576-6171-490d-a3f5-bc21254014d8/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_058385.1",
    type = "polygon", name = "Cat Kidney H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE, ylim = c(0, 4)
)
Felis_Kidney_H3K4me3_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Felis_catus/chip2/Kidney/H3K4me3/chip/4feb7e57-19c9-49da-b7fa-ac4cf7bdbaed/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_058385.1",
    type = "polygon", name = "Cat Kidney H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE, ylim = c(0, 4)
)
ot_Felis_kidney <- OverlayTrack(trackList = list(Felis_Kidney_H3K27ac_dr, Felis_Kidney_H3K4me3_dr), baseline = 2, )
Felis_Liver_H3K27ac_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Felis_catus/chip2/Liver/H3K27ac/chip/bc424a13-8c88-470c-b708-7650754a40dd/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_058385.1",
    type = "polygon", name = "Cat Liver H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE, ylim = c(0, 6)
)
Felis_Liver_H3K4me3_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Felis_catus/chip2/Liver/H3K4me3/chip/fa3bb32d-0163-48b2-913a-0c94966d773f/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_058385.1",
    type = "polygon", name = "Cat Liver H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE, ylim = c(0, 6)
)
ot_Felis_Liver <- OverlayTrack(trackList = list(Felis_Liver_H3K27ac_dr, Felis_Liver_H3K4me3_dr), baseline = 2)
Felis_Kidney_Liver_ht <- HighlightTrack(
    trackList = list(gtrack, ot_Felis_kidney, ot_Felis_Liver),
    start = 20835969, width = 543, fill = "#f5dede", alpha = 0.2,
    chromosome = "NC_058385.1", sizes = c(1, 1.5, 1.5), inBackground = FALSE
)
plotTracks(Felis_Kidney_Liver_ht,
    chr = "NC_058385.1",
    from = 20835000, to = 20837000,
    baseline = 2, legend = TRUE
)
Felis_peak_gg <- as.ggplot(grid2grob(plotTracks(c(gtrack, ot_Felis_kidney, ot_Felis_Liver),
    chr = "NC_058385.1",
    from = 20835000, to = 20837000,
    ylim = c(0, 10),
    sizes = c(1, 1.5, 1.5), baseline = 2, legend = TRUE
)))
ggsave("/media/Data/zhangz/chip/analysis/summary2/ent_gene/Lactb2_Felis_peak.pdf", Felis_peak_gg, width = 8, height = 6)
## as.ggplot functions cannot properly process the plotTracks object, so we need to save the plot manually

## Equus_caballus orthologous region of Peak_15875
## grep Peak_15875 /media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Kidney/enhancer/Mus_musculus2Equus_caballus_Kidney_enhancer_halper.bed4
## equCab.NC_009152.3      15615443        15616136        Peak_15875
## Equus_caballus Liver H3K27ac /DELL_EMC/zhangz/chip/analysis/Equus_caballus/chip2/Liver/H3K27ac/chip/2556fffc-061b-4769-9b84-dfd4562d883d/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Equus_caballus Liver H3K4me3 /DELL_EMC/zhangz/chip/analysis/Equus_caballus/chip2/Liver/H3K4me3/chip/a27f5665-e050-468f-9569-19cec156ef0c/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig

Equus_Liver_H3K27ac_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Equus_caballus/chip2/Liver/H3K27ac/chip/2556fffc-061b-4769-9b84-dfd4562d883d/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_009152.3",
    type = "polygon", name = "Horse Liver H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE
)
Equus_Liver_H3K4me3_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Equus_caballus/chip2/Liver/H3K4me3/chip/a27f5665-e050-468f-9569-19cec156ef0c/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_009152.3",
    type = "polygon", name = "Horse Liver H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE
)
ot_Equus_Liver <- OverlayTrack(trackList = list(Equus_Liver_H3K27ac_dr, Equus_Liver_H3K4me3_dr), baseline = 2)
Equus_Liver_ht <- HighlightTrack(
    trackList = list(gtrack, ot_Equus_Liver),
    start = 15615443, width = 693, fill = "#f5dede", alpha = 0.2,
    chromosome = "NC_009152.3", sizes = c(1, 1.5), inBackground = FALSE
)
plotTracks(Equus_Liver_ht,
    chr = "NC_000067.7",
    from = 15615000, to = 15617000,
    ylim = c(0, 4), baseline = 2, legend = TRUE
)
# "/media/Data/zhangz/chip/analysis/summary2/ent_gene/Lactb2_Equus_peak.pdf"

## Ovis_aries orthologous region of Peak_15875
## grep Peak_15875 /media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Kidney/enhancer/Mus_musculus2Ovis_aries_Kidney_enhancer_halper.bed4
## oviAri.NC_056062.1      47331006        47331722        Peak_15875
## Ovis_aries Liver H3K27ac /DELL_EMC/zhangz/chip/analysis/Ovis_aries/chip2/Liver/H3K27ac/chip/c61255f8-bef3-43a8-9788-d65ff2bd83ee/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Ovis_aries Liver H3K4me3 /DELL_EMC/zhangz/chip/analysis/Ovis_aries/chip2/Liver/H3K4me3/chip/70471199-3f48-4152-90c4-e1e12f2914df/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
Ovis_Liver_H3K27ac <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Ovis_aries/chip2/Liver/H3K27ac/chip/c61255f8-bef3-43a8-9788-d65ff2bd83ee/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_056062.1",
    type = "polygon", name = "Sheep Liver H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE
)
Ovis_Liver_H3K4me3 <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Ovis_aries/chip2/Liver/H3K4me3/chip/70471199-3f48-4152-90c4-e1e12f2914df/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_056062.1",
    type = "polygon", name = "Sheep Liver H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE
)
ot_Ovis_Liver <- OverlayTrack(trackList = list(Ovis_Liver_H3K27ac, Ovis_Liver_H3K4me3), baseline = 2)
Ovis_Liver_ht <- HighlightTrack(
    trackList = list(gtrack, ot_Ovis_Liver),
    start = 47331006, width = 716, fill = "#f5dede", alpha = 0.2,
    chromosome = "NC_056062.1", sizes = c(1, 1.5), inBackground = FALSE
)
plotTracks(Ovis_Liver_ht,
    chr = "NC_056062.1",
    from = 47330000, to = 47332000,
    ylim = c(0, 6), baseline = 2, legend = FALSE
)


## Sus_scrofa orthologous region of Peak_15875
## grep Peak_15875 /media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Kidney/enhancer/Mus_musculus2Sus_scrofa_Kidney_enhancer_halper.bed4
## susScr.NC_010446.5      64702542        64703505        Peak_15875
## Sus_scrofa Liver H3K27ac /DELL_EMC/zhangz/chip/analysis/Sus_scrofa/chip2/Liver/H3K27ac/chip/1f4856cd-e0ac-4824-9412-99654066dcc8/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
## Sus_scrofa Liver H3K4me3 /DELL_EMC/zhangz/chip/analysis/Sus_scrofa/chip2/Liver/H3K4me3/chip/6aa3e6e2-c3ff-4aac-b854-f22d11d115d2/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig
Sus_Liver_H3K27ac_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Sus_scrofa/chip2/Liver/H3K27ac/chip/1f4856cd-e0ac-4824-9412-99654066dcc8/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_010446.5",
    type = "polygon", name = "Pig Liver H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE
)
Sus_Liver_H3K4me3_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Sus_scrofa/chip2/Liver/H3K4me3/chip/6aa3e6e2-c3ff-4aac-b854-f22d11d115d2/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig",
    chr = "NC_010446.5",
    type = "polygon", name = "Pig Liver H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE
)
ot_Sus_Liver <- OverlayTrack(trackList = list(Sus_Liver_H3K27ac_dr, Sus_Liver_H3K4me3_dr), baseline = 2)
Sus_Liver_ht <- HighlightTrack(
    trackList = list(gtrack, ot_Sus_Liver),
    start = 64702542, width = 963, fill = "#f5dede", alpha = 0.2,
    chromosome = "NC_010446.5", sizes = c(1, 1.5), inBackground = FALSE
)
plotTracks(Sus_Liver_ht,
    chr = "NC_010446.5",
    from = 64700000, to = 64705000,
    ylim = c(0, 10), baseline = 2, legend = FALSE
)

## Rhinolophus_pusillus orthologous region of Peak_15875
## grep Peak_15875 /media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Kidney/enhancer/Mus_musculus2Rhinolophus_pusillus_Kidney_enhancer_halper.bed4
## rhiPus.chr27_rhiPus     10110493        10111223        Peak_15875
## Rhinolophus_pusillus Liver H3K27ac /DELL_EMC/zhangz/chip/analysis/Rhinolophus_pusillus/chip2/Liver/H3K27ac/chip/3e189a1a-7445-4582-920a-b8150c6bf36e/call-macs2_signal_track/shard-0/execution/CR-8-1-P_1_val_1.srt.nodup_x_CR-8-3-I_1_val_1.srt.nodup.fc.signal.bigwig
## Rhinolophus_pusillus Liver H3K4me3 /DELL_EMC/zhangz/chip/analysis/Rhinolophus_pusillus/chip2/Liver/H3K4me3/chip/42fd23dd-7aba-477f-a967-9ff0fe955efb/call-macs2_signal_track/shard-0/execution/CR-8-2-P_1_val_1.srt.nodup_x_CR-8-3-I_1_val_1.srt.nodup.fc.signal.bigwig
Rhinolophus_Liver_H3K27ac_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Rhinolophus_pusillus/chip2/Liver/H3K27ac/chip/3e189a1a-7445-4582-920a-b8150c6bf36e/call-macs2_signal_track/shard-0/execution/CR-8-1-P_1_val_1.srt.nodup_x_CR-8-3-I_1_val_1.srt.nodup.fc.signal.bigwig",
    chr = "chr27_rhiPus",
    type = "polygon", name = "Bat Liver H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE
)
Rhinolophus_Liver_H3K4me3_dr <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Rhinolophus_pusillus/chip2/Liver/H3K4me3/chip/42fd23dd-7aba-477f-a967-9ff0fe955efb/call-macs2_signal_track/shard-0/execution/CR-8-2-P_1_val_1.srt.nodup_x_CR-8-3-I_1_val_1.srt.nodup.fc.signal.bigwig",
    chr = "chr27_rhiPus",
    type = "polygon", name = "Bat Liver H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE
)
ot_Rhinolophus_Liver <- OverlayTrack(trackList = list(Rhinolophus_Liver_H3K27ac_dr, Rhinolophus_Liver_H3K4me3_dr), baseline = 2)
ht_Rhinolophus_liver <- HighlightTrack(
    trackList = list(gtrack, ot_Rhinolophus_Liver),
    start = 10110493, width = 730, fill = "#f5dede", alpha = 0.2,
    chromosome = "chr27_rhiPus", sizes = c(1, 1.5), inBackground = FALSE
)
plotTracks(ht_Rhinolophus_liver,
    chr = "chr27_rhiPus",
    from = 10110000, to = 10112000,
    ylim = c(0, 10), baseline = 2, legend = FALSE
)

## Rhinopithecus_roxellana orthologous region of Peak_15875
## grep Peak_15875 /media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Kidney/enhancer/Mus_musculus2Rhinopithecus_roxellana_Kidney_enhancer_halper.bed4
## rhiRox.NC_044557.1      130283570       130284028       Peak_15875
## Rhinopithecus_roxellana Liver H3K27ac /DELL_EMC/zhangz/chip/analysis/Rhinopithecus_roxellana/chip2/Liver/H3K27ac/chip/fda39ca3-f652-4a0c-b7af-1364bb085585/call-macs2_signal_track/shard-0/execution/CR.41.1.P_1_val_1.srt.nodup_x_CR.41.3.I_1_val_1.srt.nodup.fc.signal.bigwig
## Rhinopithecus_roxellana Liver H3K4me3 /DELL_EMC/zhangz/chip/analysis/Rhinopithecus_roxellana/chip2/Liver/H3K4me3/chip/7a1f3b45-020a-4d17-a216-9f518c5f5c4a/call-macs2_signal_track/shard-0/execution/CR.41.2.P_1_val_1.srt.nodup_x_CR.41.3.I_1_val_1.srt.nodup.fc.signal.bigwig
dr_Rhinopithecus_Liver_H3K27ac <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Rhinopithecus_roxellana/chip2/Liver/H3K27ac/chip/fda39ca3-f652-4a0c-b7af-1364bb085585/call-macs2_signal_track/shard-0/execution/CR.41.1.P_1_val_1.srt.nodup_x_CR.41.3.I_1_val_1.srt.nodup.fc.signal.bigwig",
    chr = "NC_044557.1",
    type = "polygon", name = "Golden snub-nosed monkey Liver H3K27ac", alpha = 0.9,
    fill.mountain = c("white", "#ef4343"), col.mountain = "#ef4343",
    legend = TRUE
)
dr_Rhinopithecus_Liver_H3K4me3 <- DataTrack(
    range = "/DELL_EMC/zhangz/chip/analysis/Rhinopithecus_roxellana/chip2/Liver/H3K4me3/chip/7a1f3b45-020a-4d17-a216-9f518c5f5c4a/call-macs2_signal_track/shard-0/execution/CR.41.2.P_1_val_1.srt.nodup_x_CR.41.3.I_1_val_1.srt.nodup.fc.signal.bigwig",
    chr = "NC_044557.1",
    type = "polygon", name = "Monkey Liver H3K4me3", alpha = 0.9,
    fill.mountain = c("white", "#73b8d5"), col.mountain = "#73b8d5",
    legend = TRUE
)
ot_Rhinopithecus_Liver <- OverlayTrack(trackList = list(dr_Rhinopithecus_Liver_H3K27ac, dr_Rhinopithecus_Liver_H3K4me3), baseline = 2)
ht_Rhinopithecus_liver <- HighlightTrack(
    trackList = list(gtrack, ot_Rhinopithecus_Liver),
    start = 130283570, width = 458, fill = "#f5dede", alpha = 0.2,
    chromosome = "NC_044557.1", sizes = c(1, 1.5), inBackground = FALSE
)
plotTracks(ht_Rhinopithecus_liver,
    chr = "NC_044557.1",
    from = 130283000, to = 130285000,
    ylim = c(0, 6), baseline = 2, legend = FALSE
)

## load data
sp <- "Mus_musculus"
mus_ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
mus_ent <- mus_ent[mus_ent$ent1 > 0, ]
summary(mus_ent)

ent_element_num_box <- ggplot(mus_ent, aes(x = unique_ele_num, y = ent1, group = unique_ele_num)) +
    geom_boxplot(outlier.size = 0.1, linewidth = 0.3) +
    theme_classic() +
    labs(x = "Unique CRE number", y = "Regulatory entropy") +
    # scale_x_binned(breaks = c(2, 4, 6)) +
    scale_y_continuous(breaks = c(0.25, 0.75, 1.25)) +
    theme(axis.title = element_text(size = 8, face = "bold")) +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(linewidth = 0.2))

ggsave("/media/Data/zhangz/chip/analysis/summary2/pgls/ent_cre_num.pdf", ent_element_num_box, width = 6, height = 6, units = "cm")
ent_element_num_hex <- ggplot(mus_ent, aes(x = unique_ele_num, y = ent1, group = unique_ele_num)) +
    geom_hex() +
    theme_bw()

# plot density point plot to show the relationship between regulatory entropies and gene expression
ent_density <- ggplot(mus_ent, aes(x = ent1, y = fpkm_cv)) +
    geom_hex(bins = 70) +
    theme_classic() +
    labs(x = "Regulatory entropy", y = "Gene expression CV") +
    scale_x_continuous(breaks = c(0.25, 0.75, 1.25)) +
    scale_y_continuous(breaks = c(0, 0.75, 1.5)) +
    scale_fill_continuous(low = "#a2a6dd", high = "#1010a5") +
    theme(axis.title = element_text(size = 8, face = "bold")) +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(linewidth = 0.2)) +
    theme(legend.title = element_text(size = 6)) +
    theme(legend.text = element_text(size = 5)) +
    theme(legend.box.margin = margin(1, 1, 1, 1)) +
    theme(legend.key.size = unit(0.25, "cm")) +
    theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8))
ggsave("/media/Data/zhangz/chip/analysis/summary2/pgls/ent_fpkm_cv.pdf", ent_density, width = 6, height = 6, units = "cm")
cor.test(mus_ent$ent1, mus_ent$fpkm_cv, method = "spearman")
cor.test(mus_ent$ent1, mus_ent$adj_fpkm_cv, method = "spearman")

head(cre_info[which(cre_info$specificity == "evol_pleio_e" & cre_info$ec_level == "sp" & cre_info$act_sp > 10), ])
# the second one looks good, chr: NC_000070.7, start: 57130947, end: 57131131
cre_info[cre_info$chr == "NC_000070.7" & cre_info$start > 57130000 & cre_info$end < 57132000, ]
# there is also a kidney enhancer here, dont know why it is not count as pleiotropic

# NC_000080.7 121337046 121337573
cre_info[cre_info$chr == "NC_000080.7" & cre_info$start > 121336000 & cre_info$end < 121339000, ]
# there is a kidney enhancer here, cover the brain enhancer

# NC_000075.7 63603300  63603583
cre_info[cre_info$chr == "NC_000075.7" & cre_info$start > 63602000 & cre_info$end < 63605000, ]
all_cre <- read.csv("/media/Data/zhangz/chip/analysis/summary2/anno_sum/Mus_musculus/Mus_musculus_all.csv", header = TRUE, sep = ",")
all_cre[all_cre$chr == "NC_000075.7" & all_cre$start > 63602000 & all_cre$end < 63605000, ]
cre_info[541, ]
# enh_b_102049_B_Peak_88067
mouse_brain_func_group <- read.csv("/media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Mus_musculus_Brain_enhancer_func_group.tsv", header = TRUE, sep = "\t")
mouse_brain_func_group[mouse_brain_func_group$Mus_musculus == "Peak_88067", ]
