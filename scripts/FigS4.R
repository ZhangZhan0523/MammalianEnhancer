## draw plots for fig S4
library(ggplot2)

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
ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)

# read data
pgls_res <- read.csv("/media/Data/zhangz/chip/analysis/summary2/pgls/pgls_res.csv", header = TRUE)
sig_pgls_res <- pgls_res[pgls_res[["Max_P.value"]] < 0.05, ]
sig_pgls_res <- sig_pgls_res[!is.na(sig_pgls_res$OGID), ]
sig_pgls_res <- sig_pgls_res %>% arrange(trait, P.value)
sig_pgls_res <- sig_pgls_res[, 3:11]
sig_pgls_res$gene <- sapply(sig_pgls_res$OGID, find_mus_ortho)
sig_pgls_res <- sig_pgls_res[!is.na(sig_pgls_res$Trait), ]
sig_pgls_res <- sig_pgls_res[which(!sig_pgls_res$Trait %in% c("CpG_ratio_Brain", "CpG_ratio_Kidney", "CpG_ratio_Liver", "CpGcounts")), ]
upset_matrix <- table(sig_pgls_res$OGID, sig_pgls_res$Trait)
upset_df <- as.data.frame(upset_matrix)
colnames(upset_df) <- c("OGID", "Trait", "Frequency")

# long df upset_df to wide df
wide_df <- spread(upset_df, Trait, Frequency)
# pdf("/media/Data/zhangz/chip/analysis/summary2/pgls/upset.pdf", width = 10, height = 10)
wide_df <- wide_df[, c(
    "OGID", "ML.yrs.y", "MLres.y", "FTM.d.y", "FTMres.y",
    "AW.g", "Birth_weight.g.", "Litter_size", "Weaning.days.",
    "Gestation.days.", "Reproduction_frequency", "Inter.litter_interval",
    "body_colume", "Brain_volume", "Brain_ratio", "BMR_mass.g.", "BMR.O2ml.h.", "adj.BMR.O2.ml.h.g.",
    "Temperature.K.", "Genome.size"
)]
require(UpSetR)
p1 <- upset(wide_df,
    # sets = c("ML.yrs.y", "FTM.d.y", "MLres.y", "FTMres.y", "AW.g"),
    nset = 19, nintersects = 60,
    keep.order = TRUE,
    sets = c(
        "ML.yrs.y", "MLres.y", "FTM.d.y", "FTMres.y",
        "AW.g", "Birth_weight.g.", "Litter_size", "Weaning.days.",
        "Gestation.days.", "Reproduction_frequency", "Inter.litter_interval",
        "body_colume", "Brain_volume", "Brain_ratio", "BMR_mass.g.", "BMR.O2ml.h.", "adj.BMR.O2.ml.h.g.",
        "Temperature.K.", "Genome.size"
    ),
    queries = list(
        list(
            query = intersects, params = list("ML.yrs.y", "MLres.y"),
            color = "red", active = TRUE
        ),
        list(
            query = intersects, params = list("MLres.y", "ML.yrs.y", "Litter_size"),
            color = "red", active = TRUE
        ),
        list(
            query = intersects, params = list("MLres.y", "ML.yrs.y", "Litter_size", "Weaning.days."),
            color = "red", active = TRUE
        )
    ),
    mb.ratio = c(0.5, 0.5),
    text.scale = c(2, 3, 2, 1.5, 2, 2)
)
# require(ggplotify)
# g1 <- as.ggplot(p1)
# ggsave("/media/Data/zhangz/chip/analysis/summary2/pgls/upset.pdf", g1, width = 10, height = 10, units = "cm")
pdf("/media/Data/zhangz/chip/analysis/summary2/pgls/upset_all.pdf", width = 16, height = 9)
print(p1)
dev.off()
wide_df[which(wide_df$OGID %in% c("OG0009639", "OG0006213", "OG0007610")), ]

#          OGID ML.yrs.y MLres.y FTM.d.y FTMres.y AW.g Birth_weight.g. Litter_size Weaning.days. Gestation.days. Reproduction_frequency Inter.litter_interval body_colume Brain_volume Brain_ratio BMR_mass.g. BMR.O2ml.h. adj.BMR.O2.ml.h.g. Temperature.K. Genome.size
# 44  OG0006213        1       1       0        0    0               0           1             1               0                      0                     0           0            0           0           0           0                  0              0           0
# 69  OG0007610        1       1       0        0    0               0           0             0               0                      0                     0           0            0           0           0           0                  0              0           0
# 102 OG0009639        1       1       0        0    0               0           1             0               0                      0                     0           0            0           0           0           0                  0              0           0
sig_pgls_res[which(sig_pgls_res$OGID %in% c("OG0009639", "OG0006213", "OG0007610")), ]
#          OGID         Trait Model      P.value Max_P.value           Species        R2      Slope log_trait                  gene
# 242 OG0006213   Litter_size    OU 0.0021599063 0.029533118 Rattus_norvegicus 0.4543085   8.319696     FALSE Pafah2;NM_001285872.1
# 246 OG0009639   Litter_size    OU 0.0065672061 0.035726933 Rattus_norvegicus 0.3985763   6.129171     FALSE    Lactb2;NM_145381.2
# 249 OG0006213      ML.yrs.y    OU 0.0014268592 0.007930857 Rattus_norvegicus 0.4803992 -33.785652     FALSE Pafah2;NM_001285872.1
# 251 OG0009639      ML.yrs.y    OU 0.0031272911 0.011802283      Mus_musculus 0.4516265 -34.146711     FALSE    Lactb2;NM_145381.2
# 261 OG0007610      ML.yrs.y    OU 0.0116318796 0.047457037 Rattus_norvegicus 0.3983212 -40.135788     FALSE     Abcf3;NM_013852.2
# 267 OG0009639       MLres.y    OU 0.0004782305 0.001518953 Rattus_norvegicus 0.5677678  -1.941117     FALSE    Lactb2;NM_145381.2
# 273 OG0006213       MLres.y    OU 0.0097101965 0.024729596 Rattus_norvegicus 0.3499585  -2.019702     FALSE Pafah2;NM_001285872.1
# 274 OG0007610       MLres.y    OU 0.0105678017 0.028730282 Rattus_norvegicus 0.4064363  -2.965524     FALSE     Abcf3;NM_013852.2
# 342 OG0006213 Weaning.days.    OU 0.0148069462 0.038817024      Equus_asinus 0.3179075  -2.348377      TRUE Pafah2;NM_001285872.1

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
ent_cv_density <- ggplot(mus_ent, aes(x = ent1, y = fpkm_cv)) +
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
    theme(legend.box.margin =  margin(1, 1, 1, 1)) +
    theme(legend.key.size = unit(0.25, "cm")) +
    theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8))
ggsave("/media/Data/zhangz/chip/analysis/summary2/pgls/ent_fpkm_cv.pdf", ent_density, width = 6, height = 6, units = "cm")
cor.test(mus_ent$ent1, mus_ent$fpkm_cv, method = "spearman")
cor.test(mus_ent$ent1, mus_ent$adj_fpkm_cv, method = "spearman")

# plot density point plot to show the relationship between regulatory entropies and gene expression
ent_fpkm_density <- ggplot(mus_ent, aes(x = ent1, y = log2(fpkm + 1))) +
    geom_hex(bins = 70) +
    theme_classic() +
    labs(x = "Regulatory entropy", y = "Gene expression") +
    scale_x_continuous(breaks = c(0.25, 0.75, 1.25)) +
    scale_y_continuous(breaks = c(0, 0.75, 1.5)) +
    scale_fill_continuous(low = "#a2a6dd", high = "#1010a5") +
    theme(axis.title = element_text(size = 8, face = "bold")) +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.2)) +
    theme(axis.ticks = element_line(linewidth = 0.2)) +
    theme(legend.title = element_text(size = 6)) +
    theme(legend.text = element_text(size = 5)) +
    theme(legend.box.margin =  margin(1, 1, 1, 1)) +
    theme(legend.key.size = unit(0.25, "cm")) +
    theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8))
cor.test(mus_ent$ent1, mus_ent$fpkm, method = "spearman")
