library(PerformanceAnalytics)
library(corrplot)

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

for (sp in species) {
    sp_all <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_all.csv"),
        header = TRUE
    )
    sp_enh <- sp_all[sp_all$element == "enhancer", ]
    sp_pro <- sp_all[sp_all$element == "promoter", ]
    enh_rs <- cor(sp_enh[, c("gc", "align", "consrv", "mean_phylop", "mean_fc", "motif_num", "length")], method = "spearman")
    enh_pval <- cor.mtest(sp_enh[, c("gc", "align", "consrv", "mean_phylop", "mean_fc", "motif_num", "length")], method = "spearman")
    png(paste0("/media/Data/zhangz/chip/analysis/summary2/correlation/", sp, "_enh_cor.png"), width = 8000, height = 10000, res = 600)
    corrplot(enh_rs, p.mat = enh_pval$p, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)
    dev.off()
    write.csv(enh_rs, paste0("/media/Data/zhangz/chip/analysis/summary2/correlation/", sp, "_enh_cor_rho.csv"), row.names = TRUE)
    write.csv(enh_pval$p, paste0("/media/Data/zhangz/chip/analysis/summary2/correlation/", sp, "_enh_cor_pval.csv"), row.names = TRUE)
    pro_rs <- cor(sp_pro[, c("gc", "align", "consrv", "mean_phylop", "mean_fc", "motif_num", "length")], method = "spearman")
    pro_pval <- cor.mtest(sp_pro[, c("gc", "align", "consrv", "mean_phylop", "mean_fc", "motif_num", "length")], method = "spearman")
    png(paste0("/media/Data/zhangz/chip/analysis/summary2/correlation/", sp, "_pro_cor.png"), width = 8000, height = 10000, res = 600)
    corrplot(pro_rs, p.mat = pro_pval$p, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)
    dev.off()
    write.csv(pro_rs, paste0("/media/Data/zhangz/chip/analysis/summary2/correlation/", sp, "_pro_cor_rho.csv"), row.names = TRUE)
    write.csv(pro_pval$p, paste0("/media/Data/zhangz/chip/analysis/summary2/correlation/", sp, "_pro_cor_pval.csv"), row.names = TRUE)
}
