# load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/ortho_me_adj.RData")
library(optparse)
option_list <- list(
    make_option(
        c("-s", "--species"),
        type = "character", default = "/media/Data/zhangz/chip/scripts/info/info_using.csv",
        help = "species need to summarize"
    )
)
opt_parser <- parse_args(OptionParser(option_list = option_list, usage = "%prog [options]"))
sp <- opt_parser$species
print(sp)
library(magrittr)
# st_chisq <- data.frame(
#     species = rep(unique(raw_FPKM_me_0_sep$species), 3),
#     tissue = c(rep("Brain", 21), rep("Kidney", 21), rep("Liver", 21)),
#     p_ex = 0, np_ex = 0, p_nex = 0, np_nex = 0,
#     e_ex = 0, ne_ex = 0, e_nex = 0, ne_nex = 0,
#     p_ex_chiq_p = 0, p_ex_d = 0, e_ex_chiq_p = 0, e_ex_d = 0,
#     enhancer_num_spm_p = 0, enhancer_num_spm_rho = 0,
#     promoter_num_spm_p = 0, promoter_num_spm_rho = 0,
#     enhancer_fc_spm_p = 0, enhancer_fc_spm_rho = 0,
#     enhancer_gc_spm_p = 0, enhancer_gc_spm_rho = 0,
#     enhancer_length_spm_p = 0, enhancer_length_spm_rho = 0,
#     enhancer_align_spm_p = 0, enhancer_align_spm_rho = 0,
#     enhancer_con_spm_p = 0, enhancer_con_spm_rho = 0,
#     enhancer_tss_spm_p = 0, enhancer_tss_spm_rho = 0,
#     enhancer_motif_num_spm_p = 0, enhancer_motif_num_spm_rho = 0,
#     promoter_fc_spm_p = 0, promoter_fc_spm_rho = 0,
#     promoter_gc_spm_p = 0, promoter_gc_spm_rho = 0,
#     promoter_length_spm_p = 0, promoter_length_spm_rho = 0,
#     promoter_align_spm_p = 0, promoter_align_spm_rho = 0,
#     promoter_con_spm_p = 0, promoter_con_spm_rho = 0,
#     promoter_tss_spm_p = 0, promoter_tss_spm_rho = 0,
#     promoter_motif_num_spm_p = 0, promoter_motif_num_spm_rho = 0
# )
# sp_spm <- data.frame(
#     species = unique(raw_FPKM_me_0_sep$species),
#     p_ex = 0, np_ex = 0, p_nex = 0, np_nex = 0,
#     e_ex = 0, ne_ex = 0, e_nex = 0, ne_nex = 0,
#     p_ex_chiq_p = 0, p_ex_d = 0, e_ex_chiq_p = 0, e_ex_d = 0,
#     enhancer_num_spm_p = 0, enhancer_num_spm_rho = 0,
#     promoter_num_spm_p = 0, promoter_num_spm_rho = 0,
#     enhancer_fc_spm_p = 0, enhancer_fc_spm_rho = 0,
#     enhancer_gc_spm_p = 0, enhancer_gc_spm_rho = 0,
#     enhancer_length_spm_p = 0, enhancer_length_spm_rho = 0,
#     enhancer_align_spm_p = 0, enhancer_align_spm_rho = 0,
#     enhancer_con_spm_p = 0, enhancer_con_spm_rho = 0,
#     enhancer_tss_spm_p = 0, enhancer_tss_spm_rho = 0,
#     enhancer_motif_num_spm_p = 0, enhancer_motif_num_spm_rho = 0,
#     promoter_fc_spm_p = 0, promoter_fc_spm_rho = 0,
#     promoter_gc_spm_p = 0, promoter_gc_spm_rho = 0,
#     promoter_length_spm_p = 0, promoter_length_spm_rho = 0,
#     promoter_align_spm_p = 0, promoter_align_spm_rho = 0,
#     promoter_con_spm_p = 0, promoter_con_spm_rho = 0,
#     promoter_tss_spm_p = 0, promoter_tss_spm_rho = 0,
#     promoter_motif_num_spm_p = 0, promoter_motif_num_spm_rho = 0
# )
# sp <- "Rhinolophus_pusillus"
sp_all <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_all.csv"),
    header = TRUE
)
load("/media/Data/zhangz/chip/analysis/summary2/sum_all/qn_FPKM_0.RData")
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/fpkm_ortho.RData")

raw_FPKM_me_0_sep[raw_FPKM_me_0_sep$species == "Macaca_fascicularis", "species"] <- "Macaca_mulatta"
colnames(ortho_me_adj)[which(colnames(ortho_me_adj) == "Macaca_fascicularis")] <- "Macaca_mulatta"
raw_FPKM_sp <- raw_FPKM_me_0_sep[which(raw_FPKM_me_0_sep$species == sp & !is.na(raw_FPKM_me_0_sep$Brain)), ]
sp_ortho <- ortho_me_adj[ortho_me_adj[sp] != "NULL", c("id", sp)]
# gene_all <- enh num, pro num, mean enh fc, max_enh_fc, median enh fc, mean pro fc, max pro fc, median pro fc
# mean enh motif num, max enh motif num, median enh motif num, mean pro motif num, max pro motif num, median pro motif num
brain_og_sum <- data.frame()
kidney_og_sum <- data.frame()
liver_og_sum <- data.frame()
for (og in sp_ortho$id) {
    tmp <- sp_all[grep(og, sp_all$anno_orthogroup), c("peak", "mean_fc", "motif_num", "mean_phylop", "anno_distance", "element", "tissue", "lib_fc", "align", "consrv")]
    b_tmp <- tmp[tmp$tissue == "Brain", ]
    b_sp_fpkm_mean <- mean(raw_FPKM_me_0_sep[raw_FPKM_me_0_sep$OGID == og, "Brain"], na.rm = TRUE)
    k_sp_fpkm_mean <- mean(raw_FPKM_me_0_sep[raw_FPKM_me_0_sep$OGID == og, "Kidney"], na.rm = TRUE)
    l_sp_fpkm_mean <- mean(raw_FPKM_me_0_sep[raw_FPKM_me_0_sep$OGID == og, "Liver"], na.rm = TRUE)
    b_fpkm_diff <- raw_FPKM_sp[which(raw_FPKM_sp$OGID == og), "Brain"] - b_sp_fpkm_mean
    k_fpkm_diff <- raw_FPKM_sp[which(raw_FPKM_sp$OGID == og), "Kidney"] - k_sp_fpkm_mean
    l_fpkm_diff <- raw_FPKM_sp[which(raw_FPKM_sp$OGID == og), "Liver"] - l_sp_fpkm_mean
    brain_tmp <- data.frame(
        OGID = og, fpkm = raw_FPKM_sp[which(raw_FPKM_sp$OGID == og), "Brain"],
        fpkm_diff = b_fpkm_diff,
        mean_enh_align = mean(b_tmp[which(b_tmp$element == "enhancer"), "align"]),
        mean_pro_align = mean(b_tmp[which(b_tmp$element == "promoter"), "align"]),
        min_enh_align = min(b_tmp[which(b_tmp$element == "enhancer"), "align"]),
        min_pro_align = min(b_tmp[which(b_tmp$element == "promoter"), "align"]),
        max_enh_align = max(b_tmp[which(b_tmp$element == "enhancer"), "align"]),
        max_pro_align = max(b_tmp[which(b_tmp$element == "promoter"), "align"]),
        median_enh_align = median(b_tmp[which(b_tmp$element == "enhancer"), "align"]),
        median_pro_align = median(b_tmp[which(b_tmp$element == "promoter"), "align"]),
        mean_enh_con = mean(b_tmp[which(b_tmp$element == "enhancer"), "consrv"]),
        mean_pro_con = mean(b_tmp[which(b_tmp$element == "promoter"), "consrv"]),
        min_enh_con = min(b_tmp[which(b_tmp$element == "enhancer"), "consrv"]),
        min_pro_con = min(b_tmp[which(b_tmp$element == "promoter"), "consrv"]),
        max_enh_con = max(b_tmp[which(b_tmp$element == "enhancer"), "consrv"]),
        max_pro_con = max(b_tmp[which(b_tmp$element == "promoter"), "consrv"]),
        median_enh_con = median(b_tmp[which(b_tmp$element == "enhancer"), "consrv"]),
        median_pro_con = median(b_tmp[which(b_tmp$element == "promoter"), "consrv"]),
        enhancer_num = nrow(b_tmp[which(b_tmp$element == "enhancer"), ]),
        promoter_num = nrow(b_tmp[which(b_tmp$element == "promoter"), ]),
        mean_enh_fc = mean(b_tmp[which(b_tmp$element == "enhancer"), "lib_fc"]),
        max_enh_fc = max(b_tmp[which(b_tmp$element == "enhancer"), "lib_fc"]),
        median_enh_fc = median(b_tmp[which(b_tmp$element == "enhancer"), "lib_fc"]),
        mean_pro_fc = mean(b_tmp[which(b_tmp$element == "promoter"), "lib_fc"]),
        max_pro_fc = max(b_tmp[which(b_tmp$element == "promoter"), "lib_fc"]),
        median_pro_fc = median(b_tmp[which(b_tmp$element == "promoter"), "lib_fc"]),
        mean_enh_motif_num = mean(b_tmp[which(b_tmp$element == "enhancer"), "motif_num"]),
        max_enh_motif_num = max(b_tmp[which(b_tmp$element == "enhancer"), "motif_num"]),
        median_enh_motif_num = median(b_tmp[which(b_tmp$element == "enhancer"), "motif_num"]),
        mean_pro_motif_num = mean(b_tmp[which(b_tmp$element == "promoter"), "motif_num"]),
        max_pro_motif_num = max(b_tmp[which(b_tmp$element == "promoter"), "motif_num"]),
        median_pro_motif_num = median(b_tmp[which(b_tmp$element == "promoter"), "motif_num"]),
        mean_enh_phylop = mean(b_tmp[which(b_tmp$element == "enhancer"), "mean_phylop"]),
        mean_pro_phylop = mean(b_tmp[which(b_tmp$element == "promoter"), "mean_phylop"]),
        median_enh_phylop = median(b_tmp[which(b_tmp$element == "enhancer"), "mean_phylop"]),
        median_pro_phylop = median(b_tmp[which(b_tmp$element == "promoter"), "mean_phylop"])
    )
    brain_og_sum <- rbind(brain_og_sum, brain_tmp)
    k_tmp <- tmp[tmp$tissue == "Kidney", ]
    kidney_tmp <- data.frame(
        OGID = og, fpkm = raw_FPKM_sp[which(raw_FPKM_sp$OGID == og), "Kidney"],
        fpkm_diff = k_fpkm_diff,
        mean_enh_align = mean(k_tmp[which(k_tmp$element == "enhancer"), "align"]),
        mean_pro_align = mean(k_tmp[which(k_tmp$element == "promoter"), "align"]),
        min_enh_align = min(k_tmp[which(k_tmp$element == "enhancer"), "align"]),
        min_pro_align = min(k_tmp[which(k_tmp$element == "promoter"), "align"]),
        max_enh_align = max(k_tmp[which(k_tmp$element == "enhancer"), "align"]),
        max_pro_align = max(k_tmp[which(k_tmp$element == "promoter"), "align"]),
        median_enh_align = median(k_tmp[which(k_tmp$element == "enhancer"), "align"]),
        median_pro_align = median(k_tmp[which(k_tmp$element == "promoter"), "align"]),
        mean_enh_con = mean(k_tmp[which(k_tmp$element == "enhancer"), "consrv"]),
        mean_pro_con = mean(k_tmp[which(k_tmp$element == "promoter"), "consrv"]),
        min_enh_con = min(k_tmp[which(k_tmp$element == "enhancer"), "consrv"]),
        min_pro_con = min(k_tmp[which(k_tmp$element == "promoter"), "consrv"]),
        max_enh_con = max(k_tmp[which(k_tmp$element == "enhancer"), "consrv"]),
        max_pro_con = max(k_tmp[which(k_tmp$element == "promoter"), "consrv"]),
        median_enh_con = median(k_tmp[which(k_tmp$element == "enhancer"), "consrv"]),
        median_pro_con = median(k_tmp[which(k_tmp$element == "promoter"), "consrv"]),
        enhancer_num = nrow(k_tmp[which(k_tmp$element == "enhancer"), ]),
        promoter_num = nrow(k_tmp[which(k_tmp$element == "promoter"), ]),
        mean_enh_fc = mean(k_tmp[which(k_tmp$element == "enhancer"), "lib_fc"]),
        max_enh_fc = max(k_tmp[which(k_tmp$element == "enhancer"), "lib_fc"]),
        median_enh_fc = median(k_tmp[which(k_tmp$element == "enhancer"), "lib_fc"]),
        mean_pro_fc = mean(k_tmp[which(k_tmp$element == "promoter"), "lib_fc"]),
        max_pro_fc = max(k_tmp[which(k_tmp$element == "promoter"), "lib_fc"]),
        median_pro_fc = median(k_tmp[which(k_tmp$element == "promoter"), "lib_fc"]),
        mean_enh_motif_num = mean(k_tmp[which(k_tmp$element == "enhancer"), "motif_num"]),
        max_enh_motif_num = max(k_tmp[which(k_tmp$element == "enhancer"), "motif_num"]),
        median_enh_motif_num = median(k_tmp[which(k_tmp$element == "enhancer"), "motif_num"]),
        mean_pro_motif_num = mean(k_tmp[which(k_tmp$element == "promoter"), "motif_num"]),
        max_pro_motif_num = max(k_tmp[which(k_tmp$element == "promoter"), "motif_num"]),
        median_pro_motif_num = median(k_tmp[which(k_tmp$element == "promoter"), "motif_num"]),
        mean_enh_phylop = mean(k_tmp[which(k_tmp$element == "enhancer"), "mean_phylop"]),
        mean_pro_phylop = mean(k_tmp[which(k_tmp$element == "promoter"), "mean_phylop"]),
        median_enh_phylop = median(k_tmp[which(k_tmp$element == "enhancer"), "mean_phylop"]),
        median_pro_phylop = median(k_tmp[which(k_tmp$element == "promoter"), "mean_phylop"])
    )
    kidney_og_sum <- rbind(kidney_og_sum, kidney_tmp)
    liver_tmp <- tmp[tmp$tissue == "Liver", ]
    liver_tmp <- data.frame(
        OGID = og, fpkm = raw_FPKM_sp[which(raw_FPKM_sp$OGID == og), "Liver"],
        fpkm_diff = l_fpkm_diff,
        mean_enh_align = mean(liver_tmp[which(liver_tmp$element == "enhancer"), "align"]),
        mean_pro_align = mean(liver_tmp[which(liver_tmp$element == "promoter"), "align"]),
        min_enh_align = min(liver_tmp[which(liver_tmp$element == "enhancer"), "align"]),
        min_pro_align = min(liver_tmp[which(liver_tmp$element == "promoter"), "align"]),
        max_enh_align = max(liver_tmp[which(liver_tmp$element == "enhancer"), "align"]),
        max_pro_align = max(liver_tmp[which(liver_tmp$element == "promoter"), "align"]),
        median_enh_align = median(liver_tmp[which(liver_tmp$element == "enhancer"), "align"]),
        median_pro_align = median(liver_tmp[which(liver_tmp$element == "promoter"), "align"]),
        mean_enh_con = mean(liver_tmp[which(liver_tmp$element == "enhancer"), "consrv"]),
        mean_pro_con = mean(liver_tmp[which(liver_tmp$element == "promoter"), "consrv"]),
        min_enh_con = min(liver_tmp[which(liver_tmp$element == "enhancer"), "consrv"]),
        min_pro_con = min(liver_tmp[which(liver_tmp$element == "promoter"), "consrv"]),
        max_enh_con = max(liver_tmp[which(liver_tmp$element == "enhancer"), "consrv"]),
        max_pro_con = max(liver_tmp[which(liver_tmp$element == "promoter"), "consrv"]),
        median_enh_con = median(liver_tmp[which(liver_tmp$element == "enhancer"), "consrv"]),
        median_pro_con = median(liver_tmp[which(liver_tmp$element == "promoter"), "consrv"]),
        enhancer_num = nrow(liver_tmp[which(liver_tmp$element == "enhancer"), ]),
        promoter_num = nrow(liver_tmp[which(liver_tmp$element == "promoter"), ]),
        mean_enh_fc = mean(liver_tmp[which(liver_tmp$element == "enhancer"), "lib_fc"]),
        max_enh_fc = max(liver_tmp[which(liver_tmp$element == "enhancer"), "lib_fc"]),
        median_enh_fc = median(liver_tmp[which(liver_tmp$element == "enhancer"), "lib_fc"]),
        mean_pro_fc = mean(liver_tmp[which(liver_tmp$element == "promoter"), "lib_fc"]),
        max_pro_fc = max(liver_tmp[which(liver_tmp$element == "promoter"), "lib_fc"]),
        median_pro_fc = median(liver_tmp[which(liver_tmp$element == "promoter"), "lib_fc"]),
        mean_enh_motif_num = mean(liver_tmp[which(liver_tmp$element == "enhancer"), "motif_num"]),
        max_enh_motif_num = max(liver_tmp[which(liver_tmp$element == "enhancer"), "motif_num"]),
        median_enh_motif_num = median(liver_tmp[which(liver_tmp$element == "enhancer"), "motif_num"]),
        mean_pro_motif_num = mean(liver_tmp[which(liver_tmp$element == "promoter"), "motif_num"]),
        max_pro_motif_num = max(liver_tmp[which(liver_tmp$element == "promoter"), "motif_num"]),
        median_pro_motif_num = median(liver_tmp[which(liver_tmp$element == "promoter"), "motif_num"]),
        mean_enh_phylop = mean(liver_tmp[which(liver_tmp$element == "enhancer"), "mean_phylop"]),
        mean_pro_phylop = mean(liver_tmp[which(liver_tmp$element == "promoter"), "mean_phylop"]),
        median_enh_phylop = median(liver_tmp[which(liver_tmp$element == "enhancer"), "mean_phylop"]),
        median_pro_phylop = median(liver_tmp[which(liver_tmp$element == "promoter"), "mean_phylop"])
    )
    liver_og_sum <- rbind(liver_og_sum, liver_tmp)
}
write.csv(brain_og_sum, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/", sp, "_Brain_per_gene.csv"), row.names = FALSE)
write.csv(kidney_og_sum, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/", sp, "_Kidney_per_gene.csv"), row.names = FALSE)
write.csv(liver_og_sum, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/", sp, "_Liver_per_gene.csv"), row.names = FALSE)
library(corrplot)
for (tis in c("Brain", "Kidney", "Liver")) {
    if (tis == "Brain") {
        sum <- brain_og_sum
    } else if (tis == "Kidney") {
        sum <- kidney_og_sum
    } else {
        sum <- liver_og_sum
    }
    sum$fpkm <- as.numeric(sum$fpkm)
    sum[sum == -Inf] <- NA
    rhos <- cor(sum[, -1], method = "spearman", use = "pairwise.complete.obs")
    p <- cor.mtest(sum[, -1], method = "spearman", use = "pairwise.complete.obs")$p
    # convert rhos matrix into long dataframe
    rhos[upper.tri(rhos)] <- NA
    rho_df <- reshape2::melt(rhos, na.rm = TRUE)
    rho_df <- rho_df[-which(rho_df$Var1 == rho_df$Var2), ]
    p[upper.tri(p)] <- NA
    p_df <- reshape2::melt(p, na.rm = TRUE) %>% .[-which(.$Var1 == .$Var2), ]
    res_df <- merge(rho_df, p_df, by = c("Var1", "Var2"))
    colnames(res_df)[3:4] <- c("rho", "p")
    e_ex <- sum(sum$enhancer_num > 0 & sum$fpkm > 0)
    e_nex <- sum(sum$enhancer_num > 0 & sum$fpkm == 0)
    ne_ex <- sum(sum$enhancer_num == 0 & sum$fpkm > 0)
    ne_nex <- sum(sum$enhancer_num == 0 & sum$fpkm == 0)
    e_ex_chiq <- chisq.test(matrix(c(e_ex, e_nex, ne_ex, ne_nex), nrow = 2))
    res_df <- rbind(res_df, data.frame(Var1 = "enhancer", Var2 = "express", rho = ifelse(e_ex > e_nex, 1, -1), p = e_ex_chiq$p.value))
    p_ex <- sum(sum$promoter_num > 0 & sum$fpkm > 0)
    p_nex <- sum(sum$promoter_num > 0 & sum$fpkm == 0)
    np_ex <- sum(sum$promoter_num == 0 & sum$fpkm > 0)
    np_nex <- sum(sum$promoter_num == 0 & sum$fpkm == 0)
    p_ex_chiq <- chisq.test(matrix(c(p_ex, p_nex, np_ex, np_nex), nrow = 2))
    res_df <- rbind(res_df, data.frame(Var1 = "promoter", Var2 = "express", rho = ifelse(p_ex > p_nex, 1, -1), p = p_ex_chiq$p.value))
    write.csv(res_df, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/chiq_in_sp/", sp, "_", tis, "_chiq.csv"), row.names = FALSE)
}

## use orthogroups in ortho_me_adj, and fpkm value in raw_FPKM_me_0_sep, becauce
## only the raw FPKM (or count) can tell if a gene is not expressed in a species
# for (sp in unique(raw_FPKM_me_0_sep$species)) {
#     ortho_tmp <- ortho_me_adj[which(ortho_me_adj[, sp] != "NULL"), c("id", sp)]
#     fpkm_tmp <-
#         raw_FPKM_me_0_sep[which((raw_FPKM_me_0_sep$species == sp) & (raw_FPKM_me_0_sep$OGID %in% ortho_tmp$id)), ]
#     bkl_tmp <- data.frame()
#     for (tis in c("Brain", "Kidney", "Liver")) {
#         # read the raw chipseeker annotation of promoter and enhancer
#         if (sp == "Macaca_fascicularis") {
#             enhancer_tmp <- read.csv(paste0(
#                 "/media/Data/zhangz/chip/analysis/Macaca_mulatta/compare2/Macaca_mulatta_",
#                 tis, "_enhancer_info.csv"
#             ), header = TRUE, row.names = NULL)
#             promoter_tmp <- read.csv(paste0(
#                 "/media/Data/zhangz/chip/analysis/Macaca_mulatta/compare2/Macaca_mulatta_",
#                 tis, "_promoter_info.csv"
#             ), header = TRUE, row.names = NULL)
#         } else {
#             enhancer_tmp <- read.csv(paste0(
#                 "/media/Data/zhangz/chip/analysis/", sp, "/compare2/",
#                 sp, "_", tis, "_enhancer_info.csv"
#             ), header = TRUE, row.names = NULL)
#             promoter_tmp <- read.csv(paste0(
#                 "/media/Data/zhangz/chip/analysis/", sp, "/compare2/",
#                 sp, "_", tis, "_promoter_info.csv"
#             ), header = TRUE, row.names = NULL)
#         }
#         tmp <- data.frame(
#             ORID = ortho_tmp$id, fpkm = "NULL",
#             promoter_num = 0, promoter_fc = 0, promoter_gc = 0, promoter_length = 0,
#             promoter_align = 0, promoter_con = 0, promoter_tss = 0, promoter_motif_num = 0, promoter_peak_id = "",
#             enhancer_num = 0, enhancer_fc = 0, enhancer_gc = 0, enhancer_length = 0,
#             enhancer_align = 0, enhancer_con = 0, enhancer_tss = 0, enhancer_motif_num = 0, enhancer_peak_id = ""
#         )
#         for (ogid in tmp$ORID) {
#             tmp$fpkm[which(tmp$ORID == ogid)] <-
#                 fpkm_tmp[which(fpkm_tmp$OGID == ogid), tis]
#             g_t <- ortho_tmp[which(ortho_tmp$id == ogid), sp]
#             gene <- str_split(g_t, "[|]")[[1]][1]
#             transcript <- str_split(g_t, "[|]")[[1]][2]
#             enhancer_g <- enhancer_tmp[which(enhancer_tmp$geneId == gene), ]
#             promoter_g <- promoter_tmp[which(promoter_tmp$geneId == gene), ]
#             tmp$promoter_num[which(tmp$ORID == ogid)] <-
#                 nrow(promoter_g)
#             tmp$enhancer_num[which(tmp$ORID == ogid)] <-
#                 nrow(enhancer_g)
#             tmp$promoter_fc[which(tmp$ORID == ogid)] <-
#                 mean(promoter_g$lib_fc)
#             tmp$promoter_gc[which(tmp$ORID == ogid)] <-
#                 mean(promoter_g$gc)
#             tmp$promoter_length[which(tmp$ORID == ogid)] <-
#                 sum(promoter_g$length)
#             tmp$promoter_align[which(tmp$ORID == ogid)] <-
#                 mean(promoter_g$align)
#             tmp$promoter_con[which(tmp$ORID == ogid)] <-
#                 mean(promoter_g$consrv)
#             tmp$promoter_tss[which(tmp$ORID == ogid)] <-
#                 mean(abs(promoter_g$distanceToTSS))
#             tmp$promoter_motif_num[which(tmp$ORID == ogid)] <-
#                 mean(promoter_g$motif_num)
#             tmp$promoter_peak_id[which(tmp$ORID == ogid)] <-
#                 paste(promoter_g$peak, collapse = ";")
#             tmp$enhancer_fc[which(tmp$ORID == ogid)] <-
#                 mean(enhancer_g$lib_fc)
#             tmp$enhancer_gc[which(tmp$ORID == ogid)] <-
#                 mean(enhancer_g$gc)
#             tmp$enhancer_length[which(tmp$ORID == ogid)] <-
#                 sum(enhancer_g$length)
#             tmp$enhancer_align[which(tmp$ORID == ogid)] <-
#                 mean(enhancer_g$align)
#             tmp$enhancer_con[which(tmp$ORID == ogid)] <-
#                 mean(enhancer_g$consrv)
#             tmp$enhancer_tss[which(tmp$ORID == ogid)] <-
#                 mean(abs(enhancer_g$distanceToTSS))
#             tmp$enhancer_motif_num[which(tmp$ORID == ogid)] <-
#                 mean(enhancer_g$motif_num)
#             tmp$enhancer_peak_id[which(tmp$ORID == ogid)] <-
#                 paste(enhancer_g$peak, collapse = ";")
#         }
#         ## spearman correlation, chi-square test
#         ## between fpkm and promoter/enhancer
#         tmp$fpkm <- as.numeric(tmp$fpkm)
#         bkl_tmp <- rbind(bkl_tmp, tmp)
#         enhancer_fc_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$enhancer_fc[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         enhancer_gc_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$enhancer_gc[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         enhancer_length_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$enhancer_length[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         enhancer_align_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$enhancer_align[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         enhancer_con_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$enhancer_con[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         enhancer_tss_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$enhancer_tss[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         enhancer_motif_num_spm <- cor.test(
#             x = tmp$fpkm[which((tmp$fpkm > 0) & (tmp$enhancer_motif_num > 0))],
#             y = tmp$enhancer_motif_num[which((tmp$fpkm > 0) & (tmp$enhancer_motif_num > 0))],
#             method = "spearman"
#         )
#         promoter_fc_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$promoter_fc[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         promoter_gc_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$promoter_gc[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         promoter_length_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$promoter_length[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         promoter_align_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$promoter_align[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         promoter_con_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$promoter_con[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         promoter_tss_spm <- cor.test(
#             x = tmp$fpkm[which(tmp$fpkm > 0)],
#             y = tmp$promoter_tss[which(tmp$fpkm > 0)],
#             method = "spearman"
#         )
#         promoter_motif_num_spm <- cor.test(
#             x = tmp$fpkm[which((tmp$fpkm > 0) & (tmp$promoter_motif_num > 0))],
#             y = tmp$promoter_motif_num[which((tmp$fpkm > 0) & (tmp$promoter_motif_num > 0))],
#             method = "spearman"
#         )
#         enhancer_num_spm <- cor.test(
#             x = tmp$fpkm[which((tmp$fpkm > 0) & (tmp$enhancer_num > 0))],
#             y = tmp$enhancer_num[which((tmp$fpkm > 0) & (tmp$enhancer_num > 0))],
#             method = "spearman"
#         )
#         promoter_num_spm <- cor.test(
#             x = tmp$fpkm[which((tmp$fpkm > 0) & (tmp$promoter_num > 0))],
#             y = tmp$promoter_num[which((tmp$fpkm > 0) & (tmp$promoter_num > 0))],
#             method = "spearman"
#         )
#         ## chi-square test
#         ## between fpkm and promoter, enhancer or both
#         p_ex <- nrow(tmp[which((tmp$promoter_num > 0) & (tmp$fpkm > 0)), ])
#         np_ex <- nrow(tmp[which((tmp$promoter_num == 0) & (tmp$fpkm > 0)), ])
#         p_nex <- nrow(tmp[which((tmp$promoter_num > 0) & (tmp$fpkm == 0)), ])
#         np_nex <- nrow(tmp[which((tmp$promoter_num == 0) & (tmp$fpkm == 0)), ])
#         p_ex_chiq <- chisq.test(
#             matrix(
#                 c(
#                     p_ex, np_ex,
#                     p_nex, np_nex
#                 ),
#                 nrow = 2
#             )
#         )
#         e_ex <- nrow(tmp[which((tmp$enhancer_num > 0) & (tmp$fpkm > 0)), ])
#         ne_ex <- nrow(tmp[which((tmp$enhancer_num == 0) & (tmp$fpkm > 0)), ])
#         e_nex <- nrow(tmp[which((tmp$enhancer_num > 0) & (tmp$fpkm == 0)), ])
#         ne_nex <- nrow(tmp[which((tmp$enhancer_num == 0) & (tmp$fpkm == 0)), ])
#         e_ex_chiq <- chisq.test(
#             matrix(
#                 c(
#                     e_ex, ne_ex,
#                     e_nex, ne_nex
#                 ),
#                 nrow = 2
#             )
#         )
#         ## write to st_chisq
#         st_chisq$promoter_num_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_num_spm$p.value
#         st_chisq$promoter_num_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_num_spm$estimate
#         st_chisq$enhancer_num_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_num_spm$p.value
#         st_chisq$enhancer_num_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_num_spm$estimate
#         st_chisq$p_ex[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- p_ex
#         st_chisq$np_ex[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- np_ex
#         st_chisq$p_nex[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- p_nex
#         st_chisq$np_nex[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- np_nex
#         st_chisq$p_ex_chiq_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- p_ex_chiq$p.value
#         if (p_ex > np_ex) {
#             st_chisq$p_ex_d[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- 1
#         }
#         st_chisq$e_ex[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- e_ex
#         st_chisq$ne_ex[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- ne_ex
#         st_chisq$e_nex[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- e_nex
#         st_chisq$ne_nex[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- ne_nex
#         st_chisq$e_ex_chiq_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- e_ex_chiq$p.value
#         if (e_ex > ne_ex) {
#             st_chisq$e_ex_d[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <- 1
#         }
#         st_chisq$enhancer_fc_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_fc_spm$p.value
#         st_chisq$enhancer_fc_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_fc_spm$estimate
#         st_chisq$enhancer_gc_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_gc_spm$p.value
#         st_chisq$enhancer_gc_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_gc_spm$estimate
#         st_chisq$enhancer_length_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_length_spm$p.value
#         st_chisq$enhancer_length_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_length_spm$estimate
#         st_chisq$enhancer_align_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_align_spm$p.value
#         st_chisq$enhancer_align_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_align_spm$estimate
#         st_chisq$enhancer_con_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_con_spm$p.value
#         st_chisq$enhancer_con_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_con_spm$estimate
#         st_chisq$enhancer_tss_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_tss_spm$p.value
#         st_chisq$enhancer_tss_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_tss_spm$estimate
#         st_chisq$enhancer_motif_num_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_motif_num_spm$p.value
#         st_chisq$enhancer_motif_num_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             enhancer_motif_num_spm$estimate
#         st_chisq$promoter_fc_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_fc_spm$p.value
#         st_chisq$promoter_fc_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_fc_spm$estimate
#         st_chisq$promoter_gc_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_gc_spm$p.value
#         st_chisq$promoter_gc_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_gc_spm$estimate
#         st_chisq$promoter_length_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_length_spm$p.value
#         st_chisq$promoter_length_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_length_spm$estimate
#         st_chisq$promoter_align_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_align_spm$p.value
#         st_chisq$promoter_align_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_align_spm$estimate
#         st_chisq$promoter_con_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_con_spm$p.value
#         st_chisq$promoter_con_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_con_spm$estimate
#         st_chisq$promoter_tss_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_tss_spm$p.value
#         st_chisq$promoter_tss_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_tss_spm$estimate
#         st_chisq$promoter_motif_num_spm_p[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_motif_num_spm$p.value
#         st_chisq$promoter_motif_num_spm_rho[which((st_chisq$species == sp) & (st_chisq$tissue == tis))] <-
#             promoter_motif_num_spm$estimate
#     }
#         ## write to sp_spm
#     sp_enhancer_num_spm <- cor.test(
#         x = enhancer_na_bkl_tmp$fpkm[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         y = enhancer_na_bkl_tmp$enhancer_num[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_promoter_num_spm <- cor.test(
#         x = promoter_na_bkl_tmp$fpkm[which(promoter_na_bkl_tmp$fpkm > 0)],
#         y = promoter_na_bkl_tmp$promoter_num[which(promoter_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_enhancer_fc_spm <- cor.test(
#         x = enhancer_na_bkl_tmp$fpkm[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         y = enhancer_na_bkl_tmp$enhancer_fc[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_enhancer_gc_spm <- cor.test(
#         x = enhancer_na_bkl_tmp$fpkm[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         y = enhancer_na_bkl_tmp$enhancer_gc[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_enhancer_length_spm <- cor.test(
#         x = enhancer_na_bkl_tmp$fpkm[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         y = enhancer_na_bkl_tmp$enhancer_length[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_enhancer_align_spm <- cor.test(
#         x = enhancer_na_bkl_tmp$fpkm[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         y = enhancer_na_bkl_tmp$enhancer_align[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_enhancer_con_spm <- cor.test(
#         x = enhancer_na_bkl_tmp$fpkm[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         y = enhancer_na_bkl_tmp$enhancer_con[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_enhancer_tss_spm <- cor.test(
#         x = enhancer_na_bkl_tmp$fpkm[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         y = enhancer_na_bkl_tmp$enhancer_tss[which(enhancer_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_enhancer_motif_num_spm <- cor.test(
#         x = enhancer_na_bkl_tmp$fpkm[which((enhancer_na_bkl_tmp$fpkm > 0) & (enhancer_na_bkl_tmp$enhancer_motif_num > 0))],
#         y = enhancer_na_bkl_tmp$enhancer_motif_num[which((enhancer_na_bkl_tmp$fpkm > 0) & (enhancer_na_bkl_tmp$enhancer_motif_num > 0))],
#         method = "spearman"
#     )
#     sp_promoter_fc_spm <- cor.test(
#         x = promoter_na_bkl_tmp$fpkm[which(promoter_na_bkl_tmp$fpkm > 0)],
#         y = promoter_na_bkl_tmp$promoter_fc[which(promoter_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_promoter_gc_spm <- cor.test(
#         x = promoter_na_bkl_tmp$fpkm[which(promoter_na_bkl_tmp$fpkm > 0)],
#         y = promoter_na_bkl_tmp$promoter_gc[which(promoter_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_promoter_length_spm <- cor.test(
#         x = promoter_na_bkl_tmp$fpkm[which(promoter_na_bkl_tmp$fpkm > 0)],
#         y = promoter_na_bkl_tmp$promoter_length[which(promoter_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_promoter_align_spm <- cor.test(
#         x = promoter_na_bkl_tmp$fpkm[which(promoter_na_bkl_tmp$fpkm > 0)],
#         y = promoter_na_bkl_tmp$promoter_align[which(promoter_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_promoter_con_spm <- cor.test(
#         x = promoter_na_bkl_tmp$fpkm[which(promoter_na_bkl_tmp$fpkm > 0)],
#         y = promoter_na_bkl_tmp$promoter_con[which(promoter_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_promoter_tss_spm <- cor.test(
#         x = promoter_na_bkl_tmp$fpkm[which(promoter_na_bkl_tmp$fpkm > 0)],
#         y = promoter_na_bkl_tmp$promoter_tss[which(promoter_na_bkl_tmp$fpkm > 0)],
#         method = "spearman"
#     )
#     sp_promoter_motif_num_spm <- cor.test(
#         x = promoter_na_bkl_tmp$fpkm[which((promoter_na_bkl_tmp$fpkm > 0) & (promoter_na_bkl_tmp$promoter_motif_num > 0))],
#         y = promoter_na_bkl_tmp$promoter_motif_num[which((promoter_na_bkl_tmp$fpkm > 0) & (promoter_na_bkl_tmp$promoter_motif_num > 0))],
#         method = "spearman"
#     )
#     sp_spm$enhancer_num_spm_p[which(sp_spm$species == sp)] <- sp_enhancer_num_spm$p.value
#     sp_spm$enhancer_num_spm_rho[which(sp_spm$species == sp)] <- sp_enhancer_num_spm$estimate
#     sp_spm$promoter_num_spm_p[which(sp_spm$species == sp)] <- sp_promoter_num_spm$p.value
#     sp_spm$promoter_num_spm_rho[which(sp_spm$species == sp)] <- sp_promoter_num_spm$estimate
#     sp_spm$enhancer_fc_spm_p[which(sp_spm$species == sp)] <- sp_enhancer_fc_spm$p.value
#     sp_spm$enhancer_fc_spm_rho[which(sp_spm$species == sp)] <- sp_enhancer_fc_spm$estimate
#     sp_spm$enhancer_gc_spm_p[which(sp_spm$species == sp)] <- sp_enhancer_gc_spm$p.value
#     sp_spm$enhancer_gc_spm_rho[which(sp_spm$species == sp)] <- sp_enhancer_gc_spm$estimate
#     sp_spm$enhancer_length_spm_p[which(sp_spm$species == sp)] <- sp_enhancer_length_spm$p.value
#     sp_spm$enhancer_length_spm_rho[which(sp_spm$species == sp)] <- sp_enhancer_length_spm$estimate
#     sp_spm$enhancer_align_spm_p[which(sp_spm$species == sp)] <- sp_enhancer_align_spm$p.value
#     sp_spm$enhancer_align_spm_rho[which(sp_spm$species == sp)] <- sp_enhancer_align_spm$estimate
#     sp_spm$enhancer_con_spm_p[which(sp_spm$species == sp)] <- sp_enhancer_con_spm$p.value
#     sp_spm$enhancer_con_spm_rho[which(sp_spm$species == sp)] <- sp_enhancer_con_spm$estimate
#     sp_spm$enhancer_tss_spm_p[which(sp_spm$species == sp)] <- sp_enhancer_tss_spm$p.value
#     sp_spm$enhancer_tss_spm_rho[which(sp_spm$species == sp)] <- sp_enhancer_tss_spm$estimate
#     sp_spm$enhancer_motif_num_spm_p[which(sp_spm$species == sp)] <- sp_enhancer_motif_num_spm$p.value
#     sp_spm$enhancer_motif_num_spm_rho[which(sp_spm$species == sp)] <- sp_enhancer_motif_num_spm$estimate
#     sp_spm$promoter_fc_spm_p[which(sp_spm$species == sp)] <- sp_promoter_fc_spm$p.value
#     sp_spm$promoter_fc_spm_rho[which(sp_spm$species == sp)] <- sp_promoter_fc_spm$estimate
#     sp_spm$promoter_gc_spm_p[which(sp_spm$species == sp)] <- sp_promoter_gc_spm$p.value
#     sp_spm$promoter_gc_spm_rho[which(sp_spm$species == sp)] <- sp_promoter_gc_spm$estimate
#     sp_spm$promoter_length_spm_p[which(sp_spm$species == sp)] <- sp_promoter_length_spm$p.value
#     sp_spm$promoter_length_spm_rho[which(sp_spm$species == sp)] <- sp_promoter_length_spm$estimate
#     sp_spm$promoter_align_spm_p[which(sp_spm$species == sp)] <- sp_promoter_align_spm$p.value
#     sp_spm$promoter_align_spm_rho[which(sp_spm$species == sp)] <- sp_promoter_align_spm$estimate
#     sp_spm$promoter_con_spm_p[which(sp_spm$species == sp)] <- sp_promoter_con_spm$p.value
#     sp_spm$promoter_con_spm_rho[which(sp_spm$species == sp)] <- sp_promoter_con_spm$estimate
#     sp_spm$promoter_tss_spm_p[which(sp_spm$species == sp)] <- sp_promoter_tss_spm$p.value
#     sp_spm$promoter_tss_spm_rho[which(sp_spm$species == sp)] <- sp_promoter_tss_spm$estimate
#     sp_spm$promoter_motif_num_spm_p[which(sp_spm$species == sp)] <- sp_promoter_motif_num_spm$p.value
#     sp_spm$promoter_motif_num_spm_rho[which(sp_spm$species == sp)] <- sp_promoter_motif_num_spm$estimate
#     ## chi-square test
#     sp_p_ex <- nrow(bkl_tmp[which((bkl_tmp$promoter_num > 0) & (bkl_tmp$fpkm > 0)), ])
#     sp_np_ex <- nrow(bkl_tmp[which((bkl_tmp$promoter_num == 0) & (bkl_tmp$fpkm > 0)), ])
#     sp_p_nex <- nrow(bkl_tmp[which((bkl_tmp$promoter_num > 0) & (bkl_tmp$fpkm == 0)), ])
#     sp_np_nex <- nrow(bkl_tmp[which((bkl_tmp$promoter_num == 0) & (bkl_tmp$fpkm == 0)), ])
#     sp_p_ex_chiq <- chisq.test(
#         matrix(
#             c(
#                 sp_p_ex, sp_np_ex,
#                 sp_p_nex, sp_np_nex
#             ),
#             nrow = 2
#         )
#     )
#     sp_e_ex <- nrow(bkl_tmp[which((bkl_tmp$enhancer_num > 0) & (bkl_tmp$fpkm > 0)), ])
#     sp_ne_ex <- nrow(bkl_tmp[which((bkl_tmp$enhancer_num == 0) & (bkl_tmp$fpkm > 0)), ])
#     sp_e_nex <- nrow(bkl_tmp[which((bkl_tmp$enhancer_num > 0) & (bkl_tmp$fpkm == 0)), ])
#     sp_ne_nex <- nrow(bkl_tmp[which((bkl_tmp$enhancer_num == 0) & (bkl_tmp$fpkm == 0)), ])
#     sp_e_ex_chiq <- chisq.test(
#         matrix(
#             c(
#                 sp_e_ex, sp_ne_ex,
#                 sp_e_nex, sp_ne_nex
#             ),
#             nrow = 2
#         )
#     )
#     sp_spm$p_ex[which(sp_spm$species == sp)] <- sp_p_ex
#     sp_spm$np_ex[which(sp_spm$species == sp)] <- sp_np_ex
#     sp_spm$p_nex[which(sp_spm$species == sp)] <- sp_p_nex
#     sp_spm$np_nex[which(sp_spm$species == sp)] <- sp_np_nex
#     sp_spm$p_ex_chiq_p[which(sp_spm$species == sp)] <- sp_p_ex_chiq$p.value
#     if (sp_p_ex > sp_np_ex) {
#         sp_spm$p_ex_d[which(sp_spm$species == sp)] <- 1
#     }
#     sp_spm$e_ex[which(sp_spm$species == sp)] <- sp_e_ex
#     sp_spm$ne_ex[which(sp_spm$species == sp)] <- sp_ne_ex
#     sp_spm$e_nex[which(sp_spm$species == sp)] <- sp_e_nex
#     sp_spm$ne_nex[which(sp_spm$species == sp)] <- sp_ne_nex
#     sp_spm$e_ex_chiq_p[which(sp_spm$species == sp)] <- sp_e_ex_chiq$p.value
#     if (sp_e_ex > sp_ne_ex) {
#         sp_spm$e_ex_d[which(sp_spm$species == sp)] <- 1
#     }
# }
# write.table(st_chisq, "/media/Data/zhangz/chip/analysis/expression_mean/zhangz/st_chisq.tsv",
#     sep = "\t", quote = FALSE, row.names = FALSE
# )
# write.table(sp_spm, "/media/Data/zhangz/chip/analysis/expression_mean/zhangz/sp_spm.tsv",
#     sep = "\t", quote = FALSE, row.names = FALSE
# )
