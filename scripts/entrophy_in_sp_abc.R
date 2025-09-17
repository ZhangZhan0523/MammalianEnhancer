# load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/ortho_me_adj.RData")
library(optparse)
option_list <- list(
    make_option(
        c("-s", "--species"),
        type = "character", default = "/media/Data/zhangz/chip/scripts/info/info_using.csv",
        help = "species need to summarize"
    )
    # make_option(
    #     c("-d", "--distance"),
    #     type = "character", default = "100k",
    #     help = "annotation distance threshold"
    # )
)
opt_parser <- parse_args(OptionParser(option_list = option_list, usage = "%prog [options]"))
sp <- opt_parser$species
# dis <- opt_parser$distance
print(sp)
# print(dis)
library(magrittr)

dis <- "abc"
d_cnames <- list(
    "1M" = "anno_gene_1M", "500k" = "anno_gene_500k", "200k" = "anno_gene_200k",
    "100k" = "anno_gene_100k", "50k" = "anno_gene_50k", "abc" = "overlap_gene"
)
cname <- d_cnames[[dis]]
# sp_all <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/sp_all/", sp, "/", sp, "_all_1M.csv"),
#     header = TRUE
# )
sp_all <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE)
load("/media/Data/zhangz/chip/analysis/summary2/sum_all/qn_FPKM_0.RData")
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/fpkm_ortho.RData")
## contains ortho_me_adj, raw_FPKM_me_0_sep_pattern_fl
load(file = "/media/Data/zhangz/chip/analysis/expression_mean/zhangz/anno_gene_base2.RData")
raw_FPKM_me_0_sep_pattern_fl$species <- gsub("Macaca_fascicularis", "Macaca_mulatta", raw_FPKM_me_0_sep_pattern_fl$species)
gene_pattern <- raw_FPKM_me_0_sep_pattern_fl[which(raw_FPKM_me_0_sep_pattern_fl$species == sp), ]
gene_pattern$tau[is.na(gene_pattern$tau)] <- 0
colnames(qn_FPKM_0) <- gsub("Macaca_fascicularis", "Macaca_mulatta", colnames(qn_FPKM_0))
ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)

library(doMC)
doMC::registerDoMC(cores = 6)
library(plyr)
library(corrplot)
species_nots <- c(
    "Neophocaena_asiaeorientalis", "Rhinopithecus_roxellana",
    "Rhinolophus_ferrumequinum", "Myotis_ricketti"
)

lib_factor <- read.csv("/media/Data/zhangz/chip/scripts2/info/fc_library_factor.csv", header = TRUE)
# lib_factor$rep_factor <- 49
# lib_factor$rep_factor[lib_factor$species == "Rhinolophus_ferrumequinum" & lib_factor$tissue == "Liver"] <- 98
# lib_factor$rep_factor[lib_factor$species == "Myotis_chinensis" & lib_factor$tissue == "Brain"] <- 98
# write.csv(lib_factor, "/media/Data/zhangz/chip/scripts2/info/fc_library_factor.csv", row.names = FALSE)
raw_FPKM_me_0_sep[raw_FPKM_me_0_sep$species == "Macaca_fascicularis", "species"] <- "Macaca_mulatta"
colnames(ortho_me_adj)[which(colnames(ortho_me_adj) == "Macaca_fascicularis")] <- "Macaca_mulatta"
raw_FPKM_sp <- raw_FPKM_me_0_sep[which(raw_FPKM_me_0_sep$species == sp & !is.na(raw_FPKM_me_0_sep$Brain)), ]
if (sp %in% species_nots) {
    sp_ortho <- ortho_all[ortho_all[sp] != "NULL" & ortho_all[sp] != "", c("Orthogroup", sp)]
    colnames(sp_ortho)[1] <- "id"
} else {
    sp_ortho <- ortho_me_adj[ortho_me_adj[sp] != "NULL", c("id", sp)]
}
# brain_og_sum <- data.frame()
# kidney_og_sum <- data.frame()
# liver_og_sum <- data.frame()

library(stringr)
library(magrittr)
## find all genes annotated in sp_all
if (sp %in% species_nots) {
    gene_og_raw <- sp_all[, c(cname, "overlap_og")]
    gene_og_raw <- gene_og_raw[sp_all$overlap_og != "NULL" & sp_all$overlap_og != "" & sp_all$overlap_og != "NA", ]
    gene_og <- data.frame(
        gene = str_split(paste(gene_og_raw$overlap_gene, collapse = ";"), ";")[[1]],
        og = str_split(paste(gene_og_raw$overlap_og, collapse = ";"), ";")[[1]]
    )
    gene_og <- gene_og[-which(gene_og$og == "character(0)"), ]
    gene_og <- gene_og[!duplicated(gene_og$og), ]
} else {
    genes <- unique(str_split(paste(sp_all[[cname]], collapse = ";"), ";")[[1]]) %>%
        .[. != "NULL" & . != "" & . != "NA"]
    ogs <- unique(str_split(paste(sp_all[["overlap_og"]], collapse = ";"), ";")[[1]]) %>%
        .[. != "NULL" & . != "" & . != "NA"]
    gene_og <- data.frame(gene = genes, og = ogs)
}

calc_ent <- function(g_og) {
    #  r$> head(raw_FPKM_sp)
    #           OGID            species  Brain Kidney   Liver
    # 331  OG0000015 Petaurus_breviceps 3.9329 1.4476  0.6649
    # 1402 OG0000066 Petaurus_breviceps 0.0221 0.0538  0.0790
    # 2410 OG0000114 Petaurus_breviceps 0.0366 0.0297  0.0000
    # 2473 OG0000117 Petaurus_breviceps 0.0360 0.0219  0.0000
    # 2809 OG0000133 Petaurus_breviceps 1.4442 0.3422  5.3036
    # 2893 OG0000137 Petaurus_breviceps 8.8545 2.3820 13.6706
    mean_fpkm <- mean(as.numeric(raw_FPKM_sp[which(raw_FPKM_sp$OGID == g_og[1, 2]), c("Brain", "Kidney", "Liver")]))
    # mean expression of this species and orthologous group among three tissues
    sd_fpkm <- sd(as.numeric(raw_FPKM_sp[which(raw_FPKM_sp$OGID == g_og[1, 2]), c("Brain", "Kidney", "Liver")]))
    qn_FPKM_tmp <- qn_FPKM_0[qn_FPKM_0$OGID == g_og[1, 2], ]
    #        OGID Num Atelerix_albiventris.Brain.0 Atelerix_albiventris.Kidney.0
    mean_qn_fpkm_sps <- ifelse(!sp %in% species_nots, mean(as.numeric(qn_FPKM_tmp[, as.character(interaction(sp, unique(sp_all$tissue), 0))]), na.rm = TRUE), NA)
    # mean expression of this species and orthologous group among three tissues after quantile normalization
    sd_qn_fpkm_sps <- ifelse(!sp %in% species_nots, sd(as.numeric(qn_FPKM_tmp[, as.character(interaction(sp, unique(sp_all$tissue), 0))]), na.rm = TRUE), NA)
    mean_adj_fpkm <- mean(as.numeric(adj_FPKM_me[which(adj_FPKM_me$OGID == g_og[1, 2]), "fpkm"]))
    # mean expression of this species and orthologous group after removing batch effect among three tissues and replicates
    sd_adj_fpkm <- sd(as.numeric(adj_FPKM_me[which(adj_FPKM_me$OGID == g_og[1, 2]), "fpkm"]))
    cv_adj_fpkm <- ifelse(mean_adj_fpkm == 0, 0, sd_adj_fpkm / mean_adj_fpkm * 100)
    tau <- ifelse(length(gene_pattern$tau[gene_pattern$OGID == g_og[1, 2]]) > 0, as.numeric(gene_pattern$tau[gene_pattern$OGID == g_og[1, 2]]), 0)
    tau <- ifelse(is.na(tau), 0, tau)
    tmp <- sp_all[grep(paste0("(^|;)", g_og[1, 2], "(;|$)"), sp_all[["overlap_og"]]), c("peak", "ele_pattern", "mean_fc", "motif_num", "mean_phylop", "anno_distance", "element", "tissue", "lib_fc", "id", cname)]
    if (nrow(tmp) == 0) {
        print(g_og[1, 1])
        ent_row <- data.frame(
            OGID = g_og[1, 2], gene = g_og[1, 1],
            fpkm = mean_fpkm,
            fpkm_sd = sd_fpkm,
            fpkm_cv = ifelse(mean_fpkm == 0, 0, sd_fpkm / mean_fpkm * 100),
            adj_fpkm = mean_adj_fpkm,
            adj_fpkm_sd = sd_adj_fpkm,
            adj_fpkm_cv = cv_adj_fpkm,
            qn_fpkm_sps = mean_qn_fpkm_sps,
            qn_fpkm_sd = sd_qn_fpkm_sps,
            qn_fpkm_cv = ifelse(mean_qn_fpkm_sps == 0, 0, sd_qn_fpkm_sps / mean_qn_fpkm_sps * 100),
            tau = tau,
            ent1 = 0, ent2 = 0, unique_ele_num = 0
        )
    } else {
        for (i in seq_len(nrow(tmp))) {
            tmp$T[i] <- nchar(tmp$ele_pattern[i])
            tmp$raw_p[i] <- tmp$lib_fc[i] / lib_factor[lib_factor$species == sp & lib_factor$tissue == tmp$tissue[i], "rep_factor"]
        }
        tmp <- tmp[order(tmp$raw_p, decreasing = TRUE), ]
        # state_num <- 2^length(unique(tmp$id))
        # make a table of different state, each row correspond to a switch with two
        # possible state: 1 and 0, each column correspond to a whole state of all switches
        # make sure that each column is different from each other
        # 假设我们有3个开关
        # num_switches <- 3

        # # 创建一个列表，列表中的每个元素都是c(0, 1)，表示开关的两种状态
        # list_of_states <- rep(list(c(0, 1)), num_switches)

        # # 使用expand.grid生成所有可能的组合
        # df <- expand.grid(list_of_states)

        # # 打印结果
        # print(df)
        ele_num <- length(unique(tmp$id))
        unique_num <- ele_num
        ent2 <- nrow(tmp)
        # print(ele_num)
        if (ele_num > 12) {
            tmp <- tmp[tmp$id %in% unique(tmp$id)[1:12], ]
        }
        ele_num <- length(unique(tmp$id))
        list_of_states <- rep(list(c(0, 1)), ele_num)
        state_df <- expand.grid(list_of_states)
        p_df <- data.frame(id = unique(tmp$id), p1 = 0, p0 = 0)
        for (i in seq_len(nrow(p_df))) {
            p_df$p1[i] <- sum(tmp$raw_p[tmp$id == p_df$id[i]]) / 3
        }
        p_df$p0 <- 1 - p_df$p1
        colnames(state_df) <- unique(tmp$id)
        state_df$p <- 0
        for (i in seq_len(nrow(state_df))) {
            p_tmp <- 1
            for (j in 1:ele_num) {
                p_tmp <- p_tmp * p_df[p_df$id == colnames(state_df)[j], paste0("p", state_df[i, j])]
            }
            state_df$p[i] <- p_tmp
        }
        ent1 <- sum(-state_df$p * log2(state_df$p))
        # ent2 <- nrow(tmp)
        ent_row <- data.frame(
            OGID = g_og[1, 2],
            gene = g_og[1, 1],
            fpkm = mean_fpkm,
            fpkm_sd = sd_fpkm,
            fpkm_cv = ifelse(mean_fpkm == 0, 0, sd_fpkm / mean_fpkm),
            adj_fpkm = mean_adj_fpkm,
            adj_fpkm_sd = sd_adj_fpkm,
            adj_fpkm_cv = cv_adj_fpkm,
            qn_fpkm_sps = mean_qn_fpkm_sps,
            qn_fpkm_sd = sd_qn_fpkm_sps,
            qn_fpkm_cv = ifelse(mean_qn_fpkm_sps == 0, 0, sd_qn_fpkm_sps / mean_qn_fpkm_sps),
            tau = tau,
            ent1 = ent1, ent2 = ent2, unique_ele_num = unique_num
        )
    }
    return(ent_row)
}
# ent <- rbind(ent, data.frame(
#     OGID = sp_ortho$id[!sp_ortho$id %in% ent$OGID],
#     gene = str_split(sp_ortho[[sp]][!sp_ortho$id %in% ent$OGID], "[|]")[[1]][1],
#     fpkm = 0,
#     fpkm_sd = 0,
#     fpkm_se = 0,
#     tau = 0,
#     ent1 = 0, ent2 = 0, unique_ele_num = 0
# ))
ent <- plyr::adply(gene_og, 1, calc_ent, .parallel = TRUE)
ent <- ent[, -which(colnames(ent) == "og")]
# for (og in sp_ortho$id[!sp_ortho$id %in% ent$OGID])
calc_ent0 <- function(og) {
    og <- og[[1]]
    mean_fpkm <- mean(as.numeric(raw_FPKM_sp[which(raw_FPKM_sp$OGID == og), c("Brain", "Kidney", "Liver")]))
    sd_fpkm <- sd(as.numeric(raw_FPKM_sp[which(raw_FPKM_sp$OGID == og), c("Brain", "Kidney", "Liver")]))
    mean_adj_fpkm <- mean(as.numeric(adj_FPKM_me[which(adj_FPKM_me$OGID == og), "fpkm"]))
    sd_adj_fpkm <- sd(as.numeric(adj_FPKM_me[which(adj_FPKM_me$OGID == og), "fpkm"]))
    cv_adj_fpkm <- ifelse(mean_adj_fpkm == 0, 0, sd_adj_fpkm / mean_adj_fpkm * 100)
    qn_FPKM_tmp <- qn_FPKM_0[qn_FPKM_0$OGID == og, ]
    mean_qn_fpkm_sps <- ifelse(!sp %in% species_nots, mean(as.numeric(qn_FPKM_tmp[, as.character(interaction(sp, unique(sp_all$tissue), 0))]), na.rm = TRUE), NA)
    sd_qn_fpkm_sps <- ifelse(!sp %in% species_nots, sd(as.numeric(qn_FPKM_tmp[, as.character(interaction(sp, unique(sp_all$tissue), 0))]), na.rm = TRUE), NA)
    tau <- ifelse(length(gene_pattern$tau[gene_pattern$OGID == og]) > 0, as.numeric(gene_pattern$tau[gene_pattern$OGID == og]), 0)
    tau <- ifelse(is.na(tau), 0, tau)
    ent_row <- data.frame(
        OGID = og,
        gene = str_split(sp_ortho[[sp]][sp_ortho$id == og], "[|]")[[1]][1],
        fpkm = mean_fpkm,
        fpkm_sd = sd_fpkm,
        fpkm_cv = ifelse(mean_fpkm == 0, 0, sd_fpkm / mean_fpkm),
        adj_fpkm = mean_adj_fpkm,
        adj_fpkm_sd = sd_adj_fpkm,
        adj_fpkm_cv = cv_adj_fpkm,
        qn_fpkm_sps = mean_qn_fpkm_sps,
        qn_fpkm_sd = sd_qn_fpkm_sps,
        qn_fpkm_cv = ifelse(mean_qn_fpkm_sps == 0, 0, sd_qn_fpkm_sps / mean_qn_fpkm_sps),
        tau = tau,
        ent1 = 0, ent2 = 0, unique_ele_num = 0
    )
    return(ent_row)
}
if (!sp %in% species_nots) {
    ent0 <- plyr::ldply(as.list(sp_ortho$id[!sp_ortho$id %in% ent$OGID]), calc_ent0, .parallel = TRUE)
    ent <- rbind(ent, ent0)
}
write.csv(ent, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), row.names = FALSE)
## do correlation test
rhos <- cor(ent[, -c(1, 2)], method = "spearman", use = "pairwise.complete.obs")
p <- cor.mtest(ent[, -c(1, 2)], method = "spearman", use = "pairwise.complete.obs")$p
# convert rhos matrix into long dataframe
rhos[upper.tri(rhos)] <- NA
rho_df <- reshape2::melt(rhos, na.rm = TRUE)
rho_df <- rho_df[-which(rho_df$Var1 == rho_df$Var2), ]
p[upper.tri(p)] <- NA
p_df <- reshape2::melt(p, na.rm = TRUE) %>% .[-which(.$Var1 == .$Var2), ]
res_df <- merge(rho_df, p_df, by = c("Var1", "Var2"))
colnames(res_df)[3:4] <- c("rho", "p")
write.csv(res_df, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent_cor.csv"), row.names = FALSE)

ent2 <- ent[ent$ent1 != 0 & ent$ent2 != 0, ]
colnames(ent2)[3:ncol(ent2)] <- paste0("no0_", colnames(ent2)[3:ncol(ent2)])
rhos2 <- cor(ent2[, -c(1, 2)], method = "spearman", use = "pairwise.complete.obs")
p2 <- cor.mtest(ent2[, -c(1, 2)], method = "spearman", use = "pairwise.complete.obs")$p
# convert rhos matrix into long dataframe
rhos2[upper.tri(rhos2)] <- NA
rho_df2 <- reshape2::melt(rhos2, na.rm = TRUE)
rho_df2 <- rho_df2[-which(rho_df2$Var1 == rho_df2$Var2), ]
p2[upper.tri(p2)] <- NA
p_df2 <- reshape2::melt(p2, na.rm = TRUE) %>% .[-which(.$Var1 == .$Var2), ]
res_df2 <- merge(rho_df2, p_df2, by = c("Var1", "Var2"))
colnames(res_df2)[3:4] <- c("rho", "p")
write.csv(res_df2, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent_cor_no0_adj.csv"), row.names = FALSE)
