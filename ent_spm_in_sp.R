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
library(corrplot)
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/fpkm_ortho.RData")

ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/", sp, "_ent.csv"), header = TRUE)
colnames(ent)[4] <- "fpkm_cv"
for (og in ent$OGID) {
    mean_adj_fpkm <- mean(as.numeric(adj_FPKM_me[which(adj_FPKM_me$OGID == og), "fpkm"]))
    sd_adj_fpkm <- sd(as.numeric(adj_FPKM_me[which(adj_FPKM_me$OGID == og), "fpkm"]))
    cv_adj_fpkm <- sd_adj_fpkm / mean_adj_fpkm * 100
    ent[ent$OGID == og, "adj_fpkm"] <- mean_adj_fpkm
    ent[ent$OGID == og, "adj_fpkm_sd"] <- sd_adj_fpkm
    ent[ent$OGID == og, "adj_fpkm_cv"] <- cv_adj_fpkm
}

## do correlation test

rhos <- cor(ent[, -1], method = "spearman", use = "pairwise.complete.obs")
p <- cor.mtest(ent[, -1], method = "spearman", use = "pairwise.complete.obs")$p
# convert rhos matrix into long dataframe
rhos[upper.tri(rhos)] <- NA
rho_df <- reshape2::melt(rhos, na.rm = TRUE)
rho_df <- rho_df[-which(rho_df$Var1 == rho_df$Var2), ]
p[upper.tri(p)] <- NA
p_df <- reshape2::melt(p, na.rm = TRUE) %>% .[-which(.$Var1 == .$Var2), ]
res_df <- merge(rho_df, p_df, by = c("Var1", "Var2"))
colnames(res_df)[3:4] <- c("rho", "p")
write.csv(res_df, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/", sp, "_ent_cor_adj.csv"), row.names = FALSE)

ent2 <- ent[ent$ent1 != 0 & ent$ent2 != 0, ]
colnames(ent2)[2:ncol(ent2)] <- paste0("no0_", colnames(ent2)[2:ncol(ent2)])
rhos2 <- cor(ent2[, -1], method = "spearman", use = "pairwise.complete.obs")
p2 <- cor.mtest(ent2[, -1], method = "spearman", use = "pairwise.complete.obs")$p
# convert rhos matrix into long dataframe
rhos2[upper.tri(rhos2)] <- NA
rho_df2 <- reshape2::melt(rhos2, na.rm = TRUE)
rho_df2 <- rho_df2[-which(rho_df2$Var1 == rho_df2$Var2), ]
p2[upper.tri(p2)] <- NA
p_df2 <- reshape2::melt(p2, na.rm = TRUE) %>% .[-which(.$Var1 == .$Var2), ]
res_df2 <- merge(rho_df2, p_df2, by = c("Var1", "Var2"))
colnames(res_df2)[3:4] <- c("rho", "p")
write.csv(res_df2, paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp/", sp, "_ent_cor_no0_adj.csv"), row.names = FALSE)
