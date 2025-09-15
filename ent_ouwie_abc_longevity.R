# library(optparse)
# option_list <- list(
#     make_option(
#         c("-d", "--distance"),
#         type = "character", default = "100k",
#         help = "annotation distance threshold"
#     )
# )
# opt_parser <- parse_args(OptionParser(option_list = option_list, usage = "%prog [options]"))
# # sp <- opt_parser$species
# dis <- opt_parser$distance
# # print(sp)
# print(dis)
library(magrittr)
library(phytools)
library(OUwie)
library(tidyr)
library(doMC)
library(plyr)
doMC::registerDoMC(50)

setwd("/media/Data/zhangz/chip/analysis/summary2/ent_ouwie/longevity/")
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
species_ll <- c(
    "Sus_scrofa", "Neophocaena_asiaeorientalis", "Lama_glama", "Equus_asinus",
    "Equus_caballus", "Felis_catus", "Rhinolophus_ferrumequinum",
    "Myotis_chinensis", "Rhinopithecus_roxellana", "Macaca_mulatta"
)
species_ll_30 <- c(
    "Neophocaena_asiaeorientalis", "Equus_asinus", "Felis_catus", "Rhinolophus_ferrumequinum",
    "Equus_caballus", "Macaca_mulatta"
)
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/fpkm_ortho.RData")
ogid20 <- read.csv("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/ortho_me_20_id.csv", header = TRUE, row.names = NULL)
colnames(ortho_me_20_adj)[colnames(ortho_me_20_adj) == "Macaca_fascicularis"] <- "Macaca_mulatta"

# load and concat ent data
# ent1_all <- data.frame(OGID = ogid20$x)
# ent2_all <- data.frame(OGID = ogid20$x)
# for (sp in species_ts) {
#     ent <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/sum_all/ent_in_sp_abc/", sp, "_", "abc_1M_ent.csv"), header = TRUE)
#     ent1_all[sp] <- 0
#     ent2_all[sp] <- 0
#     for (i in seq_len(nrow(ent1_all))) {
#         if (ortho_me_20_adj[[sp]][ortho_me_20_adj$id == ent1_all$OGID[i]] == "NULL") {
#             ent1_all[[sp]][i] <- NA
#             ent2_all[[sp]][i] <- NA
#         } else {
#             ent1_all[[sp]][i] <- ent[ent$OGID == ent1_all$OGID[i], "ent1"]
#             ent2_all[[sp]][i] <- ent[ent$OGID == ent2_all$OGID[i], "ent2"]
#         }
#     }
# }
# write.csv(ent1_all, "ent1_all.csv", row.names = FALSE)
# write.csv(ent2_all, "ent2_all.csv", row.names = FALSE)
ent1_all <- read.csv("ent1_all.csv", header = TRUE)
ent2_all <- read.csv("ent2_all.csv", header = TRUE)
# check phylogenetic signal of each OGID
sp_tre <- ape::read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
# ent_phylosig <- data.frame(OGID = ogid20$x)
calc_physig <- function(i, ent1_all = ent2_all, ent2_all = ent2_all) {
    ent1 <- as.numeric(ent1_all[i, 2:22])
    names(ent1) <- colnames(ent1_all)[2:22]
    ent2 <- as.numeric(ent2_all[i, 2:22])
    names(ent2) <- colnames(ent2_all)[2:22]
    if (length(which(ent1 != 0 & !is.na(ent1))) != 0) {
        ent1_k <- phylosig(sp_tre, ent1, method = "K", test = TRUE)
        ent1_lambda <- phylosig(sp_tre, ent1, method = "lambda", test = TRUE)
        ent1_K <- ent1_k$K
        ent1_K_P <- ent1_k$P
        ent1_lambda_l <- ent1_lambda$lambda
        ent1_lambda_P <- ent1_lambda$P
    } else {
        ent1_K <- NA
        ent1_K_P <- NA
        ent1_lambda_l <- NA
        ent1_lambda_P <- NA
    }
    if (length(which(ent2 != 0 & !is.na(ent2))) != 0) {
        ent2_k <- phylosig(sp_tre, ent2, method = "K", test = TRUE)
        ent2_lambda <- phylosig(sp_tre, ent2, method = "lambda", test = TRUE)
        ent2_K <- ent2_k$K
        ent2_K_P <- ent2_k$P
        ent2_lambda_l <- ent2_lambda$lambda
        ent2_lambda_P <- ent2_lambda$P
    } else {
        ent2_K <- NA
        ent2_K_P <- NA
        ent2_lambda_l <- NA
        ent2_lambda_P <- NA
    }
    return(data.frame(
        OGID = ent1_all$OGID[i],
        ent1_K = ent1_K, ent1_K_P = ent1_K_P, ent1_lambda = ent1_lambda_l, ent1_lambda_P = ent1_lambda_P,
        ent2_K = ent2_K, ent2_K_P = ent2_K_P, ent2_lambda = ent2_lambda_l, ent2_lambda_P = ent2_lambda_P
    ))
}
# ent_phylosig <- plyr::adply(seq_len(nrow(ent1_all)), 1, calc_physig, .parallel = TRUE)

# write.csv(ent_phylosig, "ent_phylosig_abc.csv", row.names = FALSE)
ent_phylosig <- read.csv("ent_phylosig_abc.csv", header = TRUE)

# ## remove ent = 0
# ent_phylosig_no0 <- data.frame(OGID = ogid20$x)

# ent1_all$zcount <- rowSums(ent1_all[, 2:22] != 0 & !is.na(ent1_all[, 2:22]), na.rm = TRUE)
# ent2_all$zcount <- rowSums(ent2_all[, 2:22] != 0 & !is.na(ent2_all[, 2:22]), na.rm = TRUE)
# ent1_no0 <- ent1_all[ent1_all$zcount > 10, ]
# ent2_no0 <- ent2_all[ent2_all$zcount > 10, ]
# ent_phylosig_no0 <- ent_phylosig[ent_phylosig$OGID %in% ent1_no0$OGID, ]
# write.csv(ent1_no0, "ent1_no0.csv", row.names = FALSE)
# write.csv(ent2_no0, "ent2_no0.csv", row.names = FALSE)
# write.csv(ent_phylosig_no0, "ent_phylosig_no0.csv", row.names = FALSE)

# sp_tre <- ape::read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
# ent_phylosig <- data.frame(OGID = ogid20$x)
# for (i in seq_len(nrow(ent1_all))) {
#     ent1 <- as.numeric(ent1_all[i, 2:22])
#     names(ent1) <- colnames(ent1_all)[2:22]
#     ent2 <- as.numeric(ent2_all[i, 2:22])
#     names(ent2) <- colnames(ent2_all)[2:22]
#     if (length(which(ent1 != 0 & !is.na(ent1))) != 0) {
#         ent1_k <- phylosig(sp_tre, ent1, method = "K", test = TRUE)
#         ent1_lambda <- phylosig(sp_tre, ent1, method = "lambda", test = TRUE)
#         ent_phylosig[i, "ent1_K"] <- ent1_k$K
#         ent_phylosig[i, "ent1_K_P"] <- ent1_k$P
#         ent_phylosig[i, "ent1_lambda"] <- ent1_lambda$lambda
#         ent_phylosig[i, "ent1_lambda_P"] <- ent1_lambda$P
#     } else {
#         ent_phylosig[i, "ent1_K"] <- NA
#         ent_phylosig[i, "ent1_K_P"] <- NA
#         ent_phylosig[i, "ent1_lambda"] <- NA
#         ent_phylosig[i, "ent1_lambda_P"] <- NA
#     }
#     if (length(which(ent2 != 0 & !is.na(ent2))) != 0) {
#         ent2_k <- phylosig(sp_tre, ent2, method = "K", test = TRUE)
#         ent2_lambda <- phylosig(sp_tre, ent2, method = "lambda", test = TRUE)
#         ent_phylosig[i, "ent2_K"] <- ent2_k$K
#         ent_phylosig[i, "ent2_K_P"] <- ent2_k$P
#         ent_phylosig[i, "ent2_lambda"] <- ent2_lambda$lambda
#         ent_phylosig[i, "ent2_lambda_P"] <- ent2_lambda$P
#     } else {
#         ent_phylosig[i, "ent2_K"] <- NA
#         ent_phylosig[i, "ent2_K_P"] <- NA
#         ent_phylosig[i, "ent2_lambda"] <- NA
#         ent_phylosig[i, "ent2_lambda_P"] <- NA
#     }
# }
# write.csv(ent_phylosig, paste0("ent_phylosig_abc.csv"), row.names = FALSE)
# ent_phylosig <- read.csv(paste0(dis, "/ent_phylosig_", dis, ".csv"), header = TRUE)
test_og <- ent_phylosig[ent_phylosig$ent1_K_P < 0.05 & ent_phylosig$ent1_K >= 1, ] %>% drop_na()
print(paste0("phylosig og num is : ", nrow(test_og)))
# head(ent_phylosig[is.na(ent_phylosig$ent1_K), ])
# ent_num <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ent_evemodel/ent_num.csv", header = TRUE)
# sp_tre <- ape::read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
sp_ts_tre <- keep.tip(sp_tre, species_ts)

myOUwie <- function(tree, trait, model, ...) {
    not_na <- !is.na(trait$ent1)
    tmp_tree <- keep.tip(tree, as.vector(trait$species[not_na]))
    tmp_trait <- trait[not_na, ]
    tmp_res <- try(OUwie(tmp_tree, tmp_trait, model = model, ...))
    return(tmp_res)
}

ent1_melt <- reshape2::melt(ent1_all, id.vars = "OGID", variable.name = "species", value.name = "ent1")
ent2_melt <- reshape2::melt(ent2_all, id.vars = "OGID", variable.name = "species", value.name = "ent1")
# ent1_test <- ent1_melt[ent1_melt$OGID == "OG0002750", c("species", "ent1")]
fit_BMOU <- function(og, ent1_melt = ent2_melt, ent2_melt = ent2_melt, sp_ts_tre = sp_ts_tre, species_ll = species_ll) {
    print(og)
    ent1_test <- ent1_melt[ent1_melt$OGID == og, c("species", "ent1")]
    ent2_test <- ent2_melt[ent2_melt$OGID == og, c("species", "ent1")]
    ent1_test$regime <- 1
    ent1_test$regime[ent1_test$species %in% species_ll] <- 2
    ent1_test <- ent1_test[, c("species", "regime", "ent1")]
    ouwie_aicc <- c()
    model_set <- list()
    test_BM1 <- myOUwie(sp_ts_tre, ent1_test, model = "BM1", algorithm = "three.point") # clade = clades,
    if (class(test_BM1) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, BM1 = test_BM1$AICc)
        model_set <- c(model_set, list(BM1_fit = test_BM1))
    }
    test_BMS <- myOUwie(sp_ts_tre, ent1_test, model = "BMS", algorithm = "three.point") # , clade = clades
    if (class(test_BMS) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, BMS = test_BMS$AICc)
        model_set <- c(model_set, list(BMS_fit = test_BMS))
    }
    test_OU1 <- myOUwie(sp_ts_tre, ent1_test, model = "OU1", algorithm = "three.point") # , clade = clades
    if (class(test_OU1) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OU1 = test_OU1$AICc)
        model_set <- c(model_set, list(OU1_fit = test_OU1))
    }
    test_OUM <- myOUwie(sp_ts_tre, ent1_test, model = "OUM", algorithm = "three.point") # , clade = clades
    if (class(test_OUM) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OUM = test_OUM$AICc)
        model_set <- c(model_set, list(OUM_fit = test_OUM))
    }
    test_OUMA <- myOUwie(sp_ts_tre, ent1_test, model = "OUMA", algorithm = "three.point") # , clade = clades
    if (class(test_OUMA) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OUMA = test_OUMA$AICc)
        model_set <- c(model_set, list(OUMA_fit = test_OUMA))
    }
    test_OUMV <- myOUwie(sp_ts_tre, ent1_test, model = "OUMV", algorithm = "three.point") # , clade = clades
    if (class(test_OUMV) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OUMV = test_OUMV$AICc)
        model_set <- c(model_set, list(OUMV_fit = test_OUMV))
    }
    test_OUMVA <- myOUwie(sp_ts_tre, ent1_test, model = "OUMVA", algorithm = "three.point") # , clade = clades
    if (class(test_OUMVA) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OUMVA = test_OUMVA$AICc)
        model_set <- c(model_set, list(OUMVA_fit = test_OUMVA))
    }
    aicwt <- aic.w(ouwie_aicc)
    ## choose the biggest aicwt model
    # aicwt_df <- rbind(aicwt_df, data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)]))
    model_set <- list(og = model_set, aicwt = data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)]))
    names(model_set) <- c(og, "aicwt")
    # res <- c(res, model_set)
    return(model_set)
}

res <- alply(test_og$OGID, 1, fit_BMOU, ent2_melt = ent2_melt, sp_ts_tre = sp_ts_tre, species_ll = species_ll, .parallel = TRUE)
aicwt_df <- ldply(res, function(x) x$aicwt)
res <- lapply(res, function(x) x[1])
names(res) <- test_og$OGID
# aicwt_df <- data.frame()
# res <- list()
# for (og in test_og$OGID) {
#     # res <- alply(test_og$OGID, 1, function(og) {
#     print(og)
#     ent1_test <- ent1_melt[ent1_melt$OGID == og, c("species", "ent1")]
#     ent2_test <- ent2_melt[ent2_melt$OGID == og, c("species", "ent1")]
#     ent1_test$regime <- 1
#     ent1_test$regime[ent1_test$species %in% species_ll] <- 2
#     ent1_test <- ent1_test[, c("species", "regime", "ent1")]
#     # if (length(which(is.na(ent1_test$ent1[ent1_test$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
#     #     next
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Myotis_chinensis"])) {
#     #     clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Rhinolophus_pusillus"])) {
#     #     clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Hipposideros_larvatus"])) {
#     #     clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
#     # } else {
#     #     clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
#     # }
#     ouwie_aicc <- c()
#     model_set <- list()
#     test_BM1 <- myOUwie(sp_ts_tre, ent1_test, model = "BM1", algorithm = "three.point") # clade = clades,
#     if (class(test_BM1) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, BM1 = test_BM1$AICc)
#         model_set <- c(model_set, list(BM1_fit = test_BM1))
#     }
#     test_BMS <- myOUwie(sp_ts_tre, ent1_test, model = "BMS", algorithm = "three.point") # , clade = clades
#     if (class(test_BMS) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, BMS = test_BMS$AICc)
#         model_set <- c(model_set, list(BMS_fit = test_BMS))
#     }
#     test_OU1 <- myOUwie(sp_ts_tre, ent1_test, model = "OU1", algorithm = "three.point") # , clade = clades
#     if (class(test_OU1) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OU1 = test_OU1$AICc)
#         model_set <- c(model_set, list(OU1_fit = test_OU1))
#     }
#     test_OUM <- myOUwie(sp_ts_tre, ent1_test, model = "OUM", algorithm = "three.point") # , clade = clades
#     if (class(test_OUM) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUM = test_OUM$AICc)
#         model_set <- c(model_set, list(OUM_fit = test_OUM))
#     }
#     test_OUMA <- myOUwie(sp_ts_tre, ent1_test, model = "OUMA", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMA) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMA = test_OUMA$AICc)
#         model_set <- c(model_set, list(OUMA_fit = test_OUMA))
#     }
#     test_OUMV <- myOUwie(sp_ts_tre, ent1_test, model = "OUMV", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMV) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMV = test_OUMV$AICc)
#         model_set <- c(model_set, list(OUMV_fit = test_OUMV))
#     }
#     test_OUMVA <- myOUwie(sp_ts_tre, ent1_test, model = "OUMVA", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMVA) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMVA = test_OUMVA$AICc)
#         model_set <- c(model_set, list(OUMVA_fit = test_OUMVA))
#     }
#     # model_set <- list(BM1_fit = test_BM1, BMS_fit = test_BMS, OU1_fit = test_OU1, OUM_fit = test_OUM, OUMV_fit = test_OUMV, OUMVA_fit = test_OUMVA)
#     # print(getModelTable(model_set))
#     # aic.w(unlist(model_set))
#     # ouwie_aicc <- c(test_BM1$AICc, test_BMS$AICc, test_OU1$AICc, test_OUM$AICc, test_OUMV$AICc, test_OUMA$AICc, test_OUMVA$AICc)
#     # names(ouwie_aicc) <- c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")[1:length(ouwie_aicc)]
#     aicwt <- aic.w(ouwie_aicc)
#     ## choose the biggest aicwt model
#     aicwt_df <- rbind(aicwt_df, data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)]))
#     # aicwt_df <- data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)])
#     # names(aicwt)[which.max(aicwt)]
#     model_set <- list(og = model_set)
#     names(model_set) <- og
#     res <- c(res, model_set)
#     # return(c(model_set, aicwt_df))
# } # , .parallel = TRUE)
if (!dir.exists(paste0("/media/Data/zhangz/chip/analysis/summary2/ent_abc_ou/longevity/"))) {
    dir.create(paste0("/media/Data/zhangz/chip/analysis/summary2/ent_abc_ou/longevity/"))
}
write.csv(aicwt_df, paste0("aicwt_df_abc_ll.csv"), row.names = FALSE)
saveRDS(res, paste0("res_abc_ll.rds"))
# res <- readRDS("res_abc.rds")


# ortho_me_adj[ortho_me_adj$id == "OG0003643", "Mus_musculus"]
# ortho_me_adj[ortho_me_adj$id == "OG0002721", "Mus_musculus"]
lapply(res, function(x) {
    print(x$test_OUMVA$AICc)
})
print("result of ou fit for ent1 and species with longevith > 26")
print(tapply(aicwt_df$model, aicwt_df$model, length))

# aicwt_df_30 <- data.frame()
# res_30 <- list()
# for (og in test_og$OGID) {
#     # res <- alply(test_og$OGID, 1, function(og) {
#     print(og)
#     ent1_test <- ent1_melt[ent1_melt$OGID == og, c("species", "ent1")]
#     ent2_test <- ent2_melt[ent2_melt$OGID == og, c("species", "ent1")]
#     ent1_test$regime <- 1
#     ent1_test$regime[ent1_test$species %in% species_ll_30] <- 2
#     ent1_test <- ent1_test[, c("species", "regime", "ent1")]
#     # if (length(which(is.na(ent1_test$ent1[ent1_test$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
#     #     next
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Myotis_chinensis"])) {
#     #     clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Rhinolophus_pusillus"])) {
#     #     clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Hipposideros_larvatus"])) {
#     #     clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
#     # } else {
#     #     clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
#     # }
#     ouwie_aicc <- c()
#     model_set <- list()
#     test_BM1 <- myOUwie(sp_ts_tre, ent1_test, model = "BM1", algorithm = "three.point") # clade = clades,
#     if (class(test_BM1) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, BM1 = test_BM1$AICc)
#         model_set <- c(model_set, list(BM1_fit = test_BM1))
#     }
#     test_BMS <- myOUwie(sp_ts_tre, ent1_test, model = "BMS", algorithm = "three.point") # , clade = clades
#     if (class(test_BMS) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, BMS = test_BMS$AICc)
#         model_set <- c(model_set, list(BMS_fit = test_BMS))
#     }
#     test_OU1 <- myOUwie(sp_ts_tre, ent1_test, model = "OU1", algorithm = "three.point") # , clade = clades
#     if (class(test_OU1) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OU1 = test_OU1$AICc)
#         model_set <- c(model_set, list(OU1_fit = test_OU1))
#     }
#     test_OUM <- myOUwie(sp_ts_tre, ent1_test, model = "OUM", algorithm = "three.point") # , clade = clades
#     if (class(test_OUM) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUM = test_OUM$AICc)
#         model_set <- c(model_set, list(OUM_fit = test_OUM))
#     }
#     test_OUMA <- myOUwie(sp_ts_tre, ent1_test, model = "OUMA", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMA) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMA = test_OUMA$AICc)
#         model_set <- c(model_set, list(OUMA_fit = test_OUMA))
#     }
#     test_OUMV <- myOUwie(sp_ts_tre, ent1_test, model = "OUMV", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMV) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMV = test_OUMV$AICc)
#         model_set <- c(model_set, list(OUMV_fit = test_OUMV))
#     }
#     test_OUMVA <- myOUwie(sp_ts_tre, ent1_test, model = "OUMVA", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMVA) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMVA = test_OUMVA$AICc)
#         model_set <- c(model_set, list(OUMVA_fit = test_OUMVA))
#     }
#     # model_set <- list(BM1_fit = test_BM1, BMS_fit = test_BMS, OU1_fit = test_OU1, OUM_fit = test_OUM, OUMV_fit = test_OUMV, OUMVA_fit = test_OUMVA)
#     # print(getModelTable(model_set))
#     # aic.w(unlist(model_set))
#     # ouwie_aicc <- c(test_BM1$AICc, test_BMS$AICc, test_OU1$AICc, test_OUM$AICc, test_OUMV$AICc, test_OUMA$AICc, test_OUMVA$AICc)
#     # names(ouwie_aicc) <- c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")[1:length(ouwie_aicc)]
#     aicwt <- aic.w(ouwie_aicc)
#     ## choose the biggest aicwt model
#     aicwt_df_30 <- rbind(aicwt_df_30, data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)]))
#     # aicwt_df <- data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)])
#     # names(aicwt)[which.max(aicwt)]
#     model_set <- list(og = model_set)
#     names(model_set) <- og
#     res_30 <- c(res_30, model_set)
#     # return(c(model_set, aicwt_df))
# }
res_30 <- alply(test_og$OGID, 1, fit_BMOU, ent2_melt = ent2_melt, sp_ts_tre = sp_ts_tre, species_ll = species_ll_30, .parallel = TRUE)
aicwt_df_30 <- ldply(res_30, function(x) x$aicwt)
res_30 <- lapply(res_30, function(x) x[1])
names(res_30) <- test_og$OGID
write.csv(aicwt_df_30, paste0("aicwt_df_abc_ll_30.csv"), row.names = FALSE)
saveRDS(res_30, paste0("res_abc_ll_30.rds"))
print("result of ou fit for ent1 and species with longevith > 30")
print(tapply(aicwt_df_30$model, aicwt_df_30$model, length))

# ent1_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ent_ouwie/longevity/ent1_all.csv", header = TRUE)
# ent2_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ent_ouwie/longevity/ent2_all.csv", header = TRUE)

# library(preprocessCore)

# ent1_all_qn <- normalize.quantiles(as.matrix(ent1_all[, -1])) %>%
#     as.data.frame() %>%
#     cbind(ent1_all$OGID, .)
# colnames(ent1_all_qn) <- colnames(ent1_all)
# ent2_all_qn <- normalize.quantiles(as.matrix(ent2_all[, -1])) %>%
#     as.data.frame() %>%
#     cbind(ent2_all$OGID, .)
# colnames(ent2_all_qn) <- colnames(ent2_all)
# write.csv(ent1_all_qn, "ent1_all_qn.csv", row.names = FALSE)
# write.csv(ent2_all_qn, "ent2_all_qn.csv", row.names = FALSE)
ent1_all_qn <- read.csv("ent1_all_qn.csv", header = TRUE)
ent2_all_qn <- read.csv("ent2_all_qn.csv", header = TRUE)


# ent_phylosig_qn <- data.frame(OGID = ogid20$x)
ent_phylosig_qn <- adply(seq_len(nrow(ent1_all_qn)), 1, calc_physig, ent1_all = ent1_all_qn, ent2_all = ent2_all_qn, .parallel = TRUE)
# for (i in seq_len(nrow(ent1_all_qn))) {
#     ent1 <- as.numeric(ent1_all_qn[i, 2:22])
#     names(ent1) <- colnames(ent1_all_qn)[2:22]
#     ent2 <- as.numeric(ent2_all_qn[i, 2:22])
#     names(ent2) <- colnames(ent2_all_qn)[2:22]
#     if (length(which(ent1 != 0 & !is.na(ent1))) != 0) {
#         ent1_k <- phylosig(sp_tre, ent1, method = "K", test = TRUE)
#         ent1_lambda <- phylosig(sp_tre, ent1, method = "lambda", test = TRUE)
#         ent_phylosig_qn[i, "ent1_K"] <- ent1_k$K
#         ent_phylosig_qn[i, "ent1_K_P"] <- ent1_k$P
#         ent_phylosig_qn[i, "ent1_lambda"] <- ent1_lambda$lambda
#         ent_phylosig_qn[i, "ent1_lambda_P"] <- ent1_lambda$P
#     } else {
#         ent_phylosig_qn[i, "ent1_K"] <- NA
#         ent_phylosig_qn[i, "ent1_K_P"] <- NA
#         ent_phylosig_qn[i, "ent1_lambda"] <- NA
#         ent_phylosig_qn[i, "ent1_lambda_P"] <- NA
#     }
#     if (length(which(ent2 != 0 & !is.na(ent2))) != 0) {
#         ent2_k <- phylosig(sp_tre, ent2, method = "K", test = TRUE)
#         ent2_lambda <- phylosig(sp_tre, ent2, method = "lambda", test = TRUE)
#         ent_phylosig_qn[i, "ent2_K"] <- ent2_k$K
#         ent_phylosig_qn[i, "ent2_K_P"] <- ent2_k$P
#         ent_phylosig_qn[i, "ent2_lambda"] <- ent2_lambda$lambda
#         ent_phylosig_qn[i, "ent2_lambda_P"] <- ent2_lambda$P
#     } else {
#         ent_phylosig_qn[i, "ent2_K"] <- NA
#         ent_phylosig_qn[i, "ent2_K_P"] <- NA
#         ent_phylosig_qn[i, "ent2_lambda"] <- NA
#         ent_phylosig_qn[i, "ent2_lambda_P"] <- NA
#     }
# }
write.csv(ent_phylosig_qn, "ent_phylosig_abc_qn.csv", row.names = FALSE)

ent1_qn_melt <- reshape2::melt(ent1_all_qn, id.vars = "OGID", variable.name = "species", value.name = "ent1")
ent2_qn_melt <- reshape2::melt(ent2_all_qn, id.vars = "OGID", variable.name = "species", value.name = "ent1")
test_og_qn <- ent_phylosig_qn[ent_phylosig_qn$ent1_K_P < 0.05 & ent_phylosig_qn$ent1_K >= 1, ] %>% drop_na()
print(paste0("phylosig og num is : ", nrow(test_og_qn)))
# aicwt_df_30_qn <- data.frame()
# res_30_qn <- list()
# for (og in test_og_qn$OGID) {
#     # res <- alply(test_og$OGID, 1, function(og) {
#     print(og)
#     ent1_test <- ent1_qn_melt[ent1_qn_melt$OGID == og, c("species", "ent1")]
#     ent2_test <- ent2_qn_melt[ent2_qn_melt$OGID == og, c("species", "ent1")]
#     ent1_test$regime <- 1
#     ent1_test$regime[ent1_test$species %in% species_ll_30] <- 2
#     ent1_test <- ent1_test[, c("species", "regime", "ent1")]
#     # if (length(which(is.na(ent1_test$ent1[ent1_test$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
#     #     next
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Myotis_chinensis"])) {
#     #     clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Rhinolophus_pusillus"])) {
#     #     clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
#     # } else if (is.na(ent1_test$ent1[ent1_test$species == "Hipposideros_larvatus"])) {
#     #     clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
#     # } else {
#     #     clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
#     # }
#     ouwie_aicc <- c()
#     model_set <- list()
#     test_BM1 <- myOUwie(sp_ts_tre, ent1_test, model = "BM1", algorithm = "three.point") # clade = clades,
#     if (class(test_BM1) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, BM1 = test_BM1$AICc)
#         model_set <- c(model_set, list(BM1_fit = test_BM1))
#     }
#     test_BMS <- myOUwie(sp_ts_tre, ent1_test, model = "BMS", algorithm = "three.point") # , clade = clades
#     if (class(test_BMS) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, BMS = test_BMS$AICc)
#         model_set <- c(model_set, list(BMS_fit = test_BMS))
#     }
#     test_OU1 <- myOUwie(sp_ts_tre, ent1_test, model = "OU1", algorithm = "three.point") # , clade = clades
#     if (class(test_OU1) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OU1 = test_OU1$AICc)
#         model_set <- c(model_set, list(OU1_fit = test_OU1))
#     }
#     test_OUM <- myOUwie(sp_ts_tre, ent1_test, model = "OUM", algorithm = "three.point") # , clade = clades
#     if (class(test_OUM) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUM = test_OUM$AICc)
#         model_set <- c(model_set, list(OUM_fit = test_OUM))
#     }
#     test_OUMA <- myOUwie(sp_ts_tre, ent1_test, model = "OUMA", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMA) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMA = test_OUMA$AICc)
#         model_set <- c(model_set, list(OUMA_fit = test_OUMA))
#     }
#     test_OUMV <- myOUwie(sp_ts_tre, ent1_test, model = "OUMV", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMV) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMV = test_OUMV$AICc)
#         model_set <- c(model_set, list(OUMV_fit = test_OUMV))
#     }
#     test_OUMVA <- myOUwie(sp_ts_tre, ent1_test, model = "OUMVA", algorithm = "three.point") # , clade = clades
#     if (class(test_OUMVA) != "try-error") {
#         ouwie_aicc <- c(ouwie_aicc, OUMVA = test_OUMVA$AICc)
#         model_set <- c(model_set, list(OUMVA_fit = test_OUMVA))
#     }
#     # model_set <- list(BM1_fit = test_BM1, BMS_fit = test_BMS, OU1_fit = test_OU1, OUM_fit = test_OUM, OUMV_fit = test_OUMV, OUMVA_fit = test_OUMVA)
#     # print(getModelTable(model_set))
#     # aic.w(unlist(model_set))
#     # ouwie_aicc <- c(test_BM1$AICc, test_BMS$AICc, test_OU1$AICc, test_OUM$AICc, test_OUMV$AICc, test_OUMA$AICc, test_OUMVA$AICc)
#     # names(ouwie_aicc) <- c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")[1:length(ouwie_aicc)]
#     aicwt <- aic.w(ouwie_aicc)
#     ## choose the biggest aicwt model
#     aicwt_df_30_qn <- rbind(aicwt_df_30_qn, data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)]))
#     # aicwt_df <- data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)])
#     # names(aicwt)[which.max(aicwt)]
#     model_set <- list(og = model_set)
#     names(model_set) <- og
#     res_30_qn <- c(res_30_qn, model_set)
#     # return(c(model_set, aicwt_df))
# }
res_30_qn <- alply(test_og_qn$OGID, 1, fit_BMOU, ent1_melt = ent1_qn_melt, ent2_melt = ent2_qn_melt, sp_ts_tre = sp_ts_tre, species_ll = species_ll_30, .parallel = TRUE)
aicwt_df_30_qn <- ldply(res_30_qn, function(x) x$aicwt)
res_30_qn <- lapply(res_30_qn, function(x) x[1])
names(res_30_qn) <- test_og_qn$OGID
write.csv(aicwt_df_30_qn, paste0("aicwt_df_abc_ll_30.csv"), row.names = FALSE)
saveRDS(res_30_qn, paste0("res_abc_ll_30.rds"))
print("result of ou fit for qn ent1 and species with longevith > 30")
print(tapply(aicwt_df_30_qn$model, aicwt_df_30_qn$model, length))

res_26_qn <- alply(test_og_qn$OGID, 1, fit_BMOU, ent1_melt = ent1_qn_melt, ent2_melt = ent2_qn_melt, sp_ts_tre = sp_ts_tre, species_ll = species_ll, .parallel = TRUE)
aicwt_df_26_qn <- ldply(res_26_qn, function(x) x$aicwt)
res_26_qn <- lapply(res_26_qn, function(x) x[1])
names(res_26_qn) <- test_og_qn$OGID
write.csv(aicwt_df_26_qn, paste0("aicwt_df_abc_ll_26.csv"), row.names = FALSE)
saveRDS(res_26_qn, paste0("res_abc_ll_26.rds"))
print("result of ou fit for qn ent1 and species with longevith > 26")
print(tapply(aicwt_df_26_qn$model, aicwt_df_26_qn$model, length))
