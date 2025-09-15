## fit OU, BM and etc. models to expression data, to find altered expressed genes in bats

## load necessary libraries
library(OUwie)
library(plyr)
library(ape)
library(doMC)
doMC::registerDoMC(cores = 40)

## load expression data
setwd("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/OUwie")
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/fpkm_ortho.RData")
adj_FPKM_me$species <- gsub("Macaca_fascicularis", "Macaca_mulatta", adj_FPKM_me$species)
adj_brain_me <- adj_FPKM_me[adj_FPKM_me$tissue == "Brain", ]
adj_kidney_me <- adj_FPKM_me[adj_FPKM_me$tissue == "Kidney", ]
adj_liver_me <- adj_FPKM_me[adj_FPKM_me$tissue == "Liver", ]
# calc mean fpkm value of different rep in adj_brain_me
brain_mean_fpkm <- aggregate(fpkm ~ OGID + species, adj_brain_me, mean)
kidney_mean_fpkm <- aggregate(fpkm ~ OGID + species, adj_kidney_me, mean)
liver_mean_fpkm <- aggregate(fpkm ~ OGID + species, adj_liver_me, mean)

adj_b_kl_me <- merge(brain_mean_fpkm, kidney_mean_fpkm, by = c("OGID", "species")) %>%
    merge(., liver_mean_fpkm, by = c("OGID", "species"))
adj_b_kl_me$b_kl <- adj_b_kl_me$fpkm.x - ((adj_b_kl_me$fpkm.y + adj_b_kl_me$fpkm) / 2)
adj_b_kl_me <- adj_b_kl_me[, c("OGID", "species", "b_kl")]
# colnames(adj_b_kl_me)[3] <- "b_kl"
# searching for genes which brain is lower than the mean of kidney and liver,
# and the difference is bigger in bats than other species
library(OUwie)
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
sp_tre <- ape::read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
sp_ts_tre <- keep.tip(sp_tre, species_ts)
myOUwie <- function(tree, trait, model, ...) {
    not_na <- !is.na(trait$fpkm)
    tmp_tree <- keep.tip(tree, as.vector(trait$species[not_na]))
    tmp_trait <- trait[not_na, ]
    tmp_res <- try(OUwie(tmp_tree, tmp_trait, model = model, ...))
    return(tmp_res)
}
# aicwt_df <- data.frame()
res <- list()
# for (og in ortho_me_20_adj$id) {
b_kl_ou_res <- adply(ortho_me_20_adj$id, 1, function(og) {
    diff_tmp <- adj_b_kl_me[adj_b_kl_me$OGID == og, c("species", "b_kl")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "b_kl")]
    if (length(which(diff_tmp$b_kl > 0)) < 10) {
        return(data.frame(og = og, model = NA, aicwt = NA))
    }
    if (length(which(is.na(diff_tmp$b_kl[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    ouwie_aicc <- c()
    model_set <- list()
    test_BM1 <- myOUwie(sp_tre, diff_tmp, model = "BM1", clade = clades, algorithm = "three.point")
    if (class(test_BM1) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, BM1 = test_BM1$AICc)
        model_set <- c(model_set, list(BM1 = test_BM1))
    }
    # test_BMS <- myOUwie(sp_tre, diff_tmp, model = "BMS", clade = clades, algorithm = "three.point")
    # if (class(test_BMS) != "try-error") {
    #     ouwie_aicc <- c(ouwie_aicc, BMS = test_BMS$AICc)
    #     model_set <- c(model_set, list(BMS = test_BMS))
    # }
    test_OU1 <- myOUwie(sp_tre, diff_tmp, model = "OU1", clade = clades, algorithm = "three.point")
    if (class(test_OU1) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OU1 = test_OU1$AICc)
        model_set <- c(model_set, list(OU1 = test_OU1))
    }
    test_OUM <- myOUwie(sp_tre, diff_tmp, model = "OUM", clade = clades, algorithm = "three.point")
    if (class(test_OUM) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OUM = test_OUM$AICc)
        model_set <- c(model_set, list(OUM = test_OUM))
    }
    # test_OUMA <- myOUwie(sp_tre, diff_tmp, model = "OUMA", clade = clades, algorithm = "three.point")
    # if (class(test_OUMA) != "try-error") {
    #     ouwie_aicc <- c(ouwie_aicc, OUMA = test_OUMA$AICc)
    #     model_set <- c(model_set, list(OUMA = test_OUMA))
    # }
    # test_OUMV <- myOUwie(sp_tre, diff_tmp, model = "OUMV", clade = clades, algorithm = "three.point")
    # if (class(test_OUMV) != "try-error") {
    #     ouwie_aicc <- c(ouwie_aicc, OUMV = test_OUMV$AICc)
    #     model_set <- c(model_set, list(OUMV = test_OUMV))
    # }
    # test_OUMVA <- myOUwie(sp_tre, diff_tmp, model = "OUMVA", clade = clades, algorithm = "three.point")
    # if (class(test_OUMVA) != "try-error") {
    #     ouwie_aicc <- c(ouwie_aicc, OUMVA = test_OUMVA$AICc)
    #     model_set <- c(model_set, list(OUMVA = test_OUMVA))
    # }
    if (length(ouwie_aicc) == 0) {
        return(data.frame(og = og, model = NA, aicwt = NA))
    }
    aicwt <- aic.w(ouwie_aicc)
    ## choose the biggest aicwt model
    # aicwt_df <- rbind(aicwt_df, data.frame(OGID = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)]))
    # names(aicwt)[which.max(aicwt)]
    # model_set <- list(og = model_set)
    # names(model_set) <- og
    # res <- c(res, model_set)
    # res <- model_set[[which.max(aicwt)]]
    new_res <- model_set[[which.max(aicwt)]]
    res <- c(res, new_res)
    names(res)[length(res)] <- og
    res_df <- data.frame(og = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)])
    return(res_df)
}, .parallel = TRUE, .id = NULL)
saveRDS(res, "ouwie_res_model.RDS")
# write.csv(ou_res, "ou_res.csv", row.names = FALSE)
write.csv(b_kl_ou_res, "b_kl_ou_res.csv", row.names = FALSE)

## check results
res <- readRDS("ouwie_res_model.RDS")
# ou_res <- read.csv("ou_res.csv", header = TRUE)

b_kl_ou_res <- read.csv("b_kl_ou_res.csv", header = TRUE)
tapply(ou_res$model, ou_res$model, length)
#  BM1  BMS  OU1  OUM OUMA OUMV
# 1076  380 3335  646   73  517

# r$> tapply(ou_res$model, ou_res$model, length)
#  BM1  BMS  OU1  OUM OUMA OUMV
#  670  252 1808  217   44  201
# brain - (kidney + liver)/2
tapply(b_kl_ou_res$model, b_kl_ou_res$model, length)
#  BM1  OU1  OUM
#  824 2068  300
head(ou_res[ou_res$model == "OUM", ])
# oumx <- ou_res[ou_res$model %in% c("OUM", "OUMA", "OUMV"), ]
b_kl_oumx <- b_kl_ou_res[which(b_kl_ou_res$model == "OUM"), ]
res <- list()
# for (i in 1:nrow(oumx)) {
res <- llply(1:nrow(oumx), function(i) {
    og <- oumx$og[i]
    model <- oumx$model[i]
    diff_tmp <- adj_b_kl_me[adj_b_kl_me$OGID == og, c("species", "b_kl")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "b_kl")]
    if (length(which(diff_tmp$b_kl > 0)) < 10) {
        next
    }
    if (length(which(is.na(diff_tmp$b_kl[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    test_model <- myOUwie(sp_tre, diff_tmp, model = model, clade = clades, algorithm = "three.point")
    # res <- c(res, test_model)
    return(test_model)
}, .parallel = TRUE)
names(res) <- oumx$og

b_kl_res <- list()
b_kl_res <- llply(1:nrow(b_kl_oumx), function(i) {
    og <- b_kl_oumx$og[i]
    model <- b_kl_oumx$model[i]
    diff_tmp <- adj_b_kl_me[adj_b_kl_me$OGID == og, c("species", "b_kl")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "b_kl")]
    if (length(which(diff_tmp$b_kl > 0)) < 10) {
        next
    }
    if (length(which(is.na(diff_tmp$b_kl[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$b_kl[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    test_model <- myOUwie(sp_tre, diff_tmp, model = model, clade = clades, algorithm = "three.point")
    # res <- c(res, test_model)
    return(test_model)
}, .parallel = TRUE)
names(b_kl_res) <- b_kl_oumx$og
saveRDS(b_kl_res, "b_kl_res.RDS")
b_kl_res <- readRDS("b_kl_res.RDS")

# fit ou and oum model to three tissues respectively
# brain_mean_fpkm, kidney_mean_fpkm, liver_mean_fpkm
b_res <- list()
b_ou_res <- adply(ortho_me_20_adj$id, 1, function(og) {
    diff_tmp <- brain_mean_fpkm[brain_mean_fpkm$OGID == og, c("species", "fpkm")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "fpkm")]
    if (length(which(diff_tmp$fpkm > 0)) < 10) {
        return(data.frame(og = og, model = NA, aicwt = NA))
    }
    if (length(which(is.na(diff_tmp$fpkm[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    ouwie_aicc <- c()
    model_set <- list()
    # test_BM1 <- myOUwie(sp_tre, diff_tmp, model = "BM1", clade = clades, algorithm = "three.point")
    # if (class(test_BM1) != "try-error") {
    #     ouwie_aicc <- c(ouwie_aicc, BM1 = test_BM1$AICc)
    #     model_set <- c(model_set, list(BM1 = test_BM1))
    # }
    test_OU1 <- myOUwie(sp_tre, diff_tmp, model = "OU1", clade = clades, algorithm = "three.point")
    if (class(test_OU1) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OU1 = test_OU1$AICc)
        model_set <- c(model_set, list(OU1 = test_OU1))
    }
    test_OUM <- myOUwie(sp_tre, diff_tmp, model = "OUM", clade = clades, algorithm = "three.point")
    if (class(test_OUM) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OUM = test_OUM$AICc)
        model_set <- c(model_set, list(OUM = test_OUM))
    }
    if (length(ouwie_aicc) == 0) {
        return(data.frame(og = og, model = NA, aicwt = NA))
    }
    aicwt <- aic.w(ouwie_aicc)
    ## choose the biggest aicwt model
    new_res <- model_set[[which.max(aicwt)]]
    b_res <- c(b_res, new_res)
    names(b_res)[length(b_res)] <- og
    res_df <- data.frame(og = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)])
    return(res_df)
}, .parallel = TRUE, .id = NULL)
# saveRDS(res, "ouwie_res_model.RDS")
write.csv(b_ou_res, "b_ou_res.csv", row.names = FALSE)
b_oumx <- b_ou_res[which(b_ou_res$model == "OUM"), ]
tapply(b_ou_res$model, b_ou_res$model, length)

k_res <- list()
k_ou_res <- adply(ortho_me_20_adj$id, 1, function(og) {
    diff_tmp <- kidney_mean_fpkm[kidney_mean_fpkm$OGID == og, c("species", "fpkm")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "fpkm")]
    if (length(which(diff_tmp$fpkm > 0)) < 10) {
        return(data.frame(og = og, model = NA, aicwt = NA))
    }
    if (length(which(is.na(diff_tmp$fpkm[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    ouwie_aicc <- c()
    model_set <- list()
    # test_BM1 <- myOUwie(sp_tre, diff_tmp, model = "BM1", clade = clades, algorithm = "three.point")
    # if (class(test_BM1) != "try-error") {
    #     ouwie_aicc <- c(ouwie_aicc, BM1 = test_BM1$AICc)
    #     model_set <- c(model_set, list(BM1 = test_BM1))
    # }
    test_OU1 <- myOUwie(sp_tre, diff_tmp, model = "OU1", clade = clades, algorithm = "three.point")
    if (class(test_OU1) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OU1 = test_OU1$AICc)
        model_set <- c(model_set, list(OU1 = test_OU1))
    }
    test_OUM <- myOUwie(sp_tre, diff_tmp, model = "OUM", clade = clades, algorithm = "three.point")
    if (class(test_OUM) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OUM = test_OUM$AICc)
        model_set <- c(model_set, list(OUM = test_OUM))
    }
    if (length(ouwie_aicc) == 0) {
        return(data.frame(og = og, model = NA, aicwt = NA))
    }
    aicwt <- aic.w(ouwie_aicc)
    ## choose the biggest aicwt model
    new_res <- model_set[[which.max(aicwt)]]
    K_res <- c(k_res, new_res)
    names(k_res)[length(k_res)] <- og
    res_df <- data.frame(og = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)])
    return(res_df)
}, .parallel = TRUE, .id = NULL)
write.csv(k_ou_res, "k_ou_res.csv", row.names = FALSE)
k_oumx <- k_ou_res[which(k_ou_res$model == "OUM"), ]

l_res <- list()
l_ou_res <- adply(ortho_me_20_adj$id, 1, function(og) {
    diff_tmp <- liver_mean_fpkm[liver_mean_fpkm$OGID == og, c("species", "fpkm")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "fpkm")]
    if (length(which(diff_tmp$fpkm > 0)) < 10) {
        return(data.frame(og = og, model = NA, aicwt = NA))
    }
    if (length(which(is.na(diff_tmp$fpkm[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    ouwie_aicc <- c()
    model_set <- list()
    # test_BM1 <- myOUwie(sp_tre, diff_tmp, model = "BM1", clade = clades, algorithm = "three.point")
    # if (class(test_BM1) != "try-error") {
    #     ouwie_aicc <- c(ouwie_aicc, BM1 = test_BM1$AICc)
    #     model_set <- c(model_set, list(BM1 = test_BM1))
    # }
    test_OU1 <- myOUwie(sp_tre, diff_tmp, model = "OU1", clade = clades, algorithm = "three.point")
    if (class(test_OU1) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OU1 = test_OU1$AICc)
        model_set <- c(model_set, list(OU1 = test_OU1))
    }
    test_OUM <- myOUwie(sp_tre, diff_tmp, model = "OUM", clade = clades, algorithm = "three.point")
    if (class(test_OUM) != "try-error") {
        ouwie_aicc <- c(ouwie_aicc, OUM = test_OUM$AICc)
        model_set <- c(model_set, list(OUM = test_OUM))
    }
    if (length(ouwie_aicc) == 0) {
        return(data.frame(og = og, model = NA, aicwt = NA))
    }
    aicwt <- aic.w(ouwie_aicc)
    ## choose the biggest aicwt model
    new_res <- model_set[[which.max(aicwt)]]
    l_res <- c(l_res, new_res)
    names(l_res)[length(l_res)] <- og
    res_df <- data.frame(og = og, model = names(aicwt)[which.max(aicwt)], aicwt = aicwt[which.max(aicwt)])
    return(res_df)
}, .parallel = TRUE, .id = NULL)
write.csv(l_ou_res, "l_ou_res.csv", row.names = FALSE)
l_oumx <- l_ou_res[which(l_ou_res$model == "OUM"), ]

b_res <- list()
b_res <- llply(1:nrow(b_oumx), function(i) {
    og <- b_oumx$og[i]
    model <- b_oumx$model[i]
    diff_tmp <- brain_mean_fpkm[brain_mean_fpkm$OGID == og, c("species", "fpkm")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "fpkm")]
    if (length(which(diff_tmp$fpkm > 0)) < 10) {
        next
    }
    if (length(which(is.na(diff_tmp$fpkm[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    test_model <- myOUwie(sp_tre, diff_tmp, model = model, clade = clades, algorithm = "three.point")
    # res <- c(res, test_model)
    return(test_model)
}, .parallel = TRUE)
names(b_res) <- b_oumx$og
saveRDS(b_res, "b_res.RDS")
b_res <- readRDS("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/OUwie/b_res.RDS")

k_res <- list()
k_res <- llply(1:nrow(k_oumx), function(i) {
    og <- k_oumx$og[i]
    model <- k_oumx$model[i]
    diff_tmp <- kidney_mean_fpkm[kidney_mean_fpkm$OGID == og, c("species", "fpkm")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "fpkm")]
    if (length(which(diff_tmp$fpkm > 0)) < 10) {
        next
    }
    if (length(which(is.na(diff_tmp$fpkm[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    test_model <- myOUwie(sp_tre, diff_tmp, model = model, clade = clades, algorithm = "three.point")
    # res <- c(res, test_model)
    return(test_model)
}, .parallel = TRUE)
names(k_res) <- k_oumx$og
saveRDS(k_res, "k_res.RDS")
k_res <- readRDS("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/OUwie/k_res.RDS")

l_res <- list()
l_res <- llply(1:nrow(l_oumx), function(i) {
    og <- l_oumx$og[i]
    model <- l_oumx$model[i]
    diff_tmp <- liver_mean_fpkm[liver_mean_fpkm$OGID == og, c("species", "fpkm")]
    diff_tmp$regime <- 1
    diff_tmp$regime[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")] <- 2
    diff_tmp <- diff_tmp[, c("species", "regime", "fpkm")]
    if (length(which(diff_tmp$fpkm > 0)) < 10) {
        next
    }
    if (length(which(is.na(diff_tmp$fpkm[diff_tmp$species %in% c("Rhinolophus_pusillus", "Hipposideros_larvatus", "Myotis_chinensis")]))) > 1) {
        next
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Myotis_chinensis"])) {
        clades <- c("Rhinolophus_pusillus", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Rhinolophus_pusillus"])) {
        clades <- c("Myotis_chinensis", "Hipposideros_larvatus")
    } else if (is.na(diff_tmp$fpkm[diff_tmp$species == "Hipposideros_larvatus"])) {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    } else {
        clades <- c("Myotis_chinensis", "Rhinolophus_pusillus")
    }
    test_model <- myOUwie(sp_tre, diff_tmp, model = model, clade = clades, algorithm = "three.point")
    # res <- c(res, test_model)
    return(test_model)
}, .parallel = TRUE)
names(l_res) <- l_oumx$og
saveRDS(l_res, "l_res.RDS")
l_res <- readRDS("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/OUwie/l_res.RDS")

## some gene are associated with ERV silencing
## H3.3, H3-3A, H3-3B: OG0000069, CENPA: OG0003786
# ASF1 : OG0001915
# CHAF1A: OG0014158
# KDM1A: OG0009518
silen_gene <- c(
    "TRIM28", "ATRX", "DAXX", "SETDB1", "ZFP91", "YY1", "CBX5", "SUV39H1",
    "EHMT2", "NSD2", "TASOR", "MPHOSPH8", "PPHLN1", "MORC2", "SUMO2", "KDM1A",
    "CHAF1A", "HDAC1", "ASF1A", "SMARCAD1", "H3-3A", "DNMT3L", "SSRP1", "USP7",
    "ZSCAN4", "DUX4", "ZFP809", "ZFP932"
)
silen_gene <- c(
    "TRIM28", "ATRX", "DAXX", "SETDB1", "ZFP91", "YY1", "CBX5", "SUV39H1",
    "EHMT2", "NSD2", "TASOR", "MPHOSPH8", "PPHLN1", "MORC2", "SUMO2", "KDM1A",
    "CHAF1A", "HDAC1", "ASF1A", "SMARCAD1", "H3-3A", "DNMT3L", "SSRP1", "USP7",
    "ZSCAN4", "DUX4", "CHD3", "Zfp708", "Zfp91", "Zfp932", "Gm14391", "ZFP809"
)
# ZFP91 = ZNF91, CBX5 = HP1, SUV39H1 = SUV39H2, EHMT2 = EHMT1
# histone methyltransferases:
# SUV39H1/SUV39H2/SETDB1/EHMT1/EHMT2/NSD2
# HUSH: TASOR/MPHOSPH8/PPHLN1
# HUSH bind: MORC2/SETDB1
# SUMO: SUMO2
# histone chaperone: CHAF1A/ATRX/DAXX
# histone demethylase: KDM1A
# histone deacetylase: HDAC1/2
# histone chaperone isoform: ASF1A/B
# SWI/SNF like remodeler: SMARCAD1
# H3.3: H3-3A/H3-3B
# DNA methylation: DNMT3L
# Histone chaperone: SSRP1 = FACT
# recruited by FACT: USP7, repress MuERVL
# ZSCAN4 = ZSCAN4C, H3K27ac,H3K4me3, enhancer activity of MuERVL LTR
# OG0000001: ZFP809, ZFP819
# OG0000015: ZFP932

silen_og <- c(
    "OG0003719", "OG0003481", "OG0014170", "OG0006238",
    "OG0005559", "OG0000431", "OG0006241", "OG0001353", "OG0001239",
    "OG0006231", "OG0006829", "OG0009696", "OG0002089", "OG0003060",
    "OG0000230", "OG0009518", "OG0014158", "OG0001341", "OG0001915",
    "OG0011616", "OG0000069", "OG0013589", "OG0005211", "OG0007937",
    "OG0002698", "OG0016569", "OG0000001", "OG0000015"
)
silen_og <- c(
    "OG0003719", "OG0003481", "OG0014170", "OG0006238",
    "OG0005559", "OG0000431", "OG0006241", "OG0001353", "OG0001239",
    "OG0006231", "OG0006829", "OG0009696", "OG0002089", "OG0003060",
    "OG0000230", "OG0009518", "OG0014158", "OG0001341", "OG0001915",
    "OG0011616", "OG0000069", "OG0013589", "OG0005211", "OG0007937",
    "OG0002698", "OG0016569", "OG0000832", "OG0000001", "OG0005559",
    "OG0000015", "OG0002937", "OG0000001"
)
names(silen_og) <- silen_gene
silen_og[silen_og %in% b_kl_oumx$og]
silen_og[silen_og %in% l_oumx$og]

# r$> silen_og[silen_og %in% b_oumx$og]
#       KDM1A    SMARCAD1
# "OG0009518" "OG0011616"
# up
# r$> silen_og[silen_og %in% k_oumx$og]
#    SMARCAD1
# "OG0011616"
# down
# r$> silen_og[silen_og %in% l_oumx$og]
#       ZFP91        NSD2       Zfp91
# "OG0005559" "OG0006231" "OG0005559"
# down

b_kl_res[which(names(b_kl_res) %in% silen_og)]

ortho_all <- read.csv("/media/Data/zhangz/chip/analysis/summary2/ortho_all.csv", header = TRUE)
grep("Zfp809", ortho_all)
ortho_all[which(grepl("Zfp809", ortho_all$Mus_musculus)), ]
ortho_all[which(grepl("Zfp932", ortho_all$Mus_musculus)), "Orthogroup"]
adj_FPKM_me[adj_FPKM_me$OGID == "OG0003719", ]
ortho_all[which(grepl("Gm15446", ortho_all)), ]
ortho_all[which(grepl("Yy1", ortho_all$Mus_musculus)), "Orthogroup"]
ortho_all[which(grepl("Lsd1", ortho_all$Mus_musculus)), "Orthogroup"]

save.image("ouwie.RData")
load("/media/Data/zhangz/chip/analysis/expression_mean/zhangz/OUwie/ouwie.RData")
