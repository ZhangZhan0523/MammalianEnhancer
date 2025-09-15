library(stringr)
library(plyr)
library(doMC)
library(magrittr)
doMC::registerDoMC(cores = 10)

test <- read.csv("/media/Data/zhangz/chip/analysis/Felis_catus/anno/peaks2/anno2/Felis_catus_Brain_enhancer_flank_anno.csv", header = T)
anno <- sapply(test$annotation, function(x) str_split(x, " ")[[1]][1])
unique(anno)
test2 <- read.csv("/media/Data/zhangz/chip/analysis/Myotis_chinensis/anno/peaks2/anno2/Myotis_chinensis_Brain_enhancer_flank_anno.csv", header = T)
anno2 <- sapply(test2$annotation, function(x) str_split(x, " ")[[1]][1])
unique(anno2)

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
classes <- c("LINE", "SINE", "DNA", "LTR")
composition <- adply(species, 1, function(sp) {
    if (sp == "Neophocaena_asiaeorientalis") {
        tissues <- c("Brain")
    } else if (sp == "Rhinopithecus_roxellana") {
        tissues <- c("Liver", "Kidney")
    } else {
        tissues <- c("Liver", "Kidney", "Brain")
    }
    rep_info <- mdply(expand.grid(tis = tissues, ele = c("enhancer", "promoter")), function(tis, ele) {
        res <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/rep_ele/", sp, "_", tis, "_", ele, "_rep_ele_class_count.csv"), header = TRUE)
        res <- res[, c("repeat_class", "ele_count", "ele_ratio", "gen_bg_count", "gen_bg_ratio")]
        other_sum <- colSums(res[!res$repeat_class %in% classes, c("ele_count", "ele_ratio", "gen_bg_count", "gen_bg_ratio")])
        res <- rbind(res[res$repeat_class %in% classes, ], c("others", other_sum))
        return(res)
    })
    rep_wide <- rep_info[, 1:5] %>%
        pivot_wider(names_from = c(tis, ele, repeat_class), values_from = c(ele_count, ele_ratio)) %>%
        as.data.frame()
    return(cbind(species = sp, rep_wide))
})
write.csv(composition, "/media/Data/zhangz/chip/analysis/summary2/rep_ele/rep_ele_composition.csv", row.names = F)
