## lets count the constituion of target gene annotation
## how many cres can be annotated by abc, how many can be annotated by
## ts-1m, how many genes can be annotated to a cre,
## and how many cre assigned to the same gene

## load libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(doMC)
library(stringr)
library(magrittr)
doMC::registerDoMC(cores = 4)

cat_list <- c("1M_w_ne", "1M_wo_ne", "abc_1M_w_ne", "abc_1M_wo_ne", "abc_w_ne", "abc_wo_ne", "nearest")
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

## read in data
# for example sp <- "Felis_catus"
anno_sum <- adply(species, 1, function(sp) {
    info <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno_new.csv"), header = TRUE, sep = ",")
    res <- mdply(expand.grid(tis = unique(info$tissue), ele = unique(info$element)), function(tis, ele) {
        info_tmp <- info[info$tissue == tis & info$element == ele, ]
        gene_num <- summary(info_tmp$overlap_num) %>%
            c() %>%
            as.data.frame() %>%
            t() %>%
            as.data.frame()
        colnames(gene_num) <- paste0("geneN_", c("min", "f_th", "median", "mean", "t_th", "max"))
        all_gene <- paste(info_tmp$overlap_gene, collapse = ";") %>%
            strsplit(";") %>%
            unlist() %>%
            unique()
        cre_num_list <- adply(all_gene, 1, function(x) {
            length(grep(paste0("(^|;)", x, "(;|$)"), info_tmp$overlap_gene))
        })
        write.csv(cre_num_list, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/anno_summary/", sp, "_", tis, "_", ele, "_cre_num_list.csv"), row.names = FALSE)
        cre_num <- summary(cre_num_list$V1) %>%
            c() %>%
            as.data.frame() %>%
            t() %>%
            as.data.frame()
        colnames(cre_num) <- paste0("creN_", c("min", "f_th", "median", "mean", "t_th", "max"))
        note_count <- count(info_tmp, "note_abc") %>%
            t() %>%
            as.data.frame()
        colnames(note_count) <- note_count[1, ]
        note_count <- note_count[-1, ]
        rdf <- data.frame(
            tis = tis, ele = ele
        )
        rdf <- cbind(rdf, gene_num)
        rdf <- cbind(rdf, cre_num)
        rdf <- cbind(rdf, note_count)
        rdf$cre_sum <- sum(as.numeric(note_count[1, ]))
        return(rdf)
    }, .parallel = TRUE)
    res$sp <- sp
    return(res)
}, .parallel = TRUE)

write.csv(anno_sum, "/media/Data/zhangz/chip/analysis/summary2/abc_1M/anno_summary.csv", row.names = FALSE)

anno_sum <- anno_sum[, -1]
anno_sum <- anno_sum[, c(23, 1:22)]
anno_sum[, 4:23] <- lapply(anno_sum[, 4:23], as.numeric)
anno_sum$abc_sum <- rowSums(anno_sum[, c("abc_1M_w_ne", "abc_1M_wo_ne", "abc_w_ne", "abc_wo_ne")])
anno_sum$abc_ratio <- anno_sum$abc_sum / anno_sum$cre_sum
