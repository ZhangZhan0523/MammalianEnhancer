## try to capture the duality of enhancers and promoters
## find the CREs that are both enhancers and promoters in different
## tissues in the same species

library(bedtoolsr)
library(plyr)
library(doMC)
library(magrittr)
library(ggplot2)

registerDoMC(cores = 12)

## basic information
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

## read cre records, for example, mouse
sp <- "Mus_musculus"

## a function judge a overlap cre in different tissue is a enhancer/ promoter or both kind
judge_ep <- function(row) {
    if (row$V5 == row$V10) {
        res <- row$V5
    } else {
        res <- "ep"
    }
    return(res)
}
duality_df <- adply(species[which(!species %in% c("Neophocaena_asiaeorientalis", "Rhinopithecus_roxellana"))], 1, function(sp) {
    cre_info <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_anno.csv"), header = TRUE, sep = ",")
    brain_cre <- cre_info[cre_info$tissue == "Brain", c("chr", "start", "end", "unique_id", "element")]
    kidney_cre <- cre_info[cre_info$tissue == "Kidney", c("chr", "start", "end", "unique_id", "element")]
    liver_cre <- cre_info[cre_info$tissue == "Liver", c("chr", "start", "end", "unique_id", "element")]
    bk_b_cre <- bt.intersect(a = brain_cre, b = kidney_cre, f = 0.5, e = TRUE, wa = TRUE, wb = TRUE)
    bk_k_cre <- bt.intersect(a = kidney_cre, b = brain_cre, f = 0.5, e = TRUE, wa = TRUE, wb = TRUE)
    bl_b_cre <- bt.intersect(a = brain_cre, b = liver_cre, f = 0.5, e = TRUE, wa = TRUE, wb = TRUE)
    bl_l_cre <- bt.intersect(a = liver_cre, b = brain_cre, f = 0.5, e = TRUE, wa = TRUE, wb = TRUE)
    kl_k_cre <- bt.intersect(a = kidney_cre, b = liver_cre, f = 0.5, e = TRUE, wa = TRUE, wb = TRUE)
    kl_l_cre <- bt.intersect(a = liver_cre, b = kidney_cre, f = 0.5, e = TRUE, wa = TRUE, wb = TRUE)
    bk_b_cre$ep <- alply(bk_b_cre, 1, judge_ep, .parallel = TRUE) %>% unlist()
    bk_k_cre$ep <- alply(bk_k_cre, 1, judge_ep, .parallel = TRUE) %>% unlist()
    bl_b_cre$ep <- alply(bl_b_cre, 1, judge_ep, .parallel = TRUE) %>% unlist()
    bl_l_cre$ep <- alply(bl_l_cre, 1, judge_ep, .parallel = TRUE) %>% unlist()
    kl_k_cre$ep <- alply(kl_k_cre, 1, judge_ep, .parallel = TRUE) %>% unlist()
    kl_l_cre$ep <- alply(kl_l_cre, 1, judge_ep, .parallel = TRUE) %>% unlist()
    unique(kl_k_cre[, c("V1", "V2", "V3", "V4", "ep")]) %>%
        .[.$ep == "ep", ] %>%
        nrow()
    ep_cre <- rbind(bk_b_cre, bk_k_cre, bl_b_cre, bl_l_cre, kl_k_cre, kl_l_cre)[, c("V1", "V2", "V3", "V4", "ep")] %>%
        unique() %>%
        .[.$ep == "ep", ]
    cre_info$duality <- "no"
    cre_info[cre_info$unique_id %in% ep_cre$V4, "duality"] <- "yes"
    nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "enhancer"), ]) / nrow(cre_info[which(cre_info$element == "enhancer"), ])
    nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "promoter"), ]) / nrow(cre_info[which(cre_info$element == "promoter"), ])
    write.csv(cre_info, file = paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_duality.csv"), row.names = FALSE)
    res <- mdply(expand.grid(tissue = c("Brain", "Kidney", "Liver"), element = c("enhancer", "promoter")), function(tissue, element) {
        rr <- nrow(cre_info[which(cre_info$duality == "yes" & cre_info$tissue == tissue & cre_info$element == element), ]) / nrow(cre_info[which(cre_info$tissue == tissue & cre_info$element == element), ])
        count <- nrow(cre_info[which(cre_info$duality == "yes" & cre_info$tissue == tissue & cre_info$element == element), ])
        return(data.frame(species = sp, tissue = tissue, element = element, duality_count = count, duality_prop = rr))
    })
    return(res)
}, .parallel = TRUE)
write.csv(duality_df, file = "/media/Data/zhangz/chip/analysis/summary2/abc_1M/duality.csv", row.names = FALSE)
summary(duality_df)

tissue_color <- c("Brain" = "#ffe0c2", "Kidney" = "#a9d9b2", "Liver" = "#a5aecb")
## plot the duality of enhancer and promoter
# duality_df <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/duality.csv", header = TRUE, sep = ",")
duality_p <- ggplot(duality_df, aes(x = species, y = duality_prop, fill = tissue)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(element ~ .) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = scales::percent, n.breaks = 3) +
    scale_fill_manual(values = tissue_color) +
    labs(x = "Species", y = "Duality proportion", fill = "Tissue") +
    ggtitle("Duality of enhancer and promoter in different tissues") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 8)) +
    theme(panel.grid.major = element_line(linewidth = 0.1, linetype = "dashed"), panel.grid.minor = element_blank())
ggsave("/media/Data/zhangz/chip/analysis/summary2/abc_1M/duality.pdf", duality_p, width = 10, height = 6)

# whether dualistic cres are ancestral or species-specific
sp <- "Mus_musculus"
cre_info <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_duality_specificity.csv"), header = TRUE, sep = ",")

cre_info$specificity <- "N"
cre_info$specificity[which(cre_info$act_time > 0 &
    (cre_info$s1_1 > 0 | cre_info$s0_1 > 0 | !cre_info$ele_pattern %in% c("B", "K", "L")) &
    cre_info$element == "enhancer")] <- "evol_pleio_e"
cre_info$specificity[which(cre_info$act_time > 0 &
    (cre_info$s1_1 > 0 | cre_info$s0_1 > 0 | !cre_info$ele_pattern %in% c("B", "K", "L")) &
    cre_info$element == "promoter")] <- "evol_pleio_p"
cre_info$specificity[which(cre_info$act_time > 0 &
    (cre_info$s1_1 == 0 & cre_info$s0_1 == 0 & cre_info$ele_pattern %in% c("B", "K", "L")) &
    cre_info$element == "enhancer")] <- "evol_ts_e"
cre_info$specificity[which(cre_info$act_time > 0 &
    (cre_info$s1_1 == 0 & cre_info$s0_1 == 0 & cre_info$ele_pattern %in% c("B", "K", "L")) &
    cre_info$element == "promoter")] <- "evol_ts_p"
cre_info$specificity[which(cre_info$act_time == 0 &
    cre_info$ele_pattern %in% c("B", "K", "L") &
    cre_info$element == "enhancer")] <- "sp_ts_e"
cre_info$specificity[which(cre_info$act_time == 0 &
    cre_info$ele_pattern %in% c("B", "K", "L") &
    cre_info$element == "promoter")] <- "sp_ts_p"
cre_info$specificity[which(cre_info$act_time == 0 &
    !cre_info$ele_pattern %in% c("B", "K", "L") &
    cre_info$element == "enhancer")] <- "sp_pleio_e"
cre_info$specificity[which(cre_info$act_time == 0 &
    !cre_info$ele_pattern %in% c("B", "K", "L") &
    cre_info$element == "promoter")] <- "sp_pleio_p"
nrow(cre_info[which(cre_info$duality == "yes" & cre_info$specificity == "evol_pleio_e"), ]) / nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "enhancer"), ])
nrow(cre_info[which(cre_info$duality == "yes" & cre_info$specificity == "evol_ts_e"), ]) / nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "enhancer"), ])
nrow(cre_info[which(cre_info$duality == "yes" & cre_info$specificity == "sp_ts_e"), ]) / nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "enhancer"), ])
nrow(cre_info[which(cre_info$duality == "yes" & cre_info$specificity == "sp_pleio_e"), ]) / nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "enhancer"), ])
nrow(cre_info[which(cre_info$duality == "yes" & cre_info$specificity == "evol_pleio_p"), ]) / nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "promoter"), ])
nrow(cre_info[which(cre_info$duality == "yes" & cre_info$specificity == "evol_ts_p"), ]) / nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "promoter"), ])
nrow(cre_info[which(cre_info$duality == "yes" & cre_info$specificity == "sp_ts_p"), ]) / nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "promoter"), ])
nrow(cre_info[which(cre_info$duality == "yes" & cre_info$specificity == "sp_pleio_p"), ]) / nrow(cre_info[which(cre_info$duality == "yes" & cre_info$element == "promoter"), ])

count(cre_info[which(cre_info$duality == "yes" & cre_info$element == "enhancer"), ], "specificity")
count(cre_info[which(cre_info$duality == "yes" & cre_info$element == "promoter"), ], "specificity")

write.table(cre_info[, c("chr", "start", "end", "unique_id", "duality", "specificity")],
    file = paste0("/media/Data/zhangz/chip/analysis/", sp, "/anno/blast/", sp, "_overlap_duality_specificity.bed"),
    row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
)

blast_out <- read.csv("/media/Data/zhangz/chip/analysis/Mus_musculus/anno/blast/test_blast.out", header = FALSE, sep = "\t")
blast_out <- blast_out[-1, ]
colnames(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast_out <- blast_out[which(blast_out$qseqid != blast_out$sseqid), ]

cre_info[which(cre_info$unique_id %in% c("enh_b_100072_B_Peak_84147", "enh_b_100073_B_Peak_84149")), c("chr", "start", "end", "unique_id", "duality", "specificity")]
cre_info[which(cre_info$unique_id %in% c("enh_kl_9690_L_Peak_291", "enh_k_179322_K_Peak_83331", "enh_kl_3378_L_Peak_23391")), c("chr", "start", "end", "unique_id", "duality", "specificity")]

blast_out <- adply(blast_out, 1, function(row) {
    row$qseq <- unlist(strsplit(row$qseqid, "_"))[1:3] %>% paste(collapse = "_")
    row$sseq <- unlist(strsplit(row$sseqid, "_"))[1:3] %>% paste(collapse = "_")
    return(row)
}, .parallel = TRUE)
dim(blast_out)
# [1] 674612     14
blast_out <- blast_out[which(blast_out$qseq != blast_out$sseq), ]
dim(blast_out)
# [1] 125186     14
nrow(blast_out[blast_out$evalue < 1e-5, ])
blast_out <- blast_out[blast_out$evalue < 1e-5, ]
summary(blast_out)
length(unique(blast_out$qseq))
length(unique(blast_out$sseq))
length(intersect(unique(blast_out$qseq), unique(blast_out$sseq)))
length(setdiff(unique(blast_out$qseq), unique(blast_out$sseq)))
length(unique(blast_out$qseqid))
length(unique(blast_out$sseqid))
cre_info$blast_hit <- 0
cre_info$blast_hit <- alply(cre_info, 1, function(row) {
    if (row$unique_id %in% blast_out$qseqid) {
        res <- length(unique(blast_out[which(blast_out$qseqid == row$unique_id), "sseq"]))
    } else {
        res <- 0
    }
    return(res)
}, .parallel = TRUE) %>% unlist()

# r$> summary(cre_info$blast_hit)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.000   0.000   0.000   1.052   2.000   4.000

# r$> summary(cre_info$blast_hit[cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.000   0.000   0.000   1.319   3.000   4.000

# r$> summary(cre_info$blast_hit[!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.0000  0.0000  0.8018  1.0000  4.0000

wilcox.test(
    cre_info$blast_hit[cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")],
    cre_info$blast_hit[!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")]
)
t.test(
    cre_info$blast_hit[cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")],
    cre_info$blast_hit[!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")]
)
nrow(cre_info[which(cre_info$blast_hit > 0), ]) / nrow(cre_info)

blast_out <- read.csv("/media/Data/zhangz/chip/analysis/Mus_musculus/anno/blast/test_blast.out", header = FALSE, sep = "\t")
## 现在已经移除了query和subject是相同序列的行，但是enhancer与promoter的序列重叠的情况还没有考虑，要想办法找出这些有duality的序列
## 将query与subject是同一段序列但是分别是enhancer与promoter，从而在之前的记录中没有表明的情况找出来去除掉


all_cre <- read.csv("/media/Data/zhangz/chip/analysis/summary2/anno_sum/Mus_musculus/Mus_musculus_all.csv", header = TRUE, sep = ",")
colnames(all_cre) <- all_cre[1, ]
all_cre <- all_cre[-1, ]
all_cre$unique_id <- paste(all_cre$unique_id, all_cre$peak, sep = "_")
write.table(all_cre[, c("chr", "start", "end", "unique_id")],
    file = "/media/Data/zhangz/chip/analysis/Mus_musculus/anno/blast2/Mus_musculus_all.bed",
    row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE
)
# bedtools getfasta -nameOnly -fi /media/Data/zhangz/chip/genomes/Mus_musculus/Mus_musculus.fa -bed Mus_musculus_all.bed -fo fastadb/Mus_musculus_all_cre.fasta
# makeblastdb -in fastadb/Mus_musculus_all_cre.fasta -dbtype nucl -out fastadb/Mus_musculus_all_cre
# nohup blastn -query fastadb/Mus_musculus_all_cre.fasta -db fastadb/Mus_musculus_all_cre -num_threads 48 -evalue 1e-5 -outfmt 6 -out /media/Data/zhangz/chip/analysis/Mus_musculus/anno/blast2/test_blast.out > /media/Data/zhangz/chip/analysis/Mus_musculus/anno/blast2/test_blast.log 2>&1 &

## blast output of all cres
blast_out <- read.csv("/media/Data/zhangz/chip/analysis/Mus_musculus/anno/blast2/test_blast_with_qseqsseq.out", header = FALSE, sep = "\t")

# blast_out <- blast_out[-1, ]
colnames(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq")
blast_out <- blast_out[which(blast_out$qseqid != blast_out$sseqid), ]
nrow(blast_out)
# [1] 24035259

cre_info[which(cre_info$unique_id %in% c("enh_b_100072_B_Peak_84147", "enh_b_100073_B_Peak_84149")), c("chr", "start", "end", "unique_id", "duality", "specificity")]
cre_info[which(cre_info$unique_id %in% c("enh_kl_9690_L_Peak_291", "enh_k_179322_K_Peak_83331", "enh_kl_3378_L_Peak_23391")), c("chr", "start", "end", "unique_id", "duality", "specificity")]

# for (i in 1:nrow(blast_out)) {
#     blast_out[i, "qseq"] <- unlist(strsplit(blast_out[i, "qseqid"], "_"))[1:3] %>% paste(collapse = "_")
#     blast_out[i, "sseq"] <- unlist(strsplit(blast_out[i, "sseqid"], "_"))[1:3] %>% paste(collapse = "_")
# }
# blast_out <- adply(blast_out, 1, function(row) {
#     row$qseq <- unlist(strsplit(row$qseqid, "_"))[1:3] %>% paste(collapse = "_")
#     row$sseq <- unlist(strsplit(row$sseqid, "_"))[1:3] %>% paste(collapse = "_")
#     return(row)
# }, .parallel = TRUE)
dim(blast_out)
# [1] 24502383     14
blast_out <- blast_out[which(blast_out$qseq != blast_out$sseq), ]
dim(blast_out)
# [1] 23934945     14
nrow(blast_out[blast_out$evalue < 1e-5, ])
blast_out <- blast_out[blast_out$evalue < 1e-5, ]
# 23934870
summary(blast_out)
nrow(unique(blast_out[, c("qseq", "sseq")]))
# [1] 19246227
length(unique(blast_out$qseq))
length(unique(blast_out$sseq))
length(intersect(unique(blast_out$qseq), unique(blast_out$sseq)))
length(setdiff(unique(blast_out$qseq), unique(blast_out$sseq)))
length(unique(blast_out$qseqid))
length(unique(blast_out$sseqid))
cre_info$blast_hit <- 0
cre_info$blast_hit <- alply(cre_info, 1, function(row) {
    if (row$unique_id %in% blast_out$qseqid) {
        res <- length(unique(blast_out[which(blast_out$qseqid == row$unique_id), "sseq"]))
    } else {
        res <- 0
    }
    return(res)
}, .parallel = TRUE) %>% unlist()
for (i in seq_len(nrow(cre_info))) {
    if (cre_info[i, "unique_id"] %in% blast_out$qseqid) {
        cre_info[i, "blast_hit"] <- length(unique(blast_out[which(blast_out$qseqid == cre_info[i, "unique_id"]), "sseq"]))
    }
}

summary(cre_info$blast_hit)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0.00    0.00    1.00   62.61    4.00  497.00
summary(cre_info$blast_hit[cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0.00    0.00    1.00   83.42   21.00  497.00
plot(density(cre_info$blast_hit[cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")]), col = "red", lwd = 2, main = "Density of blast hit count", xlab = "Blast hit count")
# how many new cres have blast hits, that is, may be related with ancestral elements
nrow(cre_info[which(cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_info$blast_hit > 0), ])
# the ratio of new cres have blast hits
nrow(cre_info[which(cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_info$blast_hit > 0), ]) / nrow(cre_info[which(cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")), ])
# the ratio of new cres related with repeat have blast hits
nrow(cre_info[which(cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") &
    cre_info$repeat_id != "" &
    cre_info$blast_hit > 0), ]) /
    nrow(cre_info[which(cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") &
        cre_info$repeat_id != ""), ])
# the ratio of new cres not related with repeat have blast hits
nrow(cre_info[which(cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") &
    cre_info$repeat_id == "" &
    cre_info$blast_hit > 0), ]) /
    nrow(cre_info[which(cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") &
        cre_info$repeat_id == ""), ])

summary(cre_info$blast_hit[!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")])

#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0.00    0.00    1.00   43.09    2.00  497.00
plot(density(cre_info$blast_hit[!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")]), col = "blue", lwd = 2, add = TRUE)
# how many old cres have blast hits, that is, may be related with new elements
nrow(cre_info[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_info$blast_hit > 0), ])
# the ratio of old cres have blast hits
nrow(cre_info[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_info$blast_hit > 0), ]) /
    nrow(cre_info[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")), ])
wilcox.test(
    cre_info$blast_hit[cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")],
    cre_info$blast_hit[!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")]
)
# the ratio of old cres related with repeat have blast hits
nrow(cre_info[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") &
    cre_info$repeat_id != "" &
    cre_info$blast_hit > 0), ]) /
    nrow(cre_info[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") &
        cre_info$repeat_id != ""), ])
# the ratio of old cres not related with repeat have blast hits
nrow(cre_info[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") &
    cre_info$repeat_id == "" &
    cre_info$blast_hit > 0), ]) /
    nrow(cre_info[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") &
        cre_info$repeat_id == ""), ])

cre_blast <- adply(cre_info, 1, function(row) {
    if (row$unique_id %in% blast_out$qseqid) {
        blast_tmp <- blast_out[which(blast_out$qseqid == row$unique_id), ]
        # whether the blast hits including ancestral elements
        if (any(blast_tmp$sseqid %in% cre_info$unique_id[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p"))])) {
            row$blast_hit_anc <- 1
            # whether the ancestral elements are related with repeat elements
            anc_source <- cre_info[which(cre_info$unique_id %in% blast_tmp$sseqid[which(blast_tmp$sseqid %in% cre_info$unique_id[which(!cre_info$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p"))])]), ]
            if (any(anc_source$repeat_id != "")) {
                row$blast_hit_anc_repeat <- 1
            } else {
                row$blast_hit_anc_repeat <- 0
            }
            # whether the  most ancestral element is pleiotropic
            if (any(anc_source$specificity %in% c("evol_pleio_e", "evol_pleio_p"))) {
                row$blast_hit_anc_pleio <- 1
            } else {
                row$blast_hit_anc_pleio <- 0
            }
        } else {
            row$blast_hit_anc <- 0
            row$blast_hit_anc_repeat <- 0
            row$blast_hit_anc_pleio <- 0
        }
    } else {
        row$blast_hit_anc <- 0
        row$blast_hit_anc_repeat <- 0
        row$blast_hit_anc_pleio <- 0
    }
    return(row)
}, .parallel = TRUE)

write.csv(cre_blast, file = "/media/Data/zhangz/chip/analysis/summary2/abc_1M/Mus_musculus_overlap_duality_specificity_blast.csv", row.names = FALSE)

cre_blast <- read.csv("/media/Data/zhangz/chip/analysis/summary2/abc_1M/Mus_musculus_overlap_duality_specificity_blast.csv", header = TRUE, sep = ",")
# how many cres have blast hits
nrow(cre_blast[which(cre_blast$blast_hit > 0), ]) / nrow(cre_blast)
# [1] 0.5653524
# how many cres with blast hit are species specific
nrow(cre_blast[which(cre_blast$blast_hit > 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")), ]) / nrow(cre_blast[which(cre_blast$blast_hit > 0), ])
# [1] 0.5219501
# how many species-specific cres have blast hits
nrow(cre_blast[which(cre_blast$blast_hit > 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")), ]) / nrow(cre_blast[which(cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")), ])
# [1] 0.6099093
# how many cres with blast hit are species specific and have ancestral elements
nrow(cre_blast[which(cre_blast$blast_hit > 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_blast$blast_hit_anc == 1), ]) / nrow(cre_blast[which(cre_blast$blast_hit > 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")), ])
# [1] 0.5923703
# how many cres with blast hit are species specific and have ancestral elements and the ancestral elements are not related with repeat elements
nrow(cre_blast[which(cre_blast$blast_hit > 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_blast$blast_hit_anc == 1 & cre_blast$blast_hit_anc_repeat == 0), ]) / nrow(cre_blast[which(cre_blast$blast_hit > 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_blast$blast_hit > 0), ])
# [1] 0.0555448
# how many cres with blast hit are species specific and have ancestral elements and the ancestral elements are related with repeat elements
nrow(cre_blast[which(cre_blast$blast_hit > 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_blast$blast_hit_anc == 1 & cre_blast$blast_hit_anc_repeat == 1), ]) / nrow(cre_blast[which(cre_blast$blast_hit > 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p") & cre_blast$blast_hit > 0), ])
# [1] 0.5368255
# how many young cres without blast hit are not related with repeat elements
nrow(cre_blast[which(cre_blast$blast_hit == 0 & cre_blast$repeat_id == "" & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")), ]) / nrow(cre_blast[which(cre_blast$blast_hit == 0 & cre_blast$specificity %in% c("sp_ts_e", "sp_ts_p", "sp_pleio_e", "sp_pleio_p")), ])
# [1] 0.5134985
