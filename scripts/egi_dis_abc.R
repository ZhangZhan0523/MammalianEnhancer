## calculate the conservation of egi
## egi were based on abc and 1Mts annotation


library(optparse)
option_list <- list(
    make_option(
        c("-s", "--species"),
        type = "character", default = "/media/Data/zhangz/chip/scripts/info/info_using.csv",
        help = "species need to summarize"
    ),
    make_option(
        c("-t", "--threshold"),
        type = "numeric", default = 10,
        help = "threshold of peak conservation egi"
    )
)
opt_parser <- parse_args(OptionParser(option_list = option_list, usage = "%prog [options]"))
sp <- opt_parser$species
td <- opt_parser$threshold
print(sp)
print(td)
library(magrittr)
library(stringr)
library(plyr)
library(doMC)
registerDoMC(4)


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
sp_target <- species[species != sp]

egi_df <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_egi_abc_new.csv"),
    header = TRUE
)
high_ec_egi <- egi_df[egi_df$consrv >= td, ]

# egi_con_df <- data.frame()
# for (sp_t in sp_target) {
egi_con_df <- adply(sp_target, 1, function(sp_t) {
    target_egi <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp_t, "_egi_abc_new.csv"),
        header = TRUE
    )
    tmp_egi_con_df <- data.frame()
    for (tis in intersect(unique(high_ec_egi$tissue), unique(target_egi$tissue))) {
        for (ele in c("enhancer", "promoter")) {
            print(paste(sp_t, tis, ele, sep = ":"))
            ec <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "_", tis, "_", ele, "_ec_group.csv"), header = TRUE)
            high_ec <- ec[ec$ec_id %in% unique(high_ec_egi$ec_id), ]
            # for (peak in unique(high_ec_egi$peak[high_ec_egi$tissue == tis & high_ec_egi$element == ele])) {
            calc_same_ratio <- function(peak) {
                # print(peak)
                # print(length(peak))
                if (is.null(high_ec[high_ec[1] == peak, sp_t]) | length(high_ec[high_ec[1] == peak, sp_t]) == 0) {
                    return(data.frame())
                }
                target_peak <-
                    str_split(high_ec[high_ec[1] == peak, sp_t], ",")[[1]]
                target_egi_t <- target_egi[which(target_egi$peak %in% target_peak & target_egi$tissue == tis & target_egi$element == ele), ]
                if (nrow(target_egi_t) == 0) {
                    return(data.frame())
                }
                source_egi <- high_ec_egi[which(high_ec_egi$peak == peak & high_ec_egi$tissue == tis & high_ec_egi$element == ele), ]
                og_s <- unique(source_egi$gene)
                og_t <- unique(target_egi_t$gene)
                ec_id <- unique(source_egi$ec_id)
                unique_id <- unique(source_egi$unique_id)
                unique_id_t <- unique(target_egi_t$unique_id)
                same_og <- intersect(unique(source_egi$og), unique(target_egi_t$og))
                if (length(same_og) == 0) {
                    same_og <- NA
                }
                # same_og <- ifelse(length(same_og) == 0, NA, same_og) # only return the first og
                # same_num <- ifelse(is.na(same_og), 0, length(same_og))
                if (TRUE %in% is.na(same_og)) {
                    same_num <- 0
                } else {
                    same_num <- length(same_og)
                }
                same_ratio <- ifelse(nrow(target_egi_t) != 0,
                    (length(same_og) / sqrt(length(og_s) * length(og_t))), NA
                )
                distance_s <- ifelse(is.na(same_og), NA, source_egi$distance[source_egi$og %in% same_og])
                distance_t <- ifelse(is.na(same_og), NA, target_egi_t$distance[target_egi_t$og %in% same_og])
                distance_t <- ifelse(same_num == 1, unique(distance_t), distance_t)
                distance_diff <- ifelse(is.na(same_og), NA, distance_s - distance_t)
                note_abc <- unique(source_egi$note_abc)[1]
                note_abc_t <- unique(target_egi_t$note_abc)[1]
                ne_same <- ifelse(unique(source_egi$nearest_og) %in% same_og, 1, 0)
                tmp <- data.frame(
                    ec_id = ec_id,
                    peak_s = peak, sp_s = sp, peak_t = paste(unique(target_egi_t$peak), collapse = ";"), sp_t = sp_t,
                    og_s = length(og_s), og_t = length(og_t),
                    same_og = same_num,
                    same_og_name = paste(same_og, collapse = ";"),
                    same_og_ratio = same_ratio,
                    unique_id_s = unique_id, unique_id_t = paste(unique_id_t, collapse = ";"),
                    distance_s = paste(distance_s, collapse = ";"), distance_t = paste(distance_t, collapse = ";"),
                    distance_diff = paste(distance_diff, collapse = ";"),
                    note_abc_s = note_abc, note_abc_t = note_abc_t,
                    align_s = unique(source_egi$align), align_t = paste(unique(target_egi_t$align[target_egi_t$og %in% same_og]), collapse = ";"),
                    consrv_s = unique(source_egi$consrv), consrv_t = paste(unique(target_egi_t$consrv[target_egi_t$og %in% same_og]), collapse = ";"),
                    tissue = tis, element = ele,
                    ne_same = ne_same
                )
                # egi_con_df <- rbind(egi_con_df, tmp)
                return(tmp)
            }
            tmp_egi_con_df <- rbind(tmp_egi_con_df, adply(unique(high_ec_egi$peak[high_ec_egi$tissue == tis & high_ec_egi$element == ele]), 1, calc_same_ratio, .parallel = TRUE))
        }
    }
    return(tmp_egi_con_df)
})
save_cols <- c("ec_id", "peak_s", "sp_s", "peak_t", "sp_t", "og_s", "og_t", "same_og", "same_og_name", "same_og_ratio", "unique_id_s", "unique_id_t", "distance_s", "distance_t", "distance_diff", "note_abc_s", "note_abc_t", "align_s", "align_t", "consrv_s", "consrv_t", "tissue", "element", "ne_same")
write.csv(egi_con_df[, save_cols], paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_abc_egi_con_t", td, ".csv"), row.names = FALSE)
egi_con_df <- egi_con_df[, save_cols]
stat_unique_id <- function(id) {
    tmp <- egi_con_df[egi_con_df$unique_id_s == id, ]
    same_ogs <- unique(tmp$same_og_name[!is.na(tmp$same_og_name)])
    same_og_freq <- table(tmp$same_og_name[!is.na(tmp$same_og_name)])
    most_same_og <- names(same_og_freq)[which.max(same_og_freq)]
    most_same_ratio <- same_og_freq[which.max(same_og_freq)] / sum(same_og_freq)
    same_og_num <- length(same_ogs)
    same_sp <- unique(tmp$sp_t[tmp$same_og != 0])
    same_sp_num <- length(same_sp)
    diff_sp <- unique(tmp$sp_t[tmp$same_og == 0])
    diff_sp_num <- length(diff_sp)
    distance_diff <- paste(tmp$sp_t, tmp$distance_diff[!is.na(tmp$same_og_name)], sep = ":")
    same_og_pair <- paste(tmp$sp_t, tmp$same_og_name, sep = ":")
    res <- data.frame(
        unique_id = id,
        peak = paste(unique(tmp$peak_s), collapse = ";"),
        sp = unique(tmp$sp_s),
        tissue = unique(tmp$tissue),
        element = unique(tmp$element),
        align = paste(unique(tmp$align_s), collapse = ";"),
        consrv = unique(tmp$consrv_s),
        distance = paste(unique(tmp$distance_s[!is.na(tmp$same_og_name)]), collapse = ";"),
        distance_diff = paste(distance_diff, collapse = ";"),
        same_og_pair = paste(same_og_pair, collapse = ";"),
        same_og_num = same_og_num,
        same_ogs = paste(same_ogs, collapse = ";"),
        most_same_og = most_same_og,
        most_same_ratio = most_same_ratio,
        same_sp_num = same_sp_num,
        same_sp = paste(same_sp, collapse = ";"),
        diff_sp_num = diff_sp_num,
        diff_sp = paste(diff_sp, collapse = ";")
    )
    return(res)
}
unique_id_df <- adply(unique(egi_con_df$unique_id_s), 1, stat_unique_id, .parallel = TRUE)
write.csv(unique_id_df, paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_abc_egi_con_t", td, "_unique_id.csv"), row.names = FALSE)
