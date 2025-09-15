## do blast of all cre sequences
## count the ratio of new cre sequences come from old cre sequences by repeat

# load packages
library(plyr)
library(doMC)
registerDoMC(cores = 4)
# packages for excecuting bash command
library(systemPipeR)
library(magrittr)

## first, find the sequences of all cre sequences, build a blast database, then blast all cre sequences against the database
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
sp <- "Mus_musculus"

a_ply(species, 1, function(sp) {
    if (sp == "Macaca_mulatta") {
        sp_g <- "Macaca_fascicularis"
    } else {
        sp_g <- sp
    }
    all_cre <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/anno_sum/", sp, "/", sp, "_all.csv"), header = TRUE, sep = ",")
    all_cre$unique_id <- paste(all_cre$unique_id, all_cre$peak, sep = "_")
    # mkdir anno/blast2
    system(paste0("mkdir -p /media/Data/zhangz/chip/analysis/", sp, "/anno/blast2"))
    all_bed <- paste0("/media/Data/zhangz/chip/analysis/", sp, "/anno/blast2/", sp, "_all.bed")
    write.table(all_cre[, c("chr", "start", "end", "unique_id")],
        file = all_bed,
        row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE
    )

    all_fa <- paste0("/media/Data/zhangz/chip/analysis/", sp, "/anno/blast2/fastadb/", sp, "_all_cre.fa")
    # mkdir fastadb
    system(paste0("mkdir -p /media/Data/zhangz/chip/analysis/", sp, "/anno/blast2/fastadb"))

    # draw sequence fasta file using bedtools
    system(paste0("bedtools getfasta -nameOnly -fi /media/Data/zhangz/chip/genomes/", sp_g, "/", sp_g, ".fa -bed ", all_bed, " -fo ", all_fa))

    db_name <- paste0("/media/Data/zhangz/chip/analysis/", sp, "/anno/blast2/fastadb/", sp, "_all_cre")
    # build blast database
    system(paste0("conda run -n R makeblastdb -in ", all_fa, " -dbtype nucl -out ", db_name))

    # blast all cre sequences against the database
    blast_out <- paste0("/media/Data/zhangz/chip/analysis/", sp, "/anno/blast2/", sp, "_blast.out")
    system(paste0("conda run -n R blastn -query ", all_fa, " -db ", db_name, " -out ", blast_out, " -evalue 1e-5 -outfmt 6 -num_threads 12"))

    blast_res <- read.csv(blast_out, header = FALSE, sep = "\t")
    colnames(blast_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    blast_res <- blast_res[which(blast_res$qseqid != blast_res$sseqid), ]
    blast_res <- adply(blast_res, 1, function(row) {
        row$qseq <- unlist(strsplit(row$qseqid, "_"))[1:3] %>% paste(collapse = "_")
        row$sseq <- unlist(strsplit(row$sseqid, "_"))[1:3] %>% paste(collapse = "_")
        return(row)
    })

    cre_info <- read.csv(paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_duality.csv"), header = TRUE, sep = ",")

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

    cre_info$blast_hit <- 0
    cre_info$blast_hit <- alply(cre_info, 1, function(row) {
        if (row$unique_id %in% blast_res$qseqid) {
            res <- length(unique(blast_res[which(blast_res$qseqid == row$unique_id), "sseq"]))
        } else {
            res <- 0
        }
        return(res)
    }, .parallel = TRUE) %>% unlist()
    for (i in seq_len(nrow(cre_info))) {
        if (cre_info[i, "unique_id"] %in% blast_res$qseqid) {
            cre_info[i, "blast_hit"] <- length(unique(blast_res[which(blast_res$qseqid == cre_info[i, "unique_id"]), "sseq"]))
        }
    }

    cre_blast <- adply(cre_info, 1, function(row) {
        if (row$unique_id %in% blast_res$qseqid) {
            blast_tmp <- blast_res[which(blast_res$qseqid == row$unique_id), ]
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

    cre_blast <- cre_blast[, c(
        "unique_id", "peak", "tissue", "element", "chr", "start", "end", "species", "act_time", "duality",
        "ele_pattern", "repeat_id", "specificity", "blast_hit", "blast_hit_anc", "blast_hit_anc_repeat", "blast_hit_anc_pleio"
    )]
    summary(cre_blast)
    write.csv(cre_blast, file = paste0("/media/Data/zhangz/chip/analysis/summary2/abc_1M/", sp, "_overlap_duality_specificity_blast.csv"), row.names = FALSE)
}, .parallel = TRUE)


cre_info[which(cre_info$chr == "chr15_rhiPus" & cre_info$start >= 61587000 & cre_info$end <= 61591317), c("unique_id", "chr", "start", "end")]
