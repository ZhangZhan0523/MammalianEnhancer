## to make conservation plots of two categories of reguleroty elements:
## enhancer and promoter,
## with the same means of cell paper

## input: alignable and conservation values of different species, tissues and
##        different categories of regulatory elements,like:
##        /media/Data/zhangz/chip/analysis/Mus_musculus/compare2/Brain/enhancer/\
##        Mus_musculus_Brain_enhancer_stats.tsv.con

library(ggplot2)
library(stringr)
library(ape)
library(ggtree)
library(PerformanceAnalytics)
library(reshape2)
library(car)
library(classInt)

setwd("/media/Data/zhangz/chip/analysis/conservation3")
## test on single species
con_table <- read.table(
    "/media/Data/zhangz/chip/analysis/Hipposideros_larvatus/compare2/Brain/promoter/Hipposideros_larvatus_Brain_promoter_stats.tsv.con",
    header = TRUE, sep = "\t"
)
plot(con_table$consrv, con_table$alignable, xlab = "conservation", ylab = "alignable")
## read in the data
enhancer_list <- read.table(
    "/media/Data/zhangz/chip/analysis/summary2/conservation3/enhancer_list.csv"
)
promoter_list <- read.table(
    "/media/Data/zhangz/chip/analysis/summary2/conservation3/promoter_list.csv"
)
enhancer_con <- data.frame()
promoter_con <- data.frame()
peak_count <- data.frame(matrix(ncol = 4))
colnames(peak_count) <- c("sp", "tis", "ele", "count")
for (i in enhancer_list$V1) {
    new_df <- read.csv(i, sep = "\t", header = TRUE, row.names = 1)
    if (!("Neophocaena_asiaeorientalis_con" %in% colnames(new_df))) {
        new_df$Neophocaena_asiaeorientalis_con <- 0
    }
    if (!("Rhinopithecus_roxellana_con" %in% colnames(new_df))) {
        new_df$Rhinopithecus_roxellana_con <- 0
    }
    condi <- strsplit(gsub(".tsv.con", "", basename(i)), "_")[[1]]
    sp <- paste(condi[1], condi[2], sep = "_")
    tis <- condi[3]
    ele <- condi[4]
    eval(parse(text = paste("new_df$", sp, "_con = 0", sep = "")))
    eval(parse(text = paste("new_df$", sp, " = 0", sep = "")))
    peak_count <- rbind(peak_count, sp, tis, ele, nrow(new_df))
    rownames(new_df) <- paste(sub("_stats.tsv.con", "", basename(i)), rownames(new_df), sep = "_")
    enhancer_con <- rbind(enhancer_con, new_df)
}
for (i in promoter_list$V1) {
    new_df <- read.csv(i, sep = "\t", header = TRUE, row.names = 1)
    if (!("Neophocaena_asiaeorientalis_con" %in% colnames(new_df))) {
        new_df$Neophocaena_asiaeorientalis_con <- 0
    }
    if (!("Rhinopithecus_roxellana_con" %in% colnames(new_df))) {
        new_df$Rhinopithecus_roxellana_con <- 0
    }
    condi <- strsplit(gsub(".tsv.con", "", basename(i)), "_")[[1]]
    sp <- paste(condi[1], condi[2], sep = "_")
    tis <- condi[3]
    ele <- condi[4]
    eval(parse(text = paste("new_df$", sp, "_con = 0", sep = "")))
    eval(parse(text = paste("new_df$", sp, " = 0", sep = "")))
    peak_count <- rbind(peak_count, c(paste(condi[1], condi[2], "_"), condi[3], condi[4], nrow(new_df)))
    rownames(new_df) <- paste(sub("_stats.tsv.con", "", basename(i)), rownames(new_df), sep = "_")
    promoter_con <- rbind(promoter_con, new_df)
}
enhancer_dist_no0 <- ggplot(enhancer_con[enhancer_con$align != 0, ], aes(x = align)) +
    geom_histogram(binwidth = 1, fill = "#ffb677", colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("enhancer_distribution") +
    xlab("align") +
    ylab("count") +
    theme(text = element_text(size = 45))
png("enhancer_dist_no0.png", width = 8000, height = 6000, res = 600)
plot(enhancer_dist_no0)
dev.off()
promoter_dist_no0 <- ggplot(promoter_con[promoter_con$align != 0, ], aes(x = align)) +
    geom_histogram(binwidth = 1, fill = "#a893c1", colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("promoter_distribution") +
    xlab("align") +
    ylab("count") +
    theme(text = element_text(size = 45))
png("promoter_dist_no0.png", width = 8000, height = 6000, res = 600)
plot(promoter_dist_no0)
dev.off()
mus_enhancer_con <- enhancer_con[grep("Mus_musculus", rownames(enhancer_con)), ]
mus_promoter_con <- promoter_con[grep("Mus_musculus", rownames(promoter_con)), ]
can_enhancer_con <- enhancer_con[grep("Canis_lupus", rownames(enhancer_con)), ]
can_promoter_con <- promoter_con[grep("Canis_lupus", rownames(promoter_con)), ]
can_enhancer_dist_no0 <- ggplot(can_enhancer_con[can_enhancer_con$align != 0, ], aes(x = align)) +
    geom_histogram(binwidth = 1, fill = "#ffb677", colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Canis lupus enhancer distribution") +
    xlab("align") +
    ylab("count") +
    theme(text = element_text(size = 40))
png("can_enhancer_dist_no0.png", width = 8000, height = 6000, res = 600)
plot(can_enhancer_dist_no0)
dev.off()
can_promoter_dist_no0 <- ggplot(can_promoter_con[can_promoter_con$align != 0, ], aes(x = align)) +
    geom_histogram(binwidth = 1, fill = "#a893c1", colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Canis lupus promoter distribution") +
    xlab("align") +
    ylab("count") +
    theme(text = element_text(size = 40))
png("can_promoter_dist_no0.png", width = 8000, height = 6000, res = 600)
plot(can_promoter_dist_no0)
dev.off()
mus_enhancer_dist_no0 <- ggplot(mus_enhancer_con[mus_enhancer_con$align != 0, ], aes(x = align)) +
    geom_histogram(binwidth = 1, fill = "#ffb677", colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Mus musculus enhancer distribution") +
    xlab("align") +
    ylab("count") +
    theme(text = element_text(size = 40))
png("mus_enhancer_dist_no0.png", width = 8000, height = 6000, res = 600)
plot(mus_enhancer_dist_no0)
dev.off()
mus_promoter_dist_no0 <- ggplot(mus_promoter_con[mus_promoter_con$align != 0, ], aes(x = align)) +
    geom_histogram(binwidth = 1, fill = "#a893c1", colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Mus musculus promoter distribution") +
    xlab("align") +
    ylab("count") +
    theme(text = element_text(size = 40))
png("mus_promoter_dist_no0.png", width = 8000, height = 6000, res = 600)
plot(mus_promoter_dist_no0)
dev.off()
can_enhancer_count <- data.frame()
for (align in 1:24) {
    for (con in 1:align) {
        can_enhancer_count <- rbind(
            can_enhancer_count,
            c(align, con, nrow(can_enhancer_con[can_enhancer_con$align == align & can_enhancer_con$consrv == con, ]))
        )
    }
}
colnames(can_enhancer_count) <- c("align", "consrv", "count")
classIntervals(can_enhancer_count$count, n = 7, style = "jenks", warnLargeN = FALSE)
can_enhancer_heat <- ggplot(can_enhancer_count, aes(x = align, y = consrv, fill = count)) +
    geom_tile() +
    theme_bw() +
    scale_fill_steps2(low = "white", high = "#ad2c00", breaks = c(0, 193, 550, 1000, 1535, 2120, 2871, 3601)) +
    ggtitle("Canis lupus enhancer conservation") +
    xlab("align") +
    ylab("conservation") +
    theme(text = element_text(size = 35), legend.text = element_text(size = 15))
png("can_enhancer_heat.png", width = 8000, height = 6000, res = 600)
plot(can_enhancer_heat)
dev.off()
can_promoter_count <- data.frame()
for (align in 1:24) {
    for (con in 1:align) {
        can_promoter_count <- rbind(
            can_promoter_count,
            c(align, con, nrow(can_promoter_con[can_promoter_con$align == align & can_promoter_con$consrv == con, ]))
        )
    }
}
colnames(can_promoter_count) <- c("align", "consrv", "count")
classIntervals(can_promoter_count$count, n = 7, style = "jenks", warnLargeN = FALSE)
can_promoter_heat <- ggplot(can_promoter_count, aes(x = align, y = consrv, fill = count)) +
    geom_tile() +
    theme_bw() +
    scale_fill_steps2(low = "white", high = "#5a1993", breaks = c(140, 310, 441, 659, 1008, 1371, 2423)) +
    ggtitle("Canis lupus promoter conservation") +
    xlab("align") +
    ylab("conservation") +
    theme(text = element_text(size = 35), legend.text = element_text(size = 15))
png("can_promoter_heat.png", width = 8000, height = 6000, res = 600)
plot(can_promoter_heat)
dev.off()
mus_enhancer_count <- data.frame()
for (align in 1:24) {
    for (con in 1:align) {
        mus_enhancer_count <- rbind(
            mus_enhancer_count,
            c(align, con, nrow(mus_enhancer_con[mus_enhancer_con$align == align & mus_enhancer_con$consrv == con, ]))
        )
    }
}
colnames(mus_enhancer_count) <- c("align", "consrv", "count")
mus_enhancer_count_no1 <- mus_enhancer_count[-which(mus_enhancer_count$align == 1 | mus_enhancer_count$align == 2), ]
classIntervals(mus_enhancer_count_no1$count, n = 7, style = "jenks", warnLargeN = FALSE)
mus_enhancer_heat_no1 <- ggplot(mus_enhancer_count_no1, aes(x = align, y = consrv, fill = count)) +
    geom_tile() +
    theme_bw() +
    scale_fill_steps2(low = "white", high = "#ad2c00", breaks = c(129, 367, 692, 1098, 1602, 2241, 3302)) +
    ggtitle("Mus musculus enhancer conservation") +
    xlab("align") +
    ylab("conservation") +
    theme(text = element_text(size = 35), legend.text = element_text(size = 15))
png("mus_enhancer_heat_no1.png", width = 8000, height = 6000, res = 600)
plot(mus_enhancer_heat_no1)
dev.off()
mus_enhancer_heat <- ggplot(mus_enhancer_count, aes(x = align, y = consrv, fill = count)) +
    geom_tile() +
    theme_bw() +
    scale_fill_steps2(low = "white", high = "#ad2c00", breaks = c(0, 153, 444, 897, 1418, 2241, 3302, 10300)) +
    ggtitle("Mus musculus enhancer conservation") +
    xlab("align") +
    ylab("conservation") +
    theme(text = element_text(size = 35), legend.text = element_text(size = 15))
png("mus_enhancer_heat.png", width = 8000, height = 6000, res = 600)
plot(mus_enhancer_heat)
dev.off()
mus_promoter_count <- data.frame()
for (align in 1:24) {
    for (con in 1:align) {
        mus_promoter_count <- rbind(
            mus_promoter_count,
            c(align, con, nrow(mus_promoter_con[mus_promoter_con$align == align & mus_promoter_con$consrv == con, ]))
        )
    }
}
colnames(mus_promoter_count) <- c("align", "consrv", "count")
mus_promoter_heat <- ggplot(mus_promoter_count, aes(x = align, y = consrv, fill = log2(count + 0.1))) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient(low = "white", high = "#5a1993") +
    ggtitle("Mus musculus promoter conservation") +
    xlab("align") +
    ylab("conservation") +
    theme(text = element_text(size = 35), legend.text = element_text(size = 15))
png("mus_promoter_heat.png", width = 8000, height = 6000, res = 600)
plot(mus_promoter_heat)
dev.off()
mus_promoter_count_no1 <- mus_promoter_count[-which(mus_promoter_count$align == 1 | mus_promoter_count$align == 2), ]
classIntervals(mus_promoter_count_no1$count, n = 7, style = "jenks", warnLargeN = FALSE)
mus_promoter_heat_no1 <- ggplot(mus_promoter_count_no1, aes(x = align, y = consrv, fill = count)) +
    geom_tile() +
    theme_bw() +
    scale_fill_steps2(low = "white", high = "#5a1993", breaks = c(0, 58, 140, 236, 315, 411, 557, 822)) +
    ggtitle("Mus musculus promoter conservation") +
    xlab("align") +
    ylab("conservation") +
    theme(text = element_text(size = 35), legend.text = element_text(size = 15))
png("mus_promoter_heat_no1.png", width = 8000, height = 6000, res = 600)
plot(mus_promoter_heat_no1)
dev.off()

enhancer_count <- data.frame()
for (align in 1:24) {
    for (con in 1:align) {
        enhancer_count <- rbind(
            enhancer_count,
            c(align, con, nrow(enhancer_con[enhancer_con$align == align & enhancer_con$consrv == con, ]))
        )
    }
}
colnames(enhancer_count) <- c("align", "consrv", "count")
enhancer_breaks <- classIntervals(enhancer_count$count, n = 7, style = "jenks", warnLargeN = FALSE)
enhancer_heat <- ggplot(enhancer_count, aes(x = align, y = consrv, fill = count)) +
    geom_tile() +
    theme_bw() +
    scale_fill_steps(low = "white", high = "#ad2c00", breaks = as.vector(enhancer_breaks[2]$brks)) +
    ggtitle("enhancer conservation") +
    xlab("align") +
    ylab("conservation") +
    theme(text = element_text(size = 35), legend.text = element_text(size = 15))
png("enhancer_heat.png", width = 8000, height = 6000, res = 600)
plot(enhancer_heat)
dev.off()
promoter_count <- data.frame()
for (align in 1:24) {
    for (con in 1:align) {
        promoter_count <- rbind(
            promoter_count,
            c(align, con, nrow(promoter_con[promoter_con$align == align & promoter_con$consrv == con, ]))
        )
    }
}
colnames(promoter_count) <- c("align", "consrv", "count")

promoter_breaks <- classIntervals(promoter_count$count, n = 7, style = "jenks", warnLargeN = FALSE)
promoter_heat <- ggplot(promoter_count, aes(x = align, y = consrv, fill = count)) +
    geom_tile() +
    theme_bw() +
    scale_fill_steps(low = "white", high = "#5a1993", breaks = as.vector(promoter_breaks[2]$brks)) +
    ggtitle("promoter conservation") +
    xlab("align") +
    ylab("conservation") +
    theme(text = element_text(size = 35), legend.text = element_text(size = 15))
png("promoter_heat.png", width = 8000, height = 6000, res = 600)
plot(promoter_heat)
dev.off()

## plot the conservation of enhancer and promoter of each species
species <- c(
    "Canis_lupus",
    "Mustela_putorius",
    "Felis_catus",
    "Equus_asinus",
    "Equus_caballus",
    "Bos_taurus",
    "Ovis_aries",
    "Sus_scrofa",
    "Lama_glama",
    "Rhinolophus_pusillus",
    "Rhinolophus_ferrumequinum",
    "Hipposideros_larvatus",
    "Myotis_chinensis",
    "Myotis_ricketti",
    "Atelerix_albiventris",
    "Rattus_norvegicus",
    "Mus_musculus",
    "Cavia_porcellus",
    "Oryctolagus_cuniculus",
    "Macaca_mulatta",
    "Tupaia_belangeri",
    "Procavia_capensis",
    "Petaurus_breviceps",
    "Rhinopithecus_roxellana",
    "Neophocaena_asiaeorientalis"
)
for (sps in species) {
    sps_enhancer_con <- enhancer_con[grep(sps, rownames(enhancer_con)), ]
    sps_promoter_con <- promoter_con[grep(sps, rownames(promoter_con)), ]
    sps_enhancer_dist_no0 <- ggplot(sps_enhancer_con[sps_enhancer_con$align != 0, ], aes(x = align)) +
        geom_histogram(binwidth = 1, fill = "#ffb677", colour = "black") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(paste(sps, "enhancer distribution", sep = " ")) +
        xlab("align") +
        ylab("count") +
        theme(text = element_text(size = 40))
    png(paste0(sps, "_enhancer_dist_no0.png"), width = 8000, height = 6000, res = 600)
    plot(sps_enhancer_dist_no0)
    dev.off()
    sps_promoter_dist_no0 <- ggplot(sps_promoter_con[sps_promoter_con$align != 0, ], aes(x = align)) +
        geom_histogram(binwidth = 1, fill = "#a893c1", colour = "black") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(paste(sps, "promoter distribution", sep = " ")) +
        xlab("align") +
        ylab("count") +
        theme(text = element_text(size = 40))
    png(paste0(sps, "_promoter_dist_no0.png"), width = 8000, height = 6000, res = 600)
    plot(sps_promoter_dist_no0)
    dev.off()
    sps_enhancer_count <- data.frame()
    for (align in 1:24) {
        for (con in 1:align) {
            sps_enhancer_count <- rbind(
                sps_enhancer_count,
                c(align, con, nrow(sps_enhancer_con[sps_enhancer_con$align == align & sps_enhancer_con$consrv == con, ]))
            )
        }
    }
    colnames(sps_enhancer_count) <- c("align", "consrv", "count")
    enhancer_brk <- classIntervals(sps_enhancer_count$count, n = 7, style = "jenks", warnLargeN = FALSE)
    sps_enhancer_heat <- ggplot(sps_enhancer_count, aes(x = align, y = consrv, fill = count)) +
        geom_tile() +
        theme_bw() +
        scale_fill_steps2(low = "white", high = "#ad2c00", breaks = as.vector(enhancer_brk[2]$brks)) +
        ggtitle(paste(sps, "enhancer conservation", sep = " ")) +
        xlab("align") +
        ylab("conservation") +
        theme(text = element_text(size = 35), legend.text = element_text(size = 15))
    png(paste0(sps, "_enhancer_heat.png"), width = 8000, height = 6000, res = 600)
    plot(sps_enhancer_heat)
    dev.off()
    sps_promoter_count <- data.frame()
    for (align in 1:24) {
        for (con in 1:align) {
            sps_promoter_count <- rbind(
                sps_promoter_count,
                c(align, con, nrow(sps_promoter_con[sps_promoter_con$align == align & sps_promoter_con$consrv == con, ]))
            )
        }
    }
    colnames(sps_promoter_count) <- c("align", "consrv", "count")
    promoter_brk <- classIntervals(sps_promoter_count$count, n = 7, style = "jenks", warnLargeN = FALSE)
    sps_promoter_heat <- ggplot(sps_promoter_count, aes(x = align, y = consrv, fill = count)) +
        geom_tile() +
        theme_bw() +
        scale_fill_steps2(low = "white", high = "#5a1993", breaks = as.vector(promoter_brk[2]$brks)) +
        ggtitle(paste(sps, "promoter conservation", sep = " ")) +
        xlab("align") +
        ylab("conservation") +
        theme(text = element_text(size = 35), legend.text = element_text(size = 15))
    png(paste0(sps, "_promoter_heat.png"), width = 8000, height = 6000, res = 600)
    plot(sps_promoter_heat)
    dev.off()
}

# plot nj tree
enhancer_con2 <- enhancer_con
enhancer_con2$ID <- row.names(enhancer_con2)
enhancer_con2$sp <- paste(str_split_fixed(enhancer_con2$ID, "_", 8)[, 1], str_split_fixed(enhancer_con2$ID, "_", 8)[, 2], sep = "_")
enhancer_con2$tis <- str_split_fixed(enhancer_con2$ID, "_", 8)[, 3]
enhancer_con2$ele <- str_split_fixed(enhancer_con2$ID, "_", 8)[, 4]

nj_species <- c(
    "Macaca_mulatta", "Rhinopithecus_roxellana",
    "Mus_musculus", "Rattus_norvegicus",
    "Oryctolagus_cuniculus", "Canis_lupus",
    "Felis_catus", "Bos_taurus",
    "Sus_scrofa", "Petaurus_breviceps"
)

# enhancer
enhancer_con_matrix <- as.data.frame(matrix(nrow = length(species), ncol = length(species)))
colnames(enhancer_con_matrix) <- species
rownames(enhancer_con_matrix) <- species
for (spe1 in species) {
    spe1_peak_num <- sum(enhancer_con2$sp == spe1)
    for (spe2 in species) {
        enhancer_con_matrix[spe1, spe2] <-
            sum(enhancer_con2[
                which(enhancer_con2$sp == spe1 & enhancer_con2[, paste0(spe2, "_con")] == 1),
                paste0(spe2, "_con")
            ]) / spe1_peak_num
    }
}
enhancer_con_matrix[is.na(enhancer_con_matrix)] <- 0
enhancer_con_matrix2 <- as.data.frame(matrix(nrow = length(species), ncol = length(species)))
colnames(enhancer_con_matrix2) <- species
rownames(enhancer_con_matrix2) <- species
for (spe1 in species) {
    for (spe2 in species) {
        enhancer_con_matrix2[spe1, spe2] <- mean(c(enhancer_con_matrix[spe1, spe2], enhancer_con_matrix[spe2, spe1]))
        enhancer_con_matrix2[spe2, spe1] <- enhancer_con_matrix2[spe1, spe2]
    }
}
enhancer_con_matrix2[is.na(enhancer_con_matrix2)] <- 0
enhancer_con_matrix_nj <- 1 - enhancer_con_matrix2
enhancer_con_nj <- nj(as.matrix(enhancer_con_matrix_nj))
enhancer_con_nj$tip.label <- species
enhancer_nj2 <- root(enhancer_con_nj, outgroup = "Petaurus_breviceps")
plot(enhancer_nj2, type = "phylogram", main = "enhancer nj tree")
enhancer_nj_p <- ggtree(enhancer_nj2, layout = "rectangular", size = 0.8) +
    geom_tiplab(size = 10) +
    theme_tree2() + xlim(NA, 1) + ggtitle("enhancer nj tree") +
    theme(plot.title = element_text(hjust = 0.5, size = 35))
enhancer_con_matrix_nj_10sp <- enhancer_con_matrix_nj[nj_species, nj_species]
enhancer_nj_10sp <- nj(as.matrix(enhancer_con_matrix_nj_10sp))
enhancer_nj_10sp$tip.label <- nj_species
enhancer_nj_10sp2 <- root(enhancer_nj_10sp, outgroup = "Petaurus_breviceps")
enhancer_nj_10sp_p <- ggtree(enhancer_nj_10sp2, layout = "rectangular", size = 0.8) +
    geom_tiplab(size = 10) +
    theme_tree2() + xlim(NA, 1) + ggtitle("10sp enhancer nj tree") +
    theme(plot.title = element_text(hjust = 0.5, size = 35))
png("enhancer_nj_10sp.png", height = 6000, width = 4000, res = 600)
plot(enhancer_nj_10sp_p)
dev.off()

# promoter
promoter_con2 <- promoter_con
promoter_con2$ID <- row.names(promoter_con2)
promoter_con2$sp <- paste(str_split_fixed(promoter_con2$ID, "_", 8)[, 1], str_split_fixed(promoter_con2$ID, "_", 8)[, 2], sep = "_")
promoter_con2$tis <- str_split_fixed(promoter_con2$ID, "_", 8)[, 3]
promoter_con2$ele <- str_split_fixed(promoter_con2$ID, "_", 8)[, 4]
promoter_con_matrix <- as.data.frame(matrix(nrow = length(species), ncol = length(species)))
colnames(promoter_con_matrix) <- species
rownames(promoter_con_matrix) <- species
for (spe1 in species) {
    spe1_peak_num <- sum(promoter_con2$sp == spe1)
    for (spe2 in species) {
        promoter_con_matrix[spe1, spe2] <- sum(promoter_con2[which(promoter_con2$sp == spe1 & promoter_con2[, paste0(spe2, "_con")] == 1), paste0(spe2, "_con")]) / spe1_peak_num
    }
}
promoter_con_matrix[is.na(promoter_con_matrix)] <- 0
promoter_con_matrix2 <- as.data.frame(matrix(nrow = length(species), ncol = length(species)))
colnames(promoter_con_matrix2) <- species
rownames(promoter_con_matrix2) <- species
for (spe1 in species) {
    for (spe2 in species) {
        promoter_con_matrix2[spe1, spe2] <- mean(c(promoter_con_matrix[spe1, spe2], promoter_con_matrix[spe2, spe1]))
        promoter_con_matrix2[spe2, spe1] <- promoter_con_matrix2[spe1, spe2]
    }
}
promoter_con_matrix2[is.na(promoter_con_matrix2)] <- 0
promoter_con_matrix_nj <- 1 - promoter_con_matrix2
promoter_con_nj <- nj(as.matrix(promoter_con_matrix_nj))
promoter_con_nj$tip.label <- species
plot(promoter_con_nj, type = "phylogram", main = "promoter nj tree")
promoter_nj2 <- root(promoter_con_nj, outgroup = "Petaurus_breviceps")
promoter_nj_p <- ggtree(promoter_nj2, layout = "rectangular", size = 0.8) +
    geom_tiplab(size = 10) +
    theme_tree2() + xlim(NA, 1) + ggtitle("promoter nj tree") +
    theme(plot.title = element_text(hjust = 0.5, size = 35))
promoter_con_matrix_nj_10sp <- promoter_con_matrix_nj[nj_species, nj_species]
promoter_nj_10sp <- nj(as.matrix(promoter_con_matrix_nj_10sp))
promoter_nj_10sp$tip.label <- nj_species
promoter_nj_10sp2 <- root(promoter_nj_10sp, outgroup = "Petaurus_breviceps")
promoter_nj_10sp_p <- ggtree(promoter_nj_10sp2, layout = "rectangular", size = 0.8) +
    geom_tiplab(size = 10) +
    theme_tree2() + xlim(NA, 1) + ggtitle("10sp promoter nj tree") +
    theme(plot.title = element_text(hjust = 0.5, size = 35))
png("promoter_nj_10sp.png", height = 6000, width = 4000, res = 600)
plot(promoter_nj_10sp_p)
dev.off()
