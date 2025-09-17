library(ggplot2)
library(reshape2, help, pos = 2, lib.loc = NULL)
library(car)

setwd("/media/usb1/chip/analysis/summary2/")
info_df <- read.csv("/media/Data/zhangz/chip/scripts2/info/info_using_ele2.csv", header = T)
info_df$H3K4me3 <- 0
info_df$H3K27ac <- 0
for (i in 1:nrow(info_df)) {
    me3 <- read.csv(paste0("/media/usb1/chip/analysis/", info_df$species[i], "/anno/peaks/", info_df$species[i], "_", info_df$tissue[i], "_H3K4me3_overlap.optimal_peak.narrowPeak"), header = F, sep = "\t")
    info_df$H3K4me3[i] <- nrow(me3)
    ac27 <- read.csv(paste0("/media/usb1/chip/analysis/", info_df$species[i], "/anno/peaks/", info_df$species[i], "_", info_df$tissue[i], "_H3K27ac_overlap.optimal_peak.narrowPeak"), header = F, sep = "\t")
    info_df$H3K27ac[i] <- nrow(ac27)
}
write.csv(info_df, file = "/media/Data/zhangz/chip/scripts2/info/info_using_ele2.csv", row.names = FALSE, quote = FALSE)

ele_stats <- info_df[, c(1, 2, 34, 35, 36, 37, 38)]
write.csv(ele_stats, file = "/media/usb1/chip/analysis/summary2/ele_stats.csv", row.names = FALSE, quote = FALSE)
ele_stats$mean_enhancer_length <- 0
ele_stats$mid_enhancer_length <- 0
ele_stats$mean_promoter_length <- 0
ele_stats$mid_promoter_length <- 0
ele_stats$mean_enhancer_align <- 0
ele_stats$mid_enhancer_align <- 0
ele_stats$mean_promoter_align <- 0
ele_stats$mid_promoter_align <- 0
ele_stats$mean_enhancer_consrv <- 0
ele_stats$mid_enhancer_consrv <- 0
ele_stats$mean_promoter_consrv <- 0
ele_stats$mid_promoter_consrv <- 0
ele_stats$mean_promoter_fc <- 0
ele_stats$mid_promoter_fc <- 0


for (i in 1:nrow(ele_stats)) {
    enhancer <- read.csv(paste0(
        "/media/Data/zhangz/chip/analysis/", ele_stats$species[i], "/compare2/",
        ele_stats$species[i], "_", ele_stats$tissue[i], "_enhancer_info.csv"
    ), header = TRUE, sep = ",")
    ele_stats$mean_enhancer_length[i] <- mean(enhancer$length)
    ele_stats$mid_enhancer_length[i] <- median(enhancer$length)
    ele_stats$mean_enhancer_align[i] <- mean(enhancer$align)
    ele_stats$mid_enhancer_align[i] <- median(enhancer$align)
    ele_stats$mean_enhancer_consrv[i] <- mean(enhancer$consrv)
    ele_stats$mid_enhancer_consrv[i] <- median(enhancer$consrv)
    ele_stats$mean_enhancer_fc[i] <- mean(enhancer$mean_fc)
    ele_stats$mid_enhancer_fc[i] <- median(enhancer$mean_fc)
    ele_stats$mean_enhancer_gc[i] <- mean(enhancer$gc)
    ele_stats$mid_enhancer_gc[i] <- median(enhancer$gc)
    ele_stats$mean_enhancer_tss[i] <- mean(enhancer$distanceToTSS)
    ele_stats$mid_enhancer_tss[i] <- median(enhancer$distanceToTSS)
    ele_stats$mean_enhancer_motif_num[i] <- mean(enhancer$motif_num)
    ele_stats$mid_enhancer_motif_num[i] <- median(enhancer$motif_num)
    promoter <- read.csv(paste0(
        "/media/Data/zhangz/chip/analysis/", ele_stats$species[i], "/compare2/",
        ele_stats$species[i], "_", ele_stats$tissue[i], "_promoter_info.csv"
    ), header = TRUE, sep = ",")
    ele_stats$mean_promoter_length[i] <- mean(promoter$length)
    ele_stats$mid_promoter_length[i] <- median(promoter$length)
    ele_stats$mean_promoter_align[i] <- mean(promoter$align)
    ele_stats$mid_promoter_align[i] <- median(promoter$align)
    ele_stats$mean_promoter_consrv[i] <- mean(promoter$consrv)
    ele_stats$mid_promoter_consrv[i] <- median(promoter$consrv)
    ele_stats$mean_promoter_fc[i] <- mean(promoter$mean_fc)
    ele_stats$mid_promoter_fc[i] <- median(promoter$mean_fc)
    ele_stats$mean_promoter_gc[i] <- mean(promoter$gc)
    ele_stats$mid_promoter_gc[i] <- median(promoter$gc)
    ele_stats$mean_promoter_tss[i] <- mean(promoter$distanceToTSS)
    ele_stats$mid_promoter_tss[i] <- median(promoter$distanceToTSS)
    ele_stats$mean_promoter_motif_num[i] <- mean(promoter$motif_num)
    ele_stats$mid_promoter_motif_num[i] <- median(promoter$motif_num)
}
for (i in c(
    "mean_enhancer_length", "mid_enhancer_length", "mean_promoter_length", "mid_promoter_length",
    "mean_enhancer_align", "mid_enhancer_align", "mean_promoter_align", "mid_promoter_align", "mean_enhancer_consrv",
    "mid_enhancer_consrv", "mean_promoter_consrv", "mid_promoter_consrv"
)) {
    dist_p <- ggplot(ele_stats, aes_string(x = i, fill = "tissue")) +
        geom_density(alpha = 0.5) +
        theme_bw() +
        labs(x = i, y = "Density", title = paste("Density plot of", i, "in different tissues", sep = " ")) +
        theme(plot.title = element_text(hjust = 0.5))
    pdf(paste(i, "density.pdf", sep = "_"), width = 6, height = 4)
    print(dist_p)
    dev.off()
}
for (i in c(
    "mean_enhancer_fc", "mid_enhancer_fc", "mean_promoter_fc", "mid_promoter_fc",
    "mean_enhancer_gc", "mid_enhancer_gc", "mean_promoter_gc", "mid_promoter_gc",
    "mean_enhancer_tss", "mid_enhancer_tss", "mean_promoter_tss", "mid_promoter_tss",
    "mean_enhancer_motif_num", "mid_enhancer_motif_num", "mean_promoter_motif_num", "mid_promoter_motif_num"
)) {
    dist_p <- ggplot(ele_stats, aes_string(x = i, fill = "tissue")) +
        geom_density(alpha = 0.5) +
        theme_bw() +
        labs(x = i, y = "Density", title = paste("Density plot of", i, "in different tissues", sep = " ")) +
        theme(plot.title = element_text(hjust = 0.5))
    pdf(paste(i, "density.pdf", sep = "_"), width = 6, height = 4)
    print(dist_p)
    dev.off()
}
write.csv(ele_stats, file = "/media/usb1/chip/analysis/summary2/ele_stats.csv", row.names = FALSE, quote = FALSE)

ele_stats_brain <- ele_stats[ele_stats$tissue == "Brain", ]
ele_stats_liver <- ele_stats[ele_stats$tissue == "Liver", ]
ele_stats_kidney <- ele_stats[ele_stats$tissue == "Kidney", ]
mean(ele_stats_brain$mean_enhancer_length)
mean(ele_stats_liver$mean_enhancer_length)
mean(ele_stats_kidney$mean_enhancer_length)
leveneTest(mean_enhancer_length ~ tissue, data = ele_stats, center = mean)
leveneTest(mean_enhancer_align ~ tissue, data = ele_stats, center = mean)
leveneTest(mean_enhancer_consrv ~ tissue, data = ele_stats, center = mean)
leveneTest(mean_promoter_length ~ tissue, data = ele_stats, center = mean)
leveneTest(mean_promoter_align ~ tissue, data = ele_stats, center = mean)
leveneTest(mean_promoter_consrv ~ tissue, data = ele_stats, center = mean)
t.test(mean_enhancer_length ~ tissue, data = ele_stats[which(ele_stats$tissue != "Liver"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_enhancer_align ~ tissue, data = ele_stats[which(ele_stats$tissue != "Liver"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_enhancer_consrv ~ tissue, data = ele_stats[which(ele_stats$tissue != "Liver"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_length ~ tissue, data = ele_stats[which(ele_stats$tissue != "Liver"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_align ~ tissue, data = ele_stats[which(ele_stats$tissue != "Liver"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_consrv ~ tissue, data = ele_stats[which(ele_stats$tissue != "Liver"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_enhancer_length ~ tissue, data = ele_stats[which(ele_stats$tissue != "Kidney"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_enhancer_align ~ tissue, data = ele_stats[which(ele_stats$tissue != "Kidney"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_enhancer_consrv ~ tissue, data = ele_stats[which(ele_stats$tissue != "Kidney"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_length ~ tissue, data = ele_stats[which(ele_stats$tissue != "Kidney"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_align ~ tissue, data = ele_stats[which(ele_stats$tissue != "Kidney"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_consrv ~ tissue, data = ele_stats[which(ele_stats$tissue != "Kidney"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_enhancer_length ~ tissue, data = ele_stats[which(ele_stats$tissue != "Brain"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_enhancer_align ~ tissue, data = ele_stats[which(ele_stats$tissue != "Brain"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_enhancer_consrv ~ tissue, data = ele_stats[which(ele_stats$tissue != "Brain"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_length ~ tissue, data = ele_stats[which(ele_stats$tissue != "Brain"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_align ~ tissue, data = ele_stats[which(ele_stats$tissue != "Brain"), ], paired = TRUE, var.equal = TRUE)
t.test(mean_promoter_consrv ~ tissue, data = ele_stats[which(ele_stats$tissue != "Brain"), ], paired = TRUE, var.equal = TRUE)
