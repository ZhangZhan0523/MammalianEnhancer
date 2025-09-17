# we summarize the overlap information of the enhanceratlas data,
# plot a barplot to show the fraction of overlapped enhancers in each tissue

library(ggplot2)
library(magrittr)

overlap_data <- read.csv("/media/Data/zhangz/chip/analysis/summary2/open_source/otherArticle/overlap_info.csv", header = TRUE)
overlap_data$research <- factor(overlap_data$research, levels = c("this_study", "ref1", "ref2", "ENCODE_ref"))
overlap_data$overlap_proportion <- gsub("%", "", overlap_data$overlap_proportion) %>% as.numeric() / 100
overlap_bar <- ggplot(overlap_data, aes(x = N, y = overlap_proportion, fill = research, color = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 0.5)) +
    labs(
        title = "Fraction of mouse enhancers from EnhancerAtlas overlapped with other datasets",
        x = "Tissue",
        y = "Fraction of overlapped enhancers"
    ) +
    scale_fill_manual(values = c("this_study" = "#ef4343", "ref1" = "grey70", "ref2" = "grey60", "ENCODE_ref" = "grey50"), labels = c("this_study" = "This study", "ref1" = "ref1", "ref2" = "ref2", "ENCODE_ref" = "ENCODE reference epigenomes")) +
    theme(legend.position = "top") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(labels = c("Liver", "Kidney", "Brain"), breaks = c(3, 7, 10.5)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(fill = guide_legend(title = "Dataset", nrow = 2), color = guide_legend(title = "Method", nrow = 2))
ggsave("/media/Data/zhangz/chip/analysis/summary2/open_source/otherArticle/overlap_barplot.pdf", overlap_bar, width = 6, height = 4)

overlap_bar2 <- ggplot(overlap_data, aes(x = N, y = overlap_proportion, fill = research, linetype = developmental_stage)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 0.5)) +
    labs(
        title = "Fraction of mouse enhancers from EnhancerAtlas overlapped with other datasets",
        x = "Tissue",
        y = "Fraction of overlapped enhancers"
    ) +
    scale_fill_manual(values = c("this_study" = "#ef4343", "ref1" = "grey90", "ref2" = "grey80", "ENCODE_ref" = "grey70"), labels = c("this_study" = "This study", "ref1" = "ref1", "ref2" = "ref2", "ENCODE_ref" = "ENCODE reference epigenomes")) +
    scale_linetype_manual(values = c("adult" = "solid", "neonate" = "longdash"), labels = c("adult" = "Adult", "neonate" = "Neonate")) +
    theme(legend.position = "top") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(labels = c("Liver", "Kidney", "Brain"), breaks = c(3, 7, 10.5)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.line = element_line(color = "grey70", linewidth = 0.25)) +
    theme(axis.ticks = element_line(color = "grey70", linewidth = 0.25)) +
    guides(fill = guide_legend(title = "Dataset", nrow = 2), color = guide_legend(title = "Developmental stage", nrow = 2))
ggsave("/media/Data/zhangz/chip/analysis/summary2/open_source/otherArticle/overlap_barplot2.pdf", overlap_bar2, width = 6, height = 4)
