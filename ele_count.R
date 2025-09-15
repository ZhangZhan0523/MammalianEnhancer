library(ggplot2)
library(reshape2) # 数据转换
require(scales) # 数据缩放
library(ggtree) # 聚类
library(aplot) # 拼图
library(dplyr)
library(patchwork)
library(cowplot)

info <- read.csv("/media/Data/zhangz/chip/scripts2/info/info_using_ele2.csv", sep = ",", header = T)
data <- melt(info, id.vars = c("species", "tissue"), measure.vars = c("enhancer", "promoter", "H3K4me3_only"))
gen_size <- data.frame(unique(data$species), 0)
colnames(gen_size) <- c("species", "size")
rownames(gen_size) <- gen_size$species
for (i in seq_len(nrow(gen_size))) {
    if (gen_size[i, 1] == "Macaca_mulatta") {
        fname <- "/data1/Genome/Macaca_fascicularis/chrom.sizes"
    } else {
        fname <- paste0("/data1/Genome/", gen_size[i, 1], "/chrom.sizes")
    }
    tmp <- read.csv(fname, sep = "\t", header = F)
    gen_size[i, 2] <- sum(tmp$V2)
}
data$id <- paste(data$species, data$tissue, sep = "_")
data$value_norm <- data$value * 1000000
for (i in seq_len(nrow(data))) {
    data[i, "value_norm"] <- data[i, "value_norm"] / gen_size[data[i, "species"], 2]
}
tree <- read.tree("/media/Data/zhangz/chip/scripts2/info/sps.nwk")
sp_info <- read.csv("/media/Data/zhangz/chip/scripts/info/species.csv", header = T, row.names = 1)
sp_info <- sp_info[, c("order", "species")]
colnames(sp_info) <- c("order", "label")
grouptree <- full_join(tree, sp_info, by = "label")
as_tibble(grouptree) %>%
    as.data.frame()
col <- paletteer::paletteer_d("ggthemes::Classic_20", n = 11)
# tree_p <- ggtree(grouptree, branch.length = "none", aes(color = order)) +
#     theme(legend.position = "none") +
#     scale_y_reverse() + coord_flip() + scale_fill_manual(values = rev(col)) #+geom_tiplab()
tree_p_2 <- ggtree(tree, branch.length = "none") +
    theme(legend.position = "none") +
    scale_y_reverse() + coord_flip() #+geom_tiplab()
# stack_bar_p <- ggplot(data = data, aes(id, value, fill = variable)) +
#     geom_bar(aes(species), position = "stack", stat = "identity") +
#     scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
#     theme_bw() +
#     labs(x = "", y = "active regulatory regions") +
#     scale_y_continuous(expand = c(0, 0)) +
#     theme(
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
#         axis.title = element_text(size = 18)
#     ) +
#     facet_wrap(~tissue, nrow = 1, scales = "fixed")
# stack_liver_p <- ggplot(data = data[which(data$tissue == "Liver"), ], aes(id, value, fill = variable)) +
#     geom_bar(aes(species), position = "stack", stat = "identity") +
#     scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
#     theme_tree2() +
#     labs(x = "", y = "active regulatory regions") +
#     scale_y_continuous(expand = c(0, 0)) +
#     theme(
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
#         axis.title = element_text(size = 18)
#     ) +
#     facet_wrap(~tissue, nrow = 1, scales = "fixed")
# stack_kidney_p <- ggplot(data = data[which(data$tissue == "Kidney"), ], aes(id, value, fill = variable)) +
#     geom_bar(aes(species), position = "stack", stat = "identity") +
#     scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
#     theme_tree2() +
#     labs(x = "", y = "active regulatory regions") +
#     scale_y_continuous(expand = c(0, 0)) +
#     theme(
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
#         axis.title = element_text(size = 18)
#     ) +
#     facet_wrap(~tissue, nrow = 1, scales = "fixed")
# stack_brain_p <- ggplot(data = data[which(data$tissue == "Brain"), ], aes(id, value, fill = variable)) +
#     geom_bar(aes(species), position = "stack", stat = "identity") +
#     scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
#     theme_tree2() +
#     labs(x = "", y = "active regulatory regions") +
#     scale_y_continuous(expand = c(0, 0)) +
#     theme(
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
#         axis.title = element_text(size = 18)
#     ) +
#     facet_wrap(~tissue, nrow = 1, scales = "fixed")
# stack_bar_p2 <- stack_bar_p + xlim2(tree_p)
# stack_brain_p2 <- stack_brain_p + xlim2(tree_p_2)
# stack_kidney_p2 <- stack_kidney_p + xlim2(tree_p_2)
# stack_liver_p2 <- stack_liver_p + xlim2(tree_p_2)
brain <- data[which(data$tissue == "Brain"), ]
brain <- dplyr::arrange(brain, species)
brain$species <- factor(brain$species, levels = c(
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
))
stack_brain_p3 <- ggplot(data = brain, aes(id, value, fill = variable)) +
    geom_bar(aes(species), position = "stack", stat = "identity") +
    scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
    theme_bw() +
    labs(x = "", y = "active regulatory regions") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
        axis.title = element_text(size = 18)
    ) +
    facet_wrap(~tissue, nrow = 1, scales = "fixed")
stack_brain_n_p <- ggplot(data = brain, aes(id, value_norm, fill = variable)) +
    geom_bar(aes(species), position = "stack", stat = "identity") +
    scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
    theme_bw() +
    labs(x = "", y = "active regulatory regions per million bp") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
        axis.title = element_text(size = 18)
    ) +
    facet_wrap(~tissue, nrow = 1, scales = "fixed")
kidney <- data[which(data$tissue == "Kidney"), ]
kidney <- dplyr::arrange(kidney, species)
kidney$species <- factor(kidney$species, levels = c(
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
))

stack_kidney_p3 <- ggplot(data = kidney, aes(id, value, fill = variable)) +
    geom_bar(aes(species), position = "stack", stat = "identity") +
    scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
    theme_bw() +
    labs(x = "", y = "active regulatory regions") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
        axis.title = element_text(size = 18)
    ) +
    facet_wrap(~tissue, nrow = 1, scales = "fixed")
stack_kidney_n_p <- ggplot(data = kidney, aes(id, value_norm, fill = variable)) +
    geom_bar(aes(species), position = "stack", stat = "identity") +
    scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
    theme_bw() +
    labs(x = "", y = "active regulatory regions per million bp") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
        axis.title = element_text(size = 18)
    ) +
    facet_wrap(~tissue, nrow = 1, scales = "fixed")
liver <- data[which(data$tissue == "Liver"), ]
liver <- dplyr::arrange(liver, species)
liver$species <- factor(liver$species, levels = c(
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
))
stack_liver_p3 <- ggplot(data = liver, aes(id, value, fill = variable)) +
    geom_bar(aes(species), position = "stack", stat = "identity") +
    scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
    theme_bw() +
    labs(x = "", y = "active regulatory regions") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
        axis.title = element_text(size = 18)
    ) +
    facet_wrap(~tissue, nrow = 1, scales = "fixed")
stack_liver_n_p <- ggplot(data = liver, aes(id, value_norm, fill = variable)) +
    geom_bar(aes(species), position = "stack", stat = "identity") +
    scale_fill_manual(values = c("#ff9740", "#9c76b0", "#00a6d3")) +
    theme_bw() +
    labs(x = "", y = "active regulatory regions per million bp") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
        axis.title = element_text(size = 18)
    ) +
    facet_wrap(~tissue, nrow = 1, scales = "fixed")
concat_p <- plot_grid(stack_liver_p3, stack_kidney_p3, stack_brain_p3,
    tree_p_2, tree_p_2, tree_p_2,
    ncol = 3,
    rel_heights = c(1, 0.3), align = "v", axis = "lr"
)
concat_n_p <- plot_grid(stack_liver_n_p, stack_kidney_n_p, stack_brain_n_p,
    tree_p_2, tree_p_2, tree_p_2,
    ncol = 3,
    rel_heights = c(1, 0.3), align = "v", axis = "lr"
)
png("/media/Data/zhangz/chip/scripts2/info/ele_count.png", height = 6000, width = 12000, res = 600)
plot(concat_p)
dev.off()
png("/media/Data/zhangz/chip/scripts2/info/ele_count_norm.png", height = 6000, width = 12000, res = 600)
plot(concat_n_p)
dev.off()
