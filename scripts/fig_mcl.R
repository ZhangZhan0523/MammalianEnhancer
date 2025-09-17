## draw pictures show basic and importance information about blast, mcl
## re-mcl clustering, alignment

library(ggplot2)


# first of all load the data, let us see how many elements are involved
mcl_res <- "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl.mcl"
clusters <- readLines(mcl_res)
stats <- lapply(seq_along(clusters), function(i) {
    elements <- unlist(strsplit(clusters[i], "\t"))
    enh_count <- sum(grepl("enh", elements)) # 假设enhancer标识包含"enh"
    pro_count <- sum(grepl("pro", elements)) # 假设promoter标识包含"pro"
    data.frame(
        cluster_id = i,
        cluster_size = length(elements),
        enhancers = enh_count,
        promoters = pro_count
    )
})

# 转换为数据框
cluster_stats <- do.call(rbind, stats)

# 查看前10个聚类的统计结果
head(cluster_stats, 10)

# save cluster_stats
save(cluster_stats, file = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl_stats.RData")
load("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl_stats.RData")
