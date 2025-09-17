#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
总体的目标是把enhancer序列转换成motif+间隔长度的形式，然后进行基于motif的enhancer比对，
从而发现一组同源enhancer的序列pattern。步骤分解为：
1. 将xstreme调用的fimo发现的denovo motif的occurrence和单独进行的known motif的fimo结果中的tsv文件进行合并，
    获得总体的enhancer中的motif的位置信息。
2. 按照cre的名字进行排序分组，对于每一个cre，判断是否出现motif的重叠，如果出现了motif的重叠，依次根据motif的q-value，
    p-value，score和长度进行排序，保留q-value最小，p-value最小，score最高，长度最大的motif。
3. 获取每个cre序列的长度，将每个cre的位置信息转换为motif+间隔长度的形式，然后进行基于motif的enhancer比对
4. 确定cre比对的打分矩阵，设计动态规划算法

这个是整理后的代码，封装了函数，主要执行对所有哺乳动物的最大cre cluster进行比对。

"""


import os
import time
import pandas as pd
from loguru import logger
import subprocess
import sys
from Bio import SeqIO
import pandas as pd
import concurrent.futures
from collections import defaultdict
import random

# from collections import OrderedDict
import pickle

# from bitarray import bitarray

# import itertools
import itertools

# 在函数顶部添加装饰器
from numba import jit

# 优化点1：使用numpy加速矩阵运算
import numpy as np
from heapq import heappush, heappop
from itertools import combinations
from functools import lru_cache
import psutil
import argparse
import hashlib

# my functions
# from enhancer_align import (
#     transform_sequence,
#     get_fasta_sequences,
#     progressive_motif_alignment,
#     motif_global_alignment,
# )


logger.add(
    sys.stdout,
    format="{time} {level} {message}",
    filter="my_module",
    level="INFO",
    diagnose=True,
    backtrace=True,
)

meme_cluster = pd.read_csv(
    "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/JASPAR_denovo_motifs_clusters.tab",
    sep="\t",
    header=0,
)
fimo_res = pd.read_csv(
    "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/fimo_filtered_res.tsv",
    sep="\t",
    header=0,
)

parser = argparse.ArgumentParser(description="Enhancer Alignment Parameters")
parser.add_argument(
    "--load-existing", action="store_true", help="Load precomputed matrix and tree"
)
parser.add_argument(
    "--save-dir",
    type=str,
    default="./cache",
    help="Directory for saving computed data",
)
parser.add_argument(
    "--outdir",
    type=str,
    help="Output directory for the results",
    default="/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/top",
)
parser.add_argument(
    "--group",
    type=str,
    help="Group name for the results",
    default="/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/top/group1/group1.txt",
)
args = parser.parse_args()

sys.setrecursionlimit(1000000)  # 设置递归深度限制


@logger.catch
def dict_to_dataframe(data_dict):
    """将字典转换为DataFrame"""
    df = pd.DataFrame.from_dict(data_dict, orient="index")
    df.reset_index(inplace=True)
    df.columns = ["Sequence Name"] + [
        f"spacer_{i//2}" if i % 2 == 1 else f"motif_{i//2}"
        for i in range(1, len(df.columns))
    ]
    return df


@logger.catch
def calculate_cluster_scores(df):
    """计算cluster列的统计指标，返回DataFrame"""
    records = []
    for col in df.columns:
        if col.startswith("motif_"):
            raw_series = df[col]
            valid_series = raw_series[raw_series != "-"]
            # 收集指标到列表
            records.append(
                {
                    "ColumnName": col,
                    "TotalRows": len(raw_series),
                    "ValidCount": len(valid_series),
                    "UniqueCount": valid_series.nunique(),
                    "TopFrequency": (
                        valid_series.value_counts(normalize=True).max()
                        if not valid_series.empty
                        else 0
                    ),
                }
            )
    # 转换为DataFrame
    return pd.DataFrame.from_records(records).set_index("ColumnName")


@logger.catch
def calculate_cohesion2(alignment):
    keys = list(alignment["cluster_alignments"].keys())
    align = (
        alignment["cluster_alignments"][keys[-1]]
        if len(keys) > 1  # Check number of alignment groups
        else alignment["cluster_alignments"][keys[0]]
    )
    align = dict_to_dataframe(align)
    # 计算比对的簇一致性分数
    score_df = calculate_cluster_scores(align)
    # 计算总分
    # TopFrequency * ValidCount
    score = (score_df["TopFrequency"] * score_df["ValidCount"]).sum()
    return score


@logger.catch
def get_fasta_sequences(fasta_file):
    """
    获取fasta文件中的序列长度
    使用dataframe存储结果
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_name = record.id
        seq_length = len(record.seq)
        sequences.append((seq_name, seq_length))

    df = pd.DataFrame(sequences, columns=["Sequence Name", "Sequence Length"])
    return df


# Add at top level (not inside any function)
@logger.catch
def kmer_similarity(a, b):
    """Top-level function for pickle compatibility"""
    return len(a & b) / max(len(a | b), 1)


# @logger.catch
# def calculate_distance(a, b, c, d):
#     return 0.3 * (1 - kmer_similarity(a, b)) + 0.7 * (1 - kmer_similarity(c, d))
@logger.catch
def calculate_distance(a, b, c, d):
    # 位运算替代集合操作
    a = int(a, 2) if isinstance(a, str) else a
    b = int(b, 2) if isinstance(b, str) else b
    c = int(c, 2) if isinstance(c, str) else c
    d = int(d, 2) if isinstance(d, str) else d
    m_intersection = bin(a & b).count("1")
    m_union = bin(a | b).count("1")
    c_intersection = bin(c & d).count("1")
    c_union = bin(c | d).count("1")
    distance = 0.3 * (1 - m_intersection / max(m_union, 1)) + 0.7 * (
        1 - c_intersection / max(c_union, 1)
    )
    return distance


def mem_monitor():
    return f"Memory used: {psutil.Process().memory_info().rss//1024//1024//1024}GB"


@logger.catch
def motif_global_alignment(
    seq_a,
    seq_b,
    cluster_a,
    cluster_b,
    cluster_tab=meme_cluster,
    match_score=10,
    cluster_match_score=8,
    mismatch_penalty=0.1,
    # mismatch和gap_penalty相同的惩罚分数，允许错配，之前mismatch是8
    # 在合并比对时会出现两条序列的motif无法比对上，于是各自加gap，合并后出现连续两个motif的情况
    # 如果mismatch罚分太高，就不会出现mismatch的情况
    same_spacer_score=2,
    sim_spacer_score=1,
    spacer_coeff=0.001,
    gap_penalty=5,
):
    """
    基于动态规划的motif序列比对算法
    输入格式示例：
    seq_a = {'enh_b_45423_B_promoter_867': [34, '22-GCGCCACCAGGGA', 0, '8-AGTGGTTAAGAATCC', 5, '7-CCAATGCAGGGGACA', 2, '2-ACCAGGGCTCGAACC', 5, '6-AAGATCCCACATGCC', 18]}
    seq_b = {'enh_b_31_B_Peak_10011_rev': [0, '22-GCGCCACCAGGGA', 0, '8-AGTGGTTAAGAATCC', 5, '7-CCAATGCAGGGGACA', 2, '2-ACCAGGGCTCGAACC', 5, '6-AAGATCCCACATGCC', 0, '5-GGAGCAACTAAGCCC', 3, '1-TCAGTAGTTGTGGCA', 25, 'GCTTCYCTKGTKG', 21, 'CTTCTCATTGCGGTG', 5, 'GTTGCKGYGMGCGGG', 1]}
    cluster_a = {'enh_b_45423_B_promoter_867': [34, 'cluster_015', 0, 'denovo_5', 5, 'denovo_3', 2, 'denovo_0', 5, 'cluster_004', 18]}
    cluster_b = {'enh_b_31_B_Peak_10011_rev': [0, 'cluster_015', 0, 'denovo_5', 5, 'denovo_3', 2, 'denovo_0', 5, 'cluster_004', 0, 'denovo_1', 3, 'cluster_109', 25, 'cluster_110', 21, 'denovo_2', 5, 'cluster_011', 1]}
    """
    seq_a_name = list(seq_a.keys())[0]
    seq_b_name = list(seq_b.keys())[0]
    seq_a = seq_a[seq_a_name][:]
    seq_b = seq_b[seq_b_name][:]
    cluster_a = cluster_a[seq_a_name][:]
    cluster_b = cluster_b[seq_b_name][:]
    # 初始化DP矩阵
    m, n = len(seq_a), len(seq_b)
    dp = np.full((m + 1, n + 1), -np.inf)
    dp[0, 0] = 0
    # 初始化边界条件优化
    dp[1:, 0] = np.arange(-gap_penalty, -gap_penalty * (m + 1), -gap_penalty)
    dp[0, 1:] = np.arange(-gap_penalty, -gap_penalty * (n + 1), -gap_penalty)
    # 优化点2：预先计算序列类型（间隔/Motif）
    type_a = [i % 2 for i in range(1, m + 1)]
    type_b = [j % 2 for j in range(1, n + 1)]
    # 对于迭代时的比对，可能出现motif连续出现的情况
    # 优化点3：向量化运算替换双重循环
    for i in range(1, m + 1):
        ti = type_a[i - 1]  # 提前计算类型
        for j in range(1, n + 1):
            tj = type_b[j - 1]
            # 合并分支判断逻辑
            if ti == 0 and tj == 0:  # Motif匹配
                # Motif簇匹配得分
                if seq_a[i - 1] == seq_b[j - 1]:
                    cluster_match = match_score
                elif cluster_a[i - 1] == cluster_b[j - 1]:
                    cluster_match = cluster_match_score
                else:
                    cluster_match = -mismatch_penalty
                dp[i][j] = max(
                    dp[i - 1][j - 1] + cluster_match,
                    dp[i - 1][j] - gap_penalty,
                    dp[i][j - 1] - gap_penalty,
                )
            elif ti == 1 and tj == 1:  # 间隔匹配
                # 间隔长度差异惩罚（使用绝对差）
                # 后续合并比对时会出现spacer与‘-’对应的情况，减法会报错，把‘-’替换为0
                spacer1 = seq_a[i - 1] if seq_a[i - 1] != "-" else 0
                spacer2 = seq_b[j - 1] if seq_b[j - 1] != "-" else 0
                spacer_penalty = (
                    -same_spacer_score
                    if spacer1 == spacer2
                    else (
                        -sim_spacer_score
                        if (spacer1 - spacer2) <= 40
                        else abs(spacer1 - spacer2) * spacer_coeff
                    )
                )
                dp[i][j] = max(
                    dp[i - 1][j - 1] - spacer_penalty,
                    dp[i - 1][j] - gap_penalty,
                    dp[i][j - 1] - gap_penalty,
                )
            else:  # 混合情况
                dp[i][j] = (
                    max(dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]) - gap_penalty
                )
    # 优化点4：使用预分配内存替代insert(0)
    # 添加最大得分计算
    max_score = dp[m, n]  # 从numpy数组获取最终得分
    # 初始化比对结果存储
    alignment_a = []
    alignment_b = []
    cluster_alignment_a = []
    cluster_alignment_b = []
    operations = []
    i, j = m, n
    # 改为append+reverse模式
    while i > 0 or j > 0:
        current = dp[i][j]
        # 处理匹配情况 (Motif或间隔)
        if i > 0 and j > 0:
            # Motif匹配处理
            if type_a[i - 1] == 0 and type_b[j - 1] == 0:
                if current == dp[i - 1][j - 1] + (
                    match_score
                    if seq_a[i - 1] == seq_b[j - 1]
                    else (
                        cluster_match_score
                        if cluster_a[i - 1] == cluster_b[j - 1]
                        else -mismatch_penalty
                    )
                ):
                    alignment_a.append(seq_a[i - 1])
                    alignment_b.append(seq_b[j - 1])
                    cluster_alignment_a.append(cluster_a[i - 1])
                    cluster_alignment_b.append(cluster_b[j - 1])
                    operations.append("M" if seq_a[i - 1] == seq_b[j - 1] else "C")
                    i -= 1
                    j -= 1
                    continue
            # 间隔匹配处理
            if type_a[i - 1] == 1 and type_b[j - 1] == 1:
                spacer1 = seq_a[i - 1] if seq_a[i - 1] != "-" else 0
                spacer2 = seq_b[j - 1] if seq_b[j - 1] != "-" else 0
                if current == dp[i - 1][j - 1] - (
                    -same_spacer_score
                    if spacer1 == spacer2
                    else (
                        -sim_spacer_score
                        if (spacer1 - spacer2) <= 40
                        else abs(spacer1 - spacer2) * spacer_coeff
                    )
                ):
                    alignment_a.append(seq_a[i - 1])
                    alignment_b.append(seq_b[j - 1])
                    cluster_alignment_a.append(cluster_a[i - 1])
                    cluster_alignment_b.append(cluster_b[j - 1])
                    operations.append("S")
                    i -= 1
                    j -= 1
                    continue
        # 处理gap情况
        if i > 0 and (j == 0 or current == dp[i - 1][j] - gap_penalty):
            alignment_a.append(seq_a[i - 1])
            alignment_b.append("-")
            cluster_alignment_a.append(cluster_a[i - 1])
            cluster_alignment_b.append("-")
            operations.append("D")
            i -= 1
        elif j > 0:
            alignment_a.append("-")
            alignment_b.append(seq_b[j - 1])
            cluster_alignment_a.append("-")
            cluster_alignment_b.append(cluster_b[j - 1])
            operations.append("I")
            j -= 1
    # 反转所有结果列表
    return {
        "alignments": {seq_a_name: alignment_a[::-1], seq_b_name: alignment_b[::-1]},
        "alignment_a": {seq_a_name: alignment_a[::-1]},
        "alignment_b": {seq_b_name: alignment_b[::-1]},
        "cluster_alignments": {
            seq_a_name: cluster_alignment_a[::-1],
            seq_b_name: cluster_alignment_b[::-1],
        },
        "cluster_alignment_a": {seq_a_name: cluster_alignment_a[::-1]},
        "cluster_alignment_b": {seq_b_name: cluster_alignment_b[::-1]},
        "operations": operations[::-1],
        "score": max_score,
        # "matrix": dp,
    }


# 新增批处理函数
@logger.catch
def calculate_distance_pair(param_chunk, motif_features, cluster_features):
    results = []
    for i, j in param_chunk:
        dist = calculate_distance(
            motif_features[i],
            motif_features[j],
            cluster_features[i],
            cluster_features[j],
        )
        results.append((i, j, dist))
    return results


# 新增块计算函数
@logger.catch
def calculate_distance_chunk(i, j, motif_i, motif_j, cluster_i, cluster_j):
    """计算块距离并返回索引"""
    return (i, j, calculate_distance(motif_i, motif_j, cluster_i, cluster_j))


@logger.catch
def process_chunk(p, motif_features, cluster_features):
    """Top-level function to process chunk of pairs"""
    return [
        (
            i,
            j,
            calculate_distance(
                motif_features[i],
                motif_features[j],
                cluster_features[i],
                cluster_features[j],
            ),
        )
        for i, j in p
    ]


# 辅助函数：计算序列间距离矩阵
# 改进1：快速k-mer距离计算（替代全局比对）
@logger.catch
def fast_kmer_distance(seq, k=5):
    """生成k-mer特征向量"""
    # kmer_counts = defaultdict(int)
    # for i in range(len(seq) - k + 1):
    #     kmer = tuple(seq[i : i + k])
    #     kmer_counts[kmer] += 1
    # return kmer_counts
    """使用内置整数作为位掩码"""
    mask = 0
    for i in range(len(seq) - k + 1):
        kmer = tuple(seq[i : i + k])
        # 使用质数哈希减少碰撞
        h = hash(kmer)
        # 取后20位 (可调整)
        pos = h & 0xFFFFF
        mask |= 1 << (pos % 64)  # 直接返回整数
    return mask  # 移除bin()转换


# 优化点：合并特征计算，减少进程间通信
@logger.catch
def process_sequence(args):
    """同时处理序列和聚类数据"""
    seq, clust = args
    return (fast_kmer_distance(seq), fast_kmer_distance(clust))


# 新增辅助函数（放在文件顶部）
def calculate_q_row(args):
    """并行计算Q矩阵的每一行"""
    distance_matrix, total_dist, n, i = args
    row = []
    for j in range(n):
        # 每处理10000列打印一次进度
        if j % 10000 == 0 and j > 0:
            logger.info(f"Processing row {i}: {j}/{n} columns. {mem_monitor()}")
        row.append((n - 2) * distance_matrix[i][j] - total_dist[i] - total_dist[j])
    return row


@logger.catch
def calculate_new_distances(args):
    """并行计算新节点到其他节点的距离"""
    k, i, j, current_distance, distance_matrix = args
    return (k, (distance_matrix[i][k] + distance_matrix[j][k] - current_distance) / 2)


# 新增5：合并序列到比对的函数
@logger.catch
def build_consensus(profile):
    """构建多数表决的一致性序列"""
    consensus = []
    for col in zip(*profile.values()):
        counts = defaultdict(int)
        for elem in col:
            if elem != "-":
                counts[elem] += 1
        consensus.append(max(counts, key=counts.get) if counts else "-")
    return consensus


@logger.catch
def build_cluster_consensus(profile):
    """构建cluster信息的一致性序列"""
    consensus = []
    for col in zip(*profile.values()):
        counts = defaultdict(int)
        for elem in col:
            if elem != "-":
                counts[elem] += 1
        consensus.append(max(counts, key=counts.get) if counts else "-")
    return consensus


# 新增4：比对一致性计算函数
@logger.catch
def calculate_cohesion(alignment):
    """计算比对的簇一致性分数"""
    score = 0
    # 遍历每个比对位置
    for col in zip(*alignment["alignments"].values()):
        # 遍历所有序列对
        for i in range(len(col)):
            for j in range(i + 1, len(col)):
                elem1, elem2 = col[i], col[j]
                # 跳过gap比较
                if elem1 == "-" or elem2 == "-":
                    continue
                # 相同簇加分
                if (
                    alignment["cluster_alignments"][i]
                    == alignment["cluster_alignments"][j]
                ):
                    score += 2
                # 相似间隔加分（差异<40）
                elif isinstance(elem1, int) and isinstance(elem2, int):
                    if abs(elem1 - elem2) <= 40:
                        score += 1
    return score


@logger.catch
def process_iteration(align1, align2, cluster1, cluster2, new_key, kwargs):
    # 生成一致性序列作为代表
    # cons1 = build_consensus(profile1["alignments"][group1_id])
    # cons2 = build_consensus(profile2["alignments"][group2_id])
    # cons1_cluster = build_cluster_consensus(
    #     profile1["cluster_alignments"][group1_id]
    # )
    # cons2_cluster = build_cluster_consensus(
    #     profile2["cluster_alignments"][group2_id]
    # )
    cons1 = build_consensus(align1)
    cons2 = build_consensus(align2)
    cons1_cluster = build_cluster_consensus(cluster1)
    cons2_cluster = build_cluster_consensus(cluster2)
    # 执行全局比对获取操作路径
    global_align = motif_global_alignment(
        {"cons1": cons1},
        {"cons2": cons2},
        {"cons1": cons1_cluster},
        {"cons2": cons2_cluster},
        **kwargs,
    )
    # 构建合并框架
    merged = {
        "alignments": {new_key: defaultdict(list)},
        "cluster_alignments": {new_key: defaultdict(list)},
        "operations": {new_key: global_align["operations"]},
    }
    # 修复点2：正确获取序列名称
    seq_a_name = list(global_align["alignment_a"].keys())[0]
    seq_b_name = list(global_align["alignment_b"].keys())[0]
    # 根据操作路径合并列
    p1_idx = p2_idx = 0
    for op in global_align["operations"]:
        # 处理profile1的列
        if op in ["M", "C", "S", "D"]:
            for seq in align1:
                # 修复点3：使用正确的序列索引
                merged["alignments"][new_key][seq].append(
                    align1[seq][p1_idx] if p1_idx < len(align1[seq]) else "-"
                )
                merged["cluster_alignments"][new_key][seq].append(
                    cluster1[seq][p1_idx] if p1_idx < len(cluster1[seq]) else "-"
                )
            p1_idx += 1
        else:  # Insert gap到profile1
            for seq in align1:
                merged["alignments"][new_key][seq].append("-")
                merged["cluster_alignments"][new_key][seq].append("-")
        # 处理profile2的列
        if op in ["M", "C", "S", "I"]:
            for seq in align2:
                # 修复点4：使用正确的序列索引
                merged["alignments"][new_key][seq].append(
                    align2[seq][p2_idx] if p2_idx < len(align2[seq]) else "-"
                )
                merged["cluster_alignments"][new_key][seq].append(
                    cluster2[seq][p2_idx] if p2_idx < len(cluster2[seq]) else "-"
                )
            p2_idx += 1
        else:  # Insert gap到profile2
            for seq in align2:
                merged["alignments"][new_key][seq].append("-")
                merged["cluster_alignments"][new_key][seq].append("-")
    current_score = calculate_cohesion2(merged)
    return merged, current_score
    # # 评估合并质量
    # current_score = calculate_cohesion(merged)
    # if current_score > best_score:
    #     best_merged = merged
    #     best_score = current_score


@logger.catch
def progressive_motif_alignment(sequences, clusters, **kwargs):
    """
    基于引导树的多序列比对算法
    参数:
        sequences: 列表，包含多个待比对序列
        clusters: 列表，包含对应的cluster信息
        **kwargs: 传递给motif_global_alignment的参数
    返回:
        字典，包含对齐后的序列和操作记录
    """
    # parser = argparse.ArgumentParser(description="Enhancer Alignment Parameters")
    # parser.add_argument(
    #     "--load-existing", action="store_true", help="Load precomputed matrix and tree"
    # )
    # parser.add_argument(
    #     "--save-dir",
    #     type=str,
    #     default="./cache",
    #     help="Directory for saving computed data",
    # )
    # args = parser.parse_args()

    # 创建缓存目录
    os.makedirs(args.save_dir, exist_ok=True)

    # 生成唯一文件标识
    data_hash = hashlib.md5(str(sequences.keys()).encode()).hexdigest()[:8]
    matrix_file = os.path.join(args.save_dir, f"dist_matrix_{data_hash}.pkl")
    matrix_file2 = os.path.join(args.save_dir, f"dist_matrix2_{data_hash}.pkl")
    tree_file = os.path.join(args.save_dir, f"nj_tree_{data_hash}.pkl")

    checkpoint_file = os.path.join(args.save_dir, f"checkpoint_{data_hash}.pkl")
    checkpoint_interval = 2 * 3600  # 6小时转换为秒
    last_save_time = time.time()

    @logger.catch
    def get_rev_sequence(seq):
        """生成反向序列（根据实际数据结构调整）"""
        return seq[::-1]  # 假设序列是列表结构，直接反转

    # 将嵌套函数移至顶层或改为静态方法
    @staticmethod
    @lru_cache(maxsize=None)
    @logger.catch
    # Add at top level (not inside any function)
    def kmer_distance(i, j, a, b, b_rev, a_cluster, b_cluster, b_rev_cluster):
        """计算k-mer相似性"""
        """Top-level function for pickle compatibility"""
        # i,j 是两个序列的编码
        # a = frozenset(seqs_values[i])
        # b = frozenset(seqs_values[j])
        # b_rev = frozenset(get_rev_sequence(seqs_values[j]))
        # a_cluster = frozenset(clusters_values[i])
        # b_cluster = frozenset(clusters_values[j])
        # b_rev_cluster = frozenset(get_rev_sequence(clusters_values[j]))
        # 计算k-mer距离，30% k-mer相似性和70% cluster相似性
        forward_kmer_dis = 0.3 * (1 - (len(a & b) / max(len(a | b), 1))) + 0.7 * (
            len(a_cluster & b_cluster) / max(len(a_cluster | b_cluster), 1)
        )
        reverse_kmer_dis = 0.3 * (
            1 - (len(a & b_rev) / max(len(a | b_rev), 1))
        ) + 0.7 * (
            len(a_cluster & b_rev_cluster) / max(len(a_cluster | b_rev_cluster), 1)
        )
        # 选择较小的距离
        if forward_kmer_dis <= reverse_kmer_dis:
            return forward_kmer_dis, True
        else:
            return reverse_kmer_dis, False

    # # Add at top level (not inside any function)
    # def kmer_similarity(a, b):
    #     """Top-level function for pickle compatibility"""
    #     return len(a & b) / max(len(a | b), 1)

    @logger.catch
    def build_distance_matrix(seqs, clusters):
        """MAFFT风格快速距离估算（优化版）"""
        n = len(seqs)

        # 预计算k-mer特征（多进程优化）
        logger.info("Precomputing k-mer features...")

        # 将生成器转为列表避免重复迭代
        seq_pairs = list(zip(seqs.values(), clusters.values()))

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=min(os.cpu_count() // 2, 40)
        ) as executor:
            # 动态调整chunksize计算公式：sqrt(N)/(4*cpu_count)
            chunk_size = max(
                1, int(len(seq_pairs) ** 0.5 / (os.cpu_count() * 4)) or 1000
            )

            # 单次map调用处理两个特征
            features = list(
                executor.map(process_sequence, seq_pairs, chunksize=chunk_size)
            )

            # 解包结果
            motif_features, cluster_features = zip(*features)

        matrix = np.zeros((n, n), dtype=np.float32)
        # 优化点：改用ThreadPool减少内存复制
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=min(os.cpu_count() // 2, 40)
        ) as executor:
            # 使用生成器代替列表存储参数
            params = ((i, j) for i in range(n) for j in range(i + 1, n))

            # 分块处理减少内存占用
            chunk_size = min(1000, n * (n - 1) // 2 // (os.cpu_count() * 2))
            chunk_size = max(chunk_size, 1000)  # 确保至少为1
            futures = [
                executor.submit(
                    process_chunk,  # 使用顶层函数替代lambda
                    list(itertools.islice(params, chunk_size)),
                    motif_features,
                    cluster_features,
                )
                for _ in range(0, n * (n - 1) // 2, chunk_size)
            ]
            # 流式处理结果
            processed_pairs_num = 0
            total_pairs = n * (n - 1) // 2
            for future in concurrent.futures.as_completed(futures):
                for i, j, dist in future.result():
                    matrix[i][j] = matrix[j][i] = dist
                    processed_pairs_num += 1
                    if processed_pairs_num % 10000 == 0:
                        logger.info(
                            f"Processed {processed_pairs_num}/{total_pairs} pairs. \n{mem_monitor()}"
                        )
            # 最终进度报告
            logger.info(
                f"Processed {processed_pairs_num}/{total_pairs} pairs. \n{mem_monitor()}"
            )

        np.fill_diagonal(matrix, 1e-6)
        return matrix

    @logger.catch
    def upgma_tree(distance_matrix, seqs=sequences):
        """UPGMA树构建实现"""
        n = distance_matrix.shape[0]
        # 初始化每个节点为一个簇
        tree_clusters = [
            {
                "id": list(seqs.keys())[i],
                "base_id": list(seqs.keys())[i].replace("_rev", ""),
                "size": 1,
                "height": 0,
            }
            for i in range(n)
        ]
        # current_dist = distance_matrix.copy()
        distance_dict = defaultdict(float)
        distance_dict = defaultdict(float)
        # 初始填充（在构建初始堆之前添加）
        for i in range(n):
            for j in range(i + 1, n):
                distance_dict[(i, j)] = distance_dict[(j, i)] = distance_matrix[i][j]
        # 初始化优先队列（堆优化关键点）
        heap = []
        # 预计算所有初始距离对
        for i in range(n):
            for j in range(i + 1, n):
                heappush(heap, (distance_matrix[i][j], i, j))
        # 维护有效索引的集合
        valid_indices = set(range(n))
        # 索引映射表（跟踪簇的位置变化）
        index_map = {i: i for i in range(n)}
        # 把每个序列的base_id存储到一个数组中，方便查找
        base_id_array = np.array(
            [cluster["base_id"] for cluster in tree_clusters]
        )  # 转换为numpy数组
        tree = []
        existing_bases = set()
        existing_ids = set()
        while len(seq_names) > 0:
            # 获取当前最小距离对（堆优化核心）
            valid_pair_found = False
            while heap:
                min_dist, i, j = heappop(heap)
                logger.info(
                    f"Processing pair: {tree_clusters[i]['id']}, {tree_clusters[j]['id']} with distance: {min_dist}"
                )
                # 验证索引有效性
                if i in valid_indices and j in valid_indices:
                    valid_pair_found = True
                    break
            if valid_pair_found:
                # 修复点1：使用index_map获取实际簇位置
                real_i = index_map[i]
                real_j = index_map[j]
                # 检查base_id是否已存在
                a_base = tree_clusters[real_i]["base_id"]
                b_base = tree_clusters[real_j]["base_id"]
                # 记录合并并更新existing_bases
                tree.append((tree_clusters[i]["id"], tree_clusters[j]["id"]))
                if isinstance(a_base, str):
                    existing_bases.add(a_base)
                if isinstance(b_base, str):
                    existing_bases.add(b_base)
                if isinstance(tree_clusters[i]["id"], str):
                    existing_ids.add(tree_clusters[i]["id"])
                if isinstance(tree_clusters[j]["id"], str):
                    existing_ids.add(tree_clusters[j]["id"])
                new_cluster = {
                    "id": (
                        tree_clusters[index_map[i]]["id"],
                        tree_clusters[index_map[j]]["id"],
                    ),
                    "size": tree_clusters[index_map[i]]["size"]
                    + tree_clusters[index_map[j]]["size"],
                    "height": min_dist / 2,
                    "base_id": (a_base, b_base),
                }
                # 更新索引映射和有效集合
                new_index = len(tree_clusters)
                if isinstance(a_base, str):
                    for t in np.where(base_id_array == a_base)[0]:
                        valid_indices.remove(t)
                        seq_names.discard(list(seqs.keys())[t])
                if isinstance(b_base, str):
                    for t in np.where(base_id_array == b_base)[0]:
                        valid_indices.remove(t)
                        seq_names.discard(list(seqs.keys())[t])
                valid_indices.add(new_index)
                index_map[new_index] = len(tree_clusters)
                # 计算新簇与其他簇的距离
                new_distances = {}
                for k in valid_indices - {new_index}:
                    orig_k = index_map[k]
                    # 从字典获取原始距离
                    dist_i = distance_dict.get((i, orig_k), 0)
                    dist_j = distance_dict.get((j, orig_k), 0)
                    new_dist = (dist_i * size_i + dist_j * size_j) / (size_i + size_j)
                    new_distances[k] = new_dist
                    distance_dict[(new_index, k)] = distance_dict[(k, new_index)] = (
                        new_dist
                    )
                # 删除不再需要的旧距离（添加在合并后）
                for k in [i, j]:
                    for m in distance_dict.copy():
                        if k in m:
                            del distance_dict[m]
                # 更新堆推送逻辑
                for k, dist in new_distances.items():
                    heappush(heap, (dist, new_index, k))
                # new_size = len(tree_clusters) - 1
                # new_dist = np.zeros((new_size, new_size))
                # row = 0
                # for k in range(len(tree_clusters)):
                #     if k != i and k != j:
                #         new_dist[row, :-1] = (
                #             current_dist[i][k] * tree_clusters[i]["size"]
                #             + current_dist[j][k] * tree_clusters[j]["size"]
                #         ) / (tree_clusters[i]["size"] + tree_clusters[j]["size"])
                #         row += 1
                # current_dist = new_dist
                # # 计算新簇与其他簇的距离并推入堆
                # for k in valid_indices - {new_index}:
                #     orig_k = index_map[k]
                #     # 计算新距离（UPGMA公式）
                #     size_i = tree_clusters[index_map[i]]["size"]
                #     size_j = tree_clusters[index_map[j]]["size"]
                #     new_dist = (
                #         current_dist[i][k] * size_i + current_dist[j][k] * size_j
                #     ) / (size_i + size_j)
                #     heappush(heap, (new_dist, new_index, k))
                tree_clusters.append(new_cluster)
        return tree, existing_ids

    @logger.catch
    def neighbor_joining(distance_matrix, seqs):
        """Neighbor-Joining算法实现 (多核优化版)"""
        n = distance_matrix.shape[0]
        active_nodes = [
            {
                "id": name,
                "age": 0,
                "base_id": name.replace("_rev", ""),  # 新增base_id字段
            }
            for name in seqs.keys()
        ]
        existing_ids = set()
        processed_nodes = 0  # 新增进度计数器
        # 预创建进程池，避免频繁创建销毁
        max_workers = min(os.cpu_count(), 80)
        process_pool = concurrent.futures.ProcessPoolExecutor(max_workers=max_workers)
        # 每5000次合并打印进度到日志，包括已合并/总节点数，以及内存占用
        valid_nodes = np.arange(n)  # 初始化有效节点索引
        while len(active_nodes) > 1:
            if processed_nodes % 5000 == 0:
                logger.info(
                    f"merged NJ {processed_nodes} of {n} nodes. {processed_nodes/n*100:.2f}% \n{mem_monitor()}"
                )
            current_n = distance_matrix.shape[0]  # 获取当前矩阵的实际维度
            # 并行计算Q矩阵
            # with concurrent.futures.ProcessPoolExecutor(
            #     max_workers=min(os.cpu_count() // 2, 40)
            # ) as executor:
            total_dist = distance_matrix.sum(axis=1)
            args_list = [
                (distance_matrix, total_dist, current_n, i) for i in range(current_n)
            ]
            # 优化点1：动态调整chunksize计算公式
            chunk_size_q = max(
                1, len(args_list) // (max_workers * 2)
            )  # 使用平方根调整块大小

            # 原代码
            # chunk_size = max(1, int(current_n / os.cpu_count()))
            # q_rows = list(executor.map(calculate_q_row, args_list, chunksize=chunk_size))

            # 优化后
            # with concurrent.futures.ProcessPoolExecutor(
            #     max_workers=min(os.cpu_count(), 80)
            # ) as executor:
            # 预生成参数列表
            # args_list = [
            #     (distance_matrix, total_dist, current_n, i)
            #     for i in range(current_n)
            # ]
            # 动态chunk_size计算
            # chunk_size_q = max(1, len(args_list) // (os.cpu_count() * 2))
            q_rows = list(
                process_pool.map(calculate_q_row, args_list, chunksize=chunk_size_q)
            )

            q_matrix = np.array(q_rows)
            np.fill_diagonal(q_matrix, np.inf)

            i, j = np.unravel_index(np.argmin(q_matrix), q_matrix.shape)
            if "id" in active_nodes[i]:
                existing_ids.add(active_nodes[i]["id"])
            if "id" in active_nodes[j]:
                existing_ids.add(active_nodes[j]["id"])
            # 创建新节点前计算必要参数
            delta = (total_dist[i] - total_dist[j]) / max((n - 2), 1)
            current_distance = distance_matrix[i, j]
            limb_i = (current_distance + delta) / 2
            limb_j = (current_distance - delta) / 2

            # 并行计算新距离
            # 优化点2：批量预生成参数
            valid_indices = np.arange(current_n)
            args_list = [
                (k, i, j, current_distance, distance_matrix) for k in valid_indices
            ]
            # 优化chunk_size计算
            chunk_size_dist = max(1, len(args_list) // (max_workers * 2))
            # with concurrent.futures.ProcessPoolExecutor(
            #     max_workers=min(os.cpu_count() // 2, 40)
            # ) as executor:
            # 使用实际有效的节点索引范围
            # valid_indices = np.arange(current_n)
            # args_list = [
            #     (k, i, j, current_distance, distance_matrix) for k in valid_indices
            # ]
            # chunk_size = max(1, int(len(args_list) / os.cpu_count()))
            new_dist_results = list(
                process_pool.map(
                    calculate_new_distances, args_list, chunksize=chunk_size_dist
                )
            )

            # 重构新距离向量
            new_dist = np.zeros(len(valid_indices))
            for k, dist in new_dist_results:
                new_dist[k] = dist

            # 处理mask和节点过滤
            merged_base_ids = {active_nodes[i]["base_id"], active_nodes[j]["base_id"]}
            exclude_indices = [
                idx
                for idx, node in enumerate(active_nodes)
                if node["base_id"] in merged_base_ids
            ]
            mask = np.ones(len(active_nodes), bool)
            mask[exclude_indices] = False
            valid_nodes = np.where(mask)[0]

            # 更新距离矩阵（使用原始矩阵的索引）
            distance_matrix = distance_matrix[np.ix_(mask, mask)]  # 先过滤现有节点
            new_dist = new_dist[mask]  # 再过滤新距离向量

            # 矩阵扩展
            if distance_matrix.size > 0:
                distance_matrix = np.vstack(
                    [
                        np.column_stack([distance_matrix, new_dist]),
                        np.append(new_dist, 0),
                    ]
                )

            # 创建新节点并更新列表
            new_node = {
                "left": active_nodes[i],
                "right": active_nodes[j],
                "age": current_distance / 2,
                "branch_i": limb_i,
                "branch_j": limb_j,
                "base_id": f"{active_nodes[i]['base_id']}|{active_nodes[j]['base_id']}",
            }
            active_nodes = [
                node for idx, node in enumerate(active_nodes) if mask[idx]
            ] + [new_node]
            n -= 1
            processed_nodes += 1  # 更新计数器
        process_pool.shutdown(wait=True)  # 关闭进程池
        return active_nodes[0], existing_ids

    @logger.catch
    def count_forward(names_chunk, labels_chunk):
        """并行统计正反向序列数"""
        cnt = [0, 0]
        for name, label in zip(names_chunk, labels_chunk):
            if "_rev" not in name:
                cnt[label] += 1
        return cnt[0], cnt[1]

    @logger.catch
    def cluster_and_filter_matrix(distance_matrix, seqs):
        # 新增内存映射优化
        if not distance_matrix.flags["C_CONTIGUOUS"]:
            distance_matrix = np.ascontiguousarray(distance_matrix)
        # 方案1：修改距离矩阵法（适用于任意聚类算法）
        # 生成base_id映射表
        seq_names = np.array(list(seqs.keys()))
        base_ids = np.array(
            [name.split("_rev")[0] for name in seq_names]
        )  # 提取base_id

        # 创建相同base_id的掩码矩阵
        same_base_mask = base_ids[:, None] == base_ids[None, :]
        np.fill_diagonal(same_base_mask, False)  # 排除对角线

        # 给相同base_id的序列设置极大距离
        constrained_matrix = distance_matrix.copy()
        constrained_matrix[same_base_mask] = 1e19

        # K-中心点优化版
        from sklearn_extra.cluster import KMedoids

        cluster = KMedoids(
            n_clusters=2,
            metric="precomputed",
            init="k-medoids++",
            max_iter=5,  # 限制迭代次数
        )

        # 随机采样1%数据初始化
        sample_idx = np.random.choice(len(distance_matrix), size=680, replace=False)
        cluster.fit(distance_matrix[sample_idx][:, sample_idx])
        # labels = cluster.fit_predict(distance_matrix)
        labels = cluster.fit_predict(constrained_matrix)  # 使用修改后的矩阵
        # 2. 并行统计正反向序列数
        seq_names = np.array(list(seqs.keys()))
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # 分块统计
            chunk_size = 5000
            futures = []
            for i in range(0, len(seq_names), chunk_size):
                chunk = seq_names[i : i + chunk_size]
                futures.append(
                    executor.submit(count_forward, chunk, labels[i : i + chunk_size])
                )

            # 合并统计结果
            counts = [0, 0]
            for future in concurrent.futures.as_completed(futures):
                c0, c1 = future.result()
                counts[0] += c0
                counts[1] += c1

        # 3. 选择正向较多的簇
        selected_label = 0 if counts[0] > counts[1] else 1
        selected_indices = np.where(labels == selected_label)[0]

        # 4. 过滤输出
        filtered_matrix = distance_matrix[selected_indices][:, selected_indices]
        filtered_seqs = {seq_names[i]: seqs[seq_names[i]] for i in selected_indices}

        return filtered_matrix, filtered_seqs

    @logger.catch
    def cluster_and_filter_matrix2(distance_matrix, seqs):
        # 新增内存映射优化
        if not distance_matrix.flags["C_CONTIGUOUS"]:
            distance_matrix = np.ascontiguousarray(distance_matrix)
        seq_names = np.array(list(seqs.keys()))
        base_ids = np.array(
            [name.split("_rev")[0] for name in seq_names]
        )  # 提取base_id

        # 方案2：聚类后过滤法（更高效）
        def enforce_base_separation(labels):
            from collections import defaultdict

            cluster_map = defaultdict(list)
            # 建立base_id到簇的映射
            for idx, label in enumerate(labels):
                base = base_ids[idx]
                cluster_map[(base, label)].append(idx)

            # 分离冲突的base_id
            new_labels = labels.copy()
            for (base, cluster), indices in cluster_map.items():
                if len(indices) > 1:  # 发现冲突
                    # 将后半部分分配到另一个簇
                    for i in indices[1:]:
                        new_labels[i] = 1 - cluster
            return new_labels

        # K-中心点优化版
        from sklearn_extra.cluster import KMedoids

        cluster = KMedoids(
            n_clusters=2,
            metric="precomputed",
            init="k-medoids++",
            max_iter=5,  # 限制迭代次数
        )

        # 随机采样1%数据初始化
        sample_idx = np.random.choice(len(distance_matrix), size=680, replace=False)
        cluster.fit(distance_matrix[sample_idx][:, sample_idx])
        # labels = cluster.fit_predict(distance_matrix)
        # 或者使用方案2的聚类后处理：
        labels = cluster.fit_predict(distance_matrix)
        labels = enforce_base_separation(labels)
        # 2. 并行统计正反向序列数
        seq_names = np.array(list(seqs.keys()))
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # 分块统计
            chunk_size = 5000
            futures = []
            for i in range(0, len(seq_names), chunk_size):
                chunk = seq_names[i : i + chunk_size]
                futures.append(
                    executor.submit(count_forward, chunk, labels[i : i + chunk_size])
                )

            # 合并统计结果
            counts = [0, 0]
            for future in concurrent.futures.as_completed(futures):
                c0, c1 = future.result()
                counts[0] += c0
                counts[1] += c1

        # 3. 选择正向较多的簇
        selected_label = 0 if counts[0] > counts[1] else 1
        selected_indices = np.where(labels == selected_label)[0]

        # 4. 过滤输出
        filtered_matrix = distance_matrix[selected_indices][:, selected_indices]
        filtered_seqs = {seq_names[i]: seqs[seq_names[i]] for i in selected_indices}

        return filtered_matrix, filtered_seqs

    @logger.catch
    def distance_matrix_to_raxml(distance_matrix, output_path):
        """将距离矩阵转换为RAxML兼容的PHYLIP格式"""
        n = distance_matrix.shape[0]
        with open(output_path, "w") as f:
            f.write(f"{n}\n")
            # 写入完整方阵
            for i in range(n):
                # 写入序列标识符（使用实际序列名替换i）
                row = " ".join(f"{distance_matrix[i,j]:.4f}" for j in range(n))
                f.write(f"{i}    {row}\n")

    @logger.catch
    def neighbor_joining_fastme(distance_matrix, seqs):
        """使用FastME算法替代NJ (性能提升10倍以上)"""
        # 生成PHYLIP格式的距离矩阵
        phylip_file = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/fastme/temp.phylip"
        n = distance_matrix.shape[0]
        with open(phylip_file, "w") as f:
            f.write(f"{n}\n")
            np.savetxt(f, distance_matrix, fmt="%.4f")

        # 调用预编译的FastME二进制文件（需提前安装）
        # 下载：http://www.atgc-montpellier.fr/fastme/
        cmd = f"fastme -i {phylip_file} -o /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/cache_fastme/output_tree.nwk -n -s"  # 并行4线程，使用SPR优化
        subprocess.run(cmd, shell=True, check=True)

        # 解析生成的newick树
        with open(
            "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/cache_fastme/output_tree.nwk"
        ) as f:
            newick_tree = f.read()

        # 转换为原有数据结构（需实现转换函数）
        return newick_to_custom_format(newick_tree), set(seqs.keys())

    @logger.catch
    def rapid_nj(distance_matrix, seqs):
        """近似邻接法实现"""
        # from rapidnj import build_tree
        # 生成PHYLIP格式的下三角矩阵
        phylip_file = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/fastme/temp.phylip"
        n = distance_matrix.shape[0]
        with open(phylip_file, "w") as f:
            f.write(f"{n}\n")
            for i in range(n):
                row = " ".join(f"{distance_matrix[i,j]:.4f}" for j in range(i + 1))
                f.write(f"{i}    {row}\n")

        # 调用预编译的rapidnj二进制
        cmd = f"rapidnj {phylip_file} -i pd -o t -c 64 > output.nwk"
        subprocess.run(cmd, shell=True, check=True)

        # 解析结果树
        with open("output.nwk") as f:
            return parse_newick_to_guide(newick_tree), set(seqs.keys())

    @logger.catch
    def convert_nj_to_guide(nj_tree):
        """将NJ树转换为UPGMA兼容的引导树结构"""
        guide_tree = []

        def traverse(node):
            if "id" in node:  # 叶子节点
                return node["id"], []

            # 递归处理左右子树
            left_id, left_pairs = traverse(node["left"])
            right_id, right_pairs = traverse(node["right"])

            # 创建合并对并生成新簇ID
            new_pair = (left_id, right_id)
            guide_tree.append(new_pair)

            # 生成新簇ID格式与UPGMA一致
            # new_cluster_id = f"({left_id},{right_id})"
            # new_cluster_id is a tuple
            new_cluster_id = (left_id, right_id)
            return new_cluster_id, left_pairs + right_pairs + [new_pair]

        # 执行遍历并返回最终结构
        final_id, all_pairs = traverse(nj_tree)
        return all_pairs

    # 把'numpy.ndarray'保存成phylip格式的矩阵
    @logger.catch
    def save_phylip_matrix(matrix, output_path, seqs):
        """保存距离矩阵为PHYLIP格式"""
        n = matrix.shape[0]
        seq_names = list(seqs.keys())
        with open(output_path, "w") as f:
            f.write(f"{n}\n")
            for i in range(n):
                row = " ".join(f"{matrix[i,j]:.4f}" for j in range(n))
                f.write(f"{seq_names[i]}    {row}\n")

    @logger.catch
    def parse_newick_to_guide(newick_str):
        """将RapidNJ生成的Newick树转换为渐进比对用的引导树结构"""
        from ete3 import Tree

        def split_by_commas(s):
            """安全分割嵌套括号字符串"""
            parts = []
            current = []
            level = 0
            for c in s:
                if c == "(":
                    level += 1
                elif c == ")":
                    level -= 1
                elif c == "," and level == 0:
                    parts.append("".join(current).strip())
                    current = []
                    continue
                current.append(c)
            parts.append("".join(current).strip())
            return parts

        def parse_node(label):
            """递归解析节点标签"""
            if not label.startswith("("):
                return label.replace("'", "")

            # 去除外层括号
            content = label[1:-1]
            children = split_by_commas(content)
            return tuple(parse_node(child) for child in children)

        # 处理文件路径或直接字符串
        if os.path.isfile(newick_str):
            with open(newick_str, "r") as f:
                newick_str = f.read().strip()

        t = Tree(newick_str)
        merge_order = []

        for node in t.traverse("postorder"):
            if not node.is_leaf():
                # 解析左右子节点
                left = (
                    parse_node(node.children[0].name)
                    if not node.children[0].is_leaf()
                    else node.children[0].name.replace("'", "")
                )
                right = (
                    parse_node(node.children[1].name)
                    if not node.children[1].is_leaf()
                    else node.children[1].name.replace("'", "")
                )

                merge_order.append((left, right))
                node.name = str((left, right))

        return merge_order
        # for node in t.traverse("postorder"):
        #     if not node.is_leaf():
        #         # 获取左右子节点的名称（叶子节点保留原名，内部节点使用元组）
        #         left = (
        #             node.children[0].name.replace("'", "")
        #             if node.children[0].is_leaf()
        #             else tuple(eval(node.children[0].name))
        #         )
        #         right = (
        #             node.children[1].name.replace("'", "")
        #             if node.children[1].is_leaf()
        #             else tuple(eval(node.children[1].name))
        #         )
        #         merge_order.append((left, right))
        #         node.name = str((left, right))  # 更新内部节点名称

        # return merge_order

    # 改进2：迭代优化的merge_alignments
    @logger.catch
    def merge_alignments(
        profile1, profile2, group1_id, group2_id, max_iter=3, **kwargs
    ):
        """
        MAFFT风格profile比对合并
        既可以合并两个比对，也可以把新序列插入到已有比对中
        # 新序列插入已有比对
        existing_alignment = {
            'alignments': {align_id: {'seq1': [0, 'M1', 5], 'seq2': [0, 'M1', 10]}},
            'cluster_alignments': {align_id: {'seq1': ['C1', 'C2'], 'seq2': ['C1', 'C3']}}
        }
        # 将新序列格式化为profile
        new_sequence = {
            'alignments': {align_id: {'new_seq': [5, 'M1', 8]}},
            'cluster_alignments': {align_id: {'new_seq': ['C1', 'C4']}}
        }
        merged = merge_alignments(existing_alignment, new_sequence)
        """
        new_key = (group1_id, group2_id)
        best_merged = {
            "alignments": {new_key: defaultdict(list)},
            "cluster_alignments": {new_key: defaultdict(list)},
            "operations": {new_key: []},
        }
        best_score = -np.inf
        # 并行处理多轮迭代优化
        with concurrent.futures.ProcessPoolExecutor() as executor:
            # 生成所有迭代参数
            iter_args = [
                (
                    profile1["alignments"][group1_id],
                    profile2["alignments"][group2_id],
                    profile1["cluster_alignments"][group1_id],
                    profile2["cluster_alignments"][group2_id],
                    new_key,
                    kwargs,
                )
                for _ in range(max_iter)
            ]

            # 并行执行迭代
            futures = [executor.submit(process_iteration, *args) for args in iter_args]
            # 收集并比较结果
            for future in concurrent.futures.as_completed(futures):
                merged, current_score = future.result()
                if current_score > best_score:
                    best_merged = merged
                    best_score = current_score
        return best_merged
        # for _ in range(max_iter):
        #     # 生成一致性序列作为代表
        #     cons1 = build_consensus(profile1["alignments"][group1_id])
        #     cons2 = build_consensus(profile2["alignments"][group2_id])
        #     cons1_cluster = build_cluster_consensus(
        #         profile1["cluster_alignments"][group1_id]
        #     )
        #     cons2_cluster = build_cluster_consensus(
        #         profile2["cluster_alignments"][group2_id]
        #     )
        #     # 执行全局比对获取操作路径
        #     global_align = motif_global_alignment(
        #         {"cons1": cons1},
        #         {"cons2": cons2},
        #         {"cons1": cons1_cluster},
        #         {"cons2": cons2_cluster},
        #         **kwargs,
        #     )
        #     # 构建合并框架
        #     merged = {
        #         "alignments": {new_key: defaultdict(list)},
        #         "cluster_alignments": {new_key: defaultdict(list)},
        #         "operations": {new_key: global_align["operations"]},
        #     }
        #     # 修复点2：正确获取序列名称
        #     seq_a_name = list(global_align["alignment_a"].keys())[0]
        #     seq_b_name = list(global_align["alignment_b"].keys())[0]
        #     # 根据操作路径合并列
        #     p1_idx = p2_idx = 0
        #     for op in global_align["operations"]:
        #         # 处理profile1的列
        #         if op in ["M", "C", "S", "D"]:
        #             for seq in profile1["alignments"][group1_id]:
        #                 # 修复点3：使用正确的序列索引
        #                 merged["alignments"][new_key][seq].append(
        #                     profile1["alignments"][group1_id][seq][p1_idx]
        #                     if p1_idx < len(profile1["alignments"][group1_id][seq])
        #                     else "-"
        #                 )
        #                 merged["cluster_alignments"][new_key][seq].append(
        #                     profile1["cluster_alignments"][group1_id][seq][p1_idx]
        #                     if p1_idx
        #                     < len(profile1["cluster_alignments"][group1_id][seq])
        #                     else "-"
        #                 )
        #             p1_idx += 1
        #         else:  # Insert gap到profile1
        #             for seq in profile1["alignments"][group1_id]:
        #                 merged["alignments"][new_key][seq].append("-")
        #                 merged["cluster_alignments"][new_key][seq].append("-")
        #         # 处理profile2的列
        #         if op in ["M", "C", "S", "I"]:
        #             for seq in profile2["alignments"][group2_id]:
        #                 # 修复点4：使用正确的序列索引
        #                 merged["alignments"][new_key][seq].append(
        #                     profile2["alignments"][group2_id][seq][p2_idx]
        #                     if p2_idx < len(profile2["alignments"][group2_id][seq])
        #                     else "-"
        #                 )
        #                 merged["cluster_alignments"][new_key][seq].append(
        #                     profile2["cluster_alignments"][group2_id][seq][p2_idx]
        #                     if p2_idx
        #                     < len(profile2["cluster_alignments"][group2_id][seq])
        #                     else "-"
        #                 )
        #             p2_idx += 1
        #         else:  # Insert gap到profile2
        #             for seq in profile2["alignments"][group2_id]:
        #                 merged["alignments"][new_key][seq].append("-")
        #                 merged["cluster_alignments"][new_key][seq].append("-")
        #     # 评估合并质量
        #     current_score = calculate_cohesion(merged)
        # if current_score > best_score:
        #     best_merged = merged
        #     best_score = current_score
        # return best_merged

    # 新增3：迭代细化函数
    @logger.catch
    def iterative_refinement(alignment, final_id, iterations=2, **kwargs):
        """MAFFT风格迭代优化"""

        # 新增分割函数
        @logger.catch
        def split_alignment(aln, final_id):
            """随机分割比对为两个子组"""
            # 获取字典的所有键并打乱顺序
            keys = list(aln["alignments"][final_id].keys())
            random.shuffle(keys)
            mid = len(keys) // 2
            iter_id1 = tuple(keys[:mid])
            iter_id2 = tuple(keys[mid:])
            # 创建两个子比对profile
            group1 = {
                "alignments": {
                    iter_id1: {k: aln["alignments"][final_id][k] for k in iter_id1}
                },
                "cluster_alignments": {
                    iter_id1: {
                        k: aln["cluster_alignments"][final_id][k] for k in iter_id1
                    }
                },
                "operations": {iter_id1: aln["operations"]},
            }
            group2 = {
                "alignments": {
                    iter_id2: {k: aln["alignments"][final_id][k] for k in keys[mid:]}
                },
                "cluster_alignments": {
                    iter_id2: {
                        k: aln["cluster_alignments"][final_id][k] for k in keys[mid:]
                    }
                },
                "operations": {iter_id2: aln["operations"][final_id]},
            }
            return group1, group2, iter_id1, iter_id2

        for _ in range(iterations):
            # 随机分割比对组进行重新对齐
            prof1, prof2, iter_id1, iter_id2 = split_alignment(alignment, final_id)
            new_alignment = merge_alignments(prof1, prof2, iter_id1, iter_id2, **kwargs)
            if calculate_cohesion2(new_alignment) > calculate_cohesion2(alignment):
                alignment = new_alignment
        return alignment

    # 主流程
    # 1. 构建距离矩阵
    logger.info("initializing distance matrix")
    # dist_matrix = build_distance_matrix(sequences, clusters)
    # 1. 构建/加载距离矩阵
    if args.load_existing and os.path.exists(matrix_file):
        logger.info(f"Loading distance matrix from {matrix_file}")
        with open(matrix_file, "rb") as f:
            dist_matrix = pickle.load(f)
    else:
        logger.info("initializing distance matrix")
        dist_matrix = build_distance_matrix(sequences, clusters)
        # 保存新计算的矩阵
        with open(matrix_file, "wb") as f:
            pickle.dump(dist_matrix, f)
        logger.info(f"Saved distance matrix to {matrix_file}")
    logger.info("distance matrix built")
    # cluster and filter
    # f_dist_matrix, f_sequences = cluster_and_filter_matrix(dist_matrix, sequences)
    if args.load_existing and os.path.exists(matrix_file2):
        logger.info(f"Loading filtered distance matrix from {matrix_file2}")
        with open(matrix_file2, "rb") as f:
            f_dist_matrix2, f_sequences2 = pickle.load(f)
        if not os.path.exists(f"{args.save_dir}/temp_half.phylip"):
            logger.info("transforming filtered distance matrix into temp_half.phylip")
            save_phylip_matrix(
                f_dist_matrix2,
                f"{args.save_dir}/temp_half.phylip",
                f_sequences2,
            )
        if not os.path.exists(f"{args.save_dir}/output_half.nwk"):
            logger.info("rapidnj building tree using temp_half.phylip")
            # 使用nj法建树并转换格式
            half_phylip = f"{args.save_dir}/temp_half.phylip"
            cmd = f"conda run -n main rapidnj {half_phylip} -i pd -o t -c 64 > {args.save_dir}/output_half.nwk"
            subprocess.run(cmd, shell=True, check=True)
    else:
        logger.info("initializing filtered distance matrix")
        # f_dist_matrix, f_sequences = cluster_and_filter_matrix(dist_matrix, sequences)
        # 使用新的聚类和过滤方法
        f_dist_matrix2, f_sequences2 = cluster_and_filter_matrix2(
            dist_matrix.copy(), sequences
        )
        # rapidnj temp_half.phylip -i pd -o t -c 72 > output_half.nwk 2>output_half.err &
        # 保存新计算的矩阵
        with open(matrix_file2, "wb") as f:
            pickle.dump((f_dist_matrix2, f_sequences2), f)
        logger.info(f"Saved filtered distance matrix to {matrix_file2}")
        save_phylip_matrix(
            f_dist_matrix2,
            f"{args.save_dir}/temp_half.phylip",
            f_sequences2,
        )
        if not os.path.exists(f"{args.save_dir}/output_half.nwk"):
            logger.info("rapidnj building tree using temp_half.phylip")
            # 使用nj法建树并转换格式
            half_phylip = f"{args.save_dir}/temp_half.phylip"
            cmd = f"conda run -n main rapidnj {half_phylip} -i pd -o t -c 64 > {args.save_dir}/output_half.nwk"
            subprocess.run(cmd, shell=True, check=True)
    # f_dist_matrix2, f_sequences2 = cluster_and_filter_matrix2(dist_matrix, sequences)
    # 2. 生成引导树
    logger.info("initializing guide tree")
    # guide_tree, existing_ids = upgma_tree(dist_matrix)
    # 使用nj法建树并转换格式
    # nj_tree, existing_ids = neighbor_joining(dist_matrix.copy(), sequences)
    # 2. 生成/加载引导树
    if args.load_existing and os.path.exists(tree_file):
        logger.info(f"Loading NJ tree from {tree_file}")
        with open(tree_file, "rb") as f:
            # nj_tree, existing_ids = pickle.load(f)
            guide_tree, existing_ids = pickle.load(f)
    else:
        logger.info("initializing guide tree")
        nj_tree = f"{args.save_dir}/output_half.nwk"
        existing_ids = set(f_sequences2.keys())
        guide_tree = parse_newick_to_guide(nj_tree)
        # 保存新计算的树
        with open(tree_file, "wb") as f:
            # pickle.dump((nj_tree, existing_ids), f)
            pickle.dump((guide_tree, existing_ids), f)
        logger.info(f"Saved NJ tree to {tree_file}")
    # guide_tree = convert_nj_to_guide(nj_tree)
    logger.info("guide tree built")
    # 内存优化， 删除dist_matrix， 过滤sequences和clusters
    del dist_matrix
    sequences = {k: sequences[k] for k in existing_ids}
    clusters = {k: clusters[k] for k in existing_ids}
    # 删除f_dist_matrix2, f_sequences2
    del f_dist_matrix2, f_sequences2
    # 3. 初始化比对列表
    logger.info("initializing alignments")
    alignments = {
        "alignments": {k: {k: v} for k, v in sequences.items() if k in existing_ids},
        "cluster_alignments": {
            k: {k: v} for k, v in clusters.items() if k in existing_ids
        },
        "operations": {k: [] for k in existing_ids},
    }
    # 4. 按树结构合并比对
    logger.info("starting progressive alignment")
    # # 分层并行合并
    # level_dict = defaultdict(list)
    # current_level = 0
    # total_alignments = len(guide_tree)
    # processed_trees_num = 0
    # # 构建层级关系（假设guide_tree是后序遍历生成的）
    # temp_tree = guide_tree[::-1]  # 反转列表得到前序遍历顺序
    # for pair in temp_tree:
    #     if not any(p in level_dict[current_level - 1] for p in pair):
    #         current_level += 1
    #     level_dict[current_level].append(pair)

    # # 按层级并行处理
    # with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
    #     for level in sorted(level_dict.keys()):
    #         # 提交当前层级所有任务
    #         futures = [
    #             executor.submit(process_pair, pair) for pair in level_dict[level]
    #         ]

    #         # 等待当前层级完成
    #         for future in concurrent.futures.as_completed(futures):
    #             merged, (g1, g2) = future.result()
    #             new_key = (g1, g2)

    #             # 更新比对词典
    #             alignments["alignments"][new_key] = merged["alignments"][new_key]
    #             alignments["cluster_alignments"][new_key] = merged[
    #                 "cluster_alignments"
    #             ][new_key]
    #             alignments["operations"][new_key] = merged["operations"][new_key]

    #             # 清理旧数据
    #             for g in [g1, g2]:
    #                 if g in alignments["alignments"]:
    #                     del alignments["alignments"][g]
    #                 if g in alignments["cluster_alignments"]:
    #                     del alignments["cluster_alignments"][g]
    #                 if g in alignments["operations"]:
    #                     del alignments["operations"][g]
    #             processed_trees_num += 1
    #             if processed_trees_num % 1000 == 0:
    #                 logger.info(
    #                     f"Merging {processed_trees_num} / {total_alignments} alignments. \n{mem_monitor()}"
    #                 )
    # 原有的单线程渐进式比对
    total_alignments = len(guide_tree)
    processed_trees_num = 0
    # 修改检查点恢复逻辑
    if args.load_existing and os.path.exists(checkpoint_file):
        logger.info(f"Loading checkpoint from {checkpoint_file}")
        with open(checkpoint_file, "rb") as f:
            checkpoint_data = pickle.load(f)
        # 加载所有必要状态
        # dist_matrix = checkpoint_data['dist_matrix']
        # guide_tree = checkpoint_data['guide_tree']
        alignments = checkpoint_data["alignments"]
        processed_pairs = checkpoint_data["processed_pairs"]  # 已处理pair列表
        remaining_pairs = checkpoint_data["remaining_pairs"]  # 剩余待处理pair列表
        last_save_time = checkpoint_data.get("last_save_time", time.time())
        logger.info(
            f"Resuming from {len(processed_pairs)} processed pairs, {len(remaining_pairs)} remaining"
        )
    # 修改渐进比对循环结构（将guide_tree转换为可续传的迭代器）
    remaining_pairs = (
        list(guide_tree)
        if not (args.load_existing and os.path.exists(checkpoint_file))
        else remaining_pairs
    )
    processed_pairs = (
        []
        if not (args.load_existing and os.path.exists(checkpoint_file))
        else processed_pairs
    )
    # for pair in guide_tree:
    p_count = 0
    while remaining_pairs:
        pair = remaining_pairs.pop(0)
        p_count += 1
        # 从引导树中获取要合并的两个组ID
        # print(type(pair[0]))
        # print(type(pair[1]))
        # group1_id, group2_id = frozenset([pair[0]]), frozenset([pair[1]])
        group1_id, group2_id = pair[0], pair[1]
        # group1_id = (
        #     group1_id.replace("'", "") if isinstance(group1_id, str) else group1_id
        # )
        # group2_id = (
        #     group2_id.replace("'", "") if isinstance(group2_id, str) else group2_id
        # )
        # print(f"Merging {group1_id} and {group2_id}")
        # 每1000个比对打印一次日志
        if (
            processed_trees_num % 1000 == 0
            or processed_trees_num == total_alignments - 1
        ):
            logger.info(
                f"Merging {processed_trees_num} / {total_alignments} alignments. \n{mem_monitor()}"
            )
        # 查找包含当前组的比对 (modified with explicit list conversion and fallback)
        merged = merge_alignments(
            {
                "alignments": {group1_id: alignments["alignments"][group1_id]},
                "cluster_alignments": {
                    group1_id: alignments["cluster_alignments"][group1_id]
                },
            },
            {
                "alignments": {group2_id: alignments["alignments"][group2_id]},
                "cluster_alignments": {
                    group2_id: alignments["cluster_alignments"][group2_id]
                },
            },
            group1_id,
            group2_id,
        )
        # 更新比对词典
        new_key = (group1_id, group2_id)
        alignments["alignments"][new_key] = merged["alignments"][new_key]
        alignments["cluster_alignments"][new_key] = merged["cluster_alignments"][
            new_key
        ]
        alignments["operations"][new_key] = merged["operations"][new_key]
        del alignments["alignments"][group1_id]
        del alignments["alignments"][group2_id]
        del alignments["cluster_alignments"][group1_id]
        del alignments["cluster_alignments"][group2_id]
        del alignments["operations"][group1_id]
        del alignments["operations"][group2_id]
        processed_trees_num += 1
        processed_pairs.append(pair)

        # 时间间隔检查
        current_time = time.time()
        if current_time - last_save_time >= checkpoint_interval:
            checkpoint_data = {
                # 'dist_matrix': dist_matrix,
                # 'guide_tree': guide_tree,
                "alignments": alignments,
                "processed_pairs": processed_pairs,
                "remaining_pairs": remaining_pairs,
                "last_save_time": current_time,
            }
            with open(checkpoint_file, "wb") as f:
                pickle.dump(checkpoint_data, f)
            last_save_time = current_time
            logger.info(f"Time-based checkpoint saved at {time.ctime(current_time)}")
    # the last pair
    group1_id = list(alignments["alignments"].keys())[0]
    group2_id = list(alignments["alignments"].keys())[1]
    merged = merge_alignments(
        {
            "alignments": {group1_id: alignments["alignments"][group1_id]},
            "cluster_alignments": {
                group1_id: alignments["cluster_alignments"][group1_id]
            },
        },
        {
            "alignments": {group2_id: alignments["alignments"][group2_id]},
            "cluster_alignments": {
                group2_id: alignments["cluster_alignments"][group2_id]
            },
        },
        group1_id,
        group2_id,
    )
    # 更新比对词典
    new_key = (group1_id, group2_id)
    alignments["alignments"][new_key] = merged["alignments"][new_key]
    alignments["cluster_alignments"][new_key] = merged["cluster_alignments"][new_key]
    alignments["operations"][new_key] = merged["operations"][new_key]
    del alignments["alignments"][group1_id]
    del alignments["alignments"][group2_id]
    del alignments["cluster_alignments"][group1_id]
    del alignments["cluster_alignments"][group2_id]
    del alignments["operations"][group1_id]
    del alignments["operations"][group2_id]
    logger.info(
        f"Merging {processed_trees_num} / {total_alignments} alignments. \n{mem_monitor()}"
    )
    # 最终返回前添加迭代优化
    logger.info("starting iterative refinement")
    final_alignment = alignments
    # final_id = (
    #     list(final_alignment["alignments"].keys())[-1]
    #     if len(list(final_alignment["alignments"].keys())[-1]) > 1
    #     else list(final_alignment["alignments"].keys())[0]
    # )
    final_id = list(final_alignment["alignments"].keys())[0]
    if len(sequences) > 3:
        final_alignment = iterative_refinement(final_alignment, final_id=final_id)
    # 最终清理时
    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)
    return final_alignment


@logger.catch
def transform_sequence(fimo_res, seq_lengths):
    """
    将FIMO结果和序列长度信息转换为motif+间隔的结构
    参数:
        fimo_res: DataFrame包含过滤后的motif位置信息
        seq_lengths: DataFrame包含每个序列的总长度
    返回:
        seq_dict: 字典，键为序列名，值为motif和间隔交替组成的列表
    """
    # 初始化存储结构的字典
    seq_dict = dict()
    cluster_dict = dict()
    # 遍历每个序列的长度信息
    for idx, row in seq_lengths.iterrows():
        seq_name = row["Sequence Name"]
        seq = []  # 存储当前序列的motif+间隔结构
        clusters = []
        seq_length = row["Sequence Length"]
        # 获取当前序列的所有motif并排序
        seq_motifs = fimo_res[fimo_res["sequence_name"] == seq_name]
        seq_motifs = seq_motifs.sort_values(
            by=["start", "stop"], ascending=[True, True]
        ).reset_index()
        last_stop = 0  # 跟踪上一个motif的结束位置
        # 遍历当前序列的每个motif
        for idx2, row2 in seq_motifs.iterrows():
            # 计算并添加前导间隔
            seq.append(row2["start"] - last_stop - 1)
            # 添加当前motif ID
            seq.append(row2["motif_id"])
            clusters.append(row2["start"] - last_stop - 1)
            clusters.append(row2["motif_cluster"])
            # 如果是最后一个motif，添加尾部间隔
            if idx2 == seq_motifs.shape[0] - 1:
                seq.append(seq_length - row2["stop"])
                clusters.append(seq_length - row2["stop"])
            last_stop = row2["stop"]  # 更新最后结束位置
        # 将当前序列的结构存入字典
        seq_dict = seq_dict | {seq_name: seq}
        # 将当前序列的结构反向后存入字典
        seq_dict = seq_dict | {seq_name + "_rev": seq[::-1]}
        # # 将当前序列的cluster结构存入词典
        cluster_dict = cluster_dict | {seq_name: clusters}
        # # 将当前序列的cluster结构反向后存入词典
        cluster_dict = cluster_dict | {seq_name + "_rev": clusters[::-1]}
    return seq_dict, cluster_dict


# define a function to save alignment result
def save_alignment_result(alignment_result, output_file):
    with open(output_file, "w") as f:
        for key, value in alignment_result.items():
            f.write(f">{key}\n")
            f.write(f"{value}\n")
    logger.info(f"save alignment result to {output_file}")
    return output_file


def save_alignment_operation(alignment_result, output_file):
    with open(output_file, "w") as f:
        for key, value in alignment_result.items():
            f.write(f">{key}\n")
            f.write(f"{value}\n")
    logger.info(f"save alignment result to {output_file}")
    return output_file


if __name__ == "__main__":
    # meme_cluster = pd.read_csv(
    #     "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/JASPAR_denovo_motifs_clusters.tab",
    #     sep="\t",
    #     header=0,
    # )
    # fimo_res = pd.read_csv(
    #     "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/fimo_filtered_res.tsv",
    #     sep="\t",
    #     header=0,
    # )
    # fa_file = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/top_cres.fa"
    # seq_lengths = get_fasta_sequences(fa_file)
    # seq_lengths.to_csv("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/seq_lengths.csv", index=False)
    # seq_lengths = pd.read_csv(
    #     "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/seq_lengths.csv",
    #     sep=",",
    #     header=0,
    # )
    # seq_dict, cluster_dict = transform_sequence(fimo_res, seq_lengths)
    # # use pickle to save seq_dict and cluster_dict
    # with open(
    #     "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/seq_dict.pkl",
    #     "wb",
    # ) as f:
    #     pickle.dump(seq_dict, f)
    # with open(
    #     "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/cluster_dict.pkl",
    #     "wb",
    # ) as f:
    #     pickle.dump(cluster_dict, f)
    # 读取motif和cluster信息
    with open(
        "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/seq_dict.pkl",
        "rb",
    ) as f:
        seq_dict = pickle.load(f)
    with open(
        "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/cluster_dict.pkl",
        "rb",
    ) as f:
        cluster_dict = pickle.load(f)
    seq_dict = {k: v for k, v in seq_dict.items() if v != []}
    cluster_dict = {k: v for k, v in cluster_dict.items() if v != []}
    if args.group:
        # 读取group，根据group筛选seq_dict和cluster_dict
        group = pd.read_csv(args.group, sep="\t", header=None, names=["Sequence_Name"])
        group = group["Sequence_Name"].tolist()
        # 增加‘_rev’后缀
        group_rev = [i + "_rev" for i in group]
        group = group + group_rev
        seq_dict = {k: v for k, v in seq_dict.items() if k in group}
        cluster_dict = {k: v for k, v in cluster_dict.items() if k in group}
    all_msa_result = progressive_motif_alignment(seq_dict, cluster_dict)
    save_alignment_result(
        list(all_msa_result["alignments"].values())[0],
        f"{args.outdir}/all_msa_motif_alignment.txt",  # /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/para_conti/all_msa_motif_alignment.txt
    )
    save_alignment_result(
        list(all_msa_result["cluster_alignments"].values())[0],
        f"{args.outdir}/all_msa_cluster_alignment.txt",  # /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/para_conti/all_msa_cluster_alignment.txt
    )
    save_alignment_operation(
        all_msa_result["operations"],
        f"{args.outdir}/all_msa_operation.txt",  # /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/para_conti/all_msa_operation.txt
    )
    # use pickle to save all_msa_result
    with open(
        f"{args.outdir}/all_msa_result.pkl",  # /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/para_conti
        "wb",
    ) as f:
        pickle.dump(all_msa_result, f)


# dict to pd.dataframe
def dict_to_dataframe(data_dict):
    """将字典转换为DataFrame"""
    df = pd.DataFrame.from_dict(data_dict, orient="index")
    df.reset_index(inplace=True)
    df.columns = ["Sequence Name"] + [
        f"Position {i}" for i in range(1, len(df.columns))
    ]
    return df
