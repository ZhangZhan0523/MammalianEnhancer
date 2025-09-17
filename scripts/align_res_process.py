#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
time: 2025/04/25
author: Zhang Zhan
we have aligned motif-based cre sequences from the most prevelent cre cluster,
now we need to process the alignment result,
1. count apearence of motifs of spacer lengths in each column,
2. draw distribution of each postion,
3. find the presentative motif cluster in each column,
4. build a consensus pattern of the whole alignment
5. randomly produce a sequence according to the alignment and distribution
6. blast the sequence against the database to keep enough de novo sequences
"""

import os
import sys
import pickle
import re
import subprocess
import numpy as np
import pandas as pd
from loguru import logger
import random
import concurrent.futures
import tempfile
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse

logger.add(
    sys.stdout,
    format="{time} {level} {message}",
    filter="my_module",
    level="INFO",
    backtrace=True,
    diagnose=True,
)

MOTIF_DB = "/media/Data/zhangz/chip/genomes/motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
MOTIF_THRESHOLD = 3


# functions
# dict to pd.dataframe
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


# a function to count the number of each motif in each column
@logger.catch
def count_motif(df):
    """统计每个列中的motif数量"""
    motif_counts = {}
    for col in df.columns:
        if col.startswith("motif_"):
            motif_counts[col] = df[col].value_counts().to_dict()
    return motif_counts


# a function to draw distribution of each postion
@logger.catch
def draw_distribution(df):
    """绘制每个位置的分布"""
    for col in df.columns:
        if col.startswith("motif_"):
            df[col].value_counts().plot(kind="bar")


# 计算cluster比对的每一列的分数，记录每一列的行数，去除’-‘后的有效值的数量，
# 去除’-‘后的unique数量，以及频率最高的unique值的频率，
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
def calculate_spacer_distribution(df):
    """计算spacer长度分布"""
    spacer_lengths = df.filter(like="spacer_").applymap(len)
    return spacer_lengths.describe()


@logger.catch
def parse_meme_file(meme_path):
    """
    解析MEME文件，返回基序信息（名称、长度、PWM矩阵）
    :param meme_path: MEME文件路径
    :return: 字典列表，每个字典包含 'name'（基序名）、'width'（长度）、'pwm'（PWM矩阵）
    """
    motifs = {}
    current_motif = None

    with open(meme_path, "r") as f:
        for line in f:
            line = line.strip()

            # 匹配MOTIF起始行（如 "MOTIF M001_1.0001"）
            if line.startswith("MOTIF"):
                if current_motif:
                    motifs[motif_name] = current_motif  # 保存上一个基序
                motif_name = (
                    line.split()[2]
                    if os.path.basename(meme_path) == "combined.meme"
                    else line.split()[1]
                )
                if motif_name.startswith("MA"):
                    motif_name = motif_name.split("-")[0]
                if "MEME" in motif_name or "STREME" in motif_name:
                    motif_name = "-".join(motif_name.split("-")[:-2])
                current_motif = {"name": motif_name, "pwm": []}

            # 匹配PWM矩阵行（如 "letter-probability matrix: alength= 4 w= 6 ..."）
            elif (
                current_motif
                and "pwm" in current_motif
                and line.startswith("letter-probability matrix")
            ):
                # 提取基序长度（w=后面的数值）
                width_match = re.search(r"w= (\d+)", line)
                if width_match:
                    current_motif["width"] = int(width_match.group(1))

            # 读取PWM矩阵数据（每行对应一个位置的A/C/G/T概率）
            elif (
                current_motif
                and "width" in current_motif
                and len(current_motif["pwm"]) < current_motif["width"]
            ):
                # 跳过非数值行（如注释）
                if all(token.replace(".", "").isdigit() for token in line.split()):
                    probabilities = list(map(float, line.split()))
                    current_motif["pwm"].append(probabilities)

        # 保存最后一个基序
        if current_motif:
            motifs[motif_name] = current_motif

    return motifs


@logger.catch
def generate_motif_sequences(motif, num_sequences=1):
    """
    根据基序的PWM生成随机序列
    :param motif: 解析后的基序字典（包含'pwm'和'width'）
    :param num_sequences: 生成序列数量
    :return: 生成的序列列表
    """
    nucleotides = ["A", "C", "G", "T"]  # 对应PWM矩阵的列顺序（A/C/G/T）
    # sequences = []
    for _ in range(num_sequences):
        sequence = []
        for position_pwm in motif["pwm"]:
            # 按概率加权抽样（确保概率和为1，避免浮点误差）
            prob = np.array(position_pwm) / np.sum(position_pwm)
            chosen_nuc = np.random.choice(nucleotides, p=prob)
            sequence.append(chosen_nuc)
        # sequences.append(''.join(sequence))
    return "".join(sequence)


# 整理比对的motif_df和cluster_df，选择cluster_df中ValidCount大于1000（参数）的motif比对列，将每两个比对列中间的spacer列合并成同一个分布
# 去除'-'的行，调整存储结构，方便对每一列进行随机抽样
@logger.catch
def process_alignment_dict(motif_df, cluster_stat, threshold=1000):
    """选择cluster_df中ValidCount大于1000的motif比对列，将每两个比对列中间的spacer列合并成同一个分布，去除'-'的行，调整存储结构，用字典存储"""
    # 过滤有效列
    # 获取有效motif列 (基于cluster统计)
    valid_cols = cluster_stat[cluster_stat["ValidCount"] > threshold].index.tolist()

    # 如果列名以motif_开头，且不在有效列中，则删除，以spacer_开头的列则保留
    motif_cols = [
        col
        for col in motif_df.columns
        if col.startswith("motif_") and col in valid_cols
    ]
    spacer_cols = [col for col in motif_df.columns if col.startswith("spacer_")]
    filtered_cols = motif_cols + spacer_cols
    filtered_cols = sorted(filtered_cols, key=lambda x: int(x.split("_")[1]))
    filtered_df = motif_df[filtered_cols].copy()

    # 遍历每一列，去除'-'的行，如果这一列和下一列都是spacer列，则合并成同一个分布，把同一列的数据转换成list，把列名作为key，list作为value
    processed_data = {}
    last_col = "spacer_0"
    processed_data[last_col] = motif_df[last_col].replace("-", np.nan).dropna().tolist()
    for col in filtered_df.columns:
        if col == last_col:
            continue
        if col.startswith("spacer_"):
            if last_col.startswith("spacer_"):
                # 如果前一列也是spacer列，则合并
                processed_data[last_col] += (
                    filtered_df[col].replace("-", np.nan).dropna().tolist()
                )
            else:
                filtered_series = filtered_df[col].replace("-", np.nan).dropna()
                processed_data[col] = filtered_series.tolist()
                last_col = col
        elif col.startswith("motif_"):
            # 如果是motif列，去除'-'的行
            filtered_series = filtered_df[col].replace("-", np.nan).dropna()
            processed_data[col] = filtered_series.tolist()
            last_col = col
    return processed_data


@logger.catch
def merge_meme_files(meme_file1, meme_file2):
    """合并两个MEME文件中的motif"""
    """meme_file1是原始jaspar文件，meme_file2是经过处理的combined.meme文件"""
    motifs1 = parse_meme_file(meme_file1)
    motifs2 = parse_meme_file(meme_file2)
    # motif1 和motif2都是字典，将motif2中不存在于motif1中的motif添加到motif1中
    for motif_name, motif in motifs2.items():
        if motif_name not in motifs1:
            motifs1[motif_name] = motif
    return motifs1


@logger.catch
def run_fimo(input_fasta, motif_db, seq_seed):
    """执行fimo分析并返回结果路径"""
    output_dir = f"/media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/exp_seq/spacer_tmp/{seq_seed}/fimo_out"
    denovo_out_dir = f"/media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/exp_seq/spacer_tmp/{seq_seed}/denovo_fimo_out"
    # 创建输出目录
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    os.makedirs(os.path.dirname(denovo_out_dir), exist_ok=True)
    fimo_denovo_cmd = f"conda run -n meme fimo --verbosity 1 --oc {denovo_out_dir} --motif 2-GGCAAGAATACTGGA \
        --motif GGAGTGGGTTGCCAT --motif 1-GGAAGATCCCCTGGA --motif 4-ACCCAGGGATCGAAC --motif GGGATTYTCCAGGCA \
        --motif GACAGAGGAGCCTGG  --motif 14-GAAGGGAWWGGCWAC --motif AGGCAGATTCTTTAC  --motif 6-GCTACAGTCCATGGG \
        --motif 7-GTCTCCTGCATTGCA --motif 9-AGAGTCRGACACGAC --motif CTGAGCCACCAGGGA --motif 23-GGCTTCCCWGGTG \
        --motif 36-CCRCCAGGCTC --motif 11-TGAGCGACTAAMMHW --motif TCTTTGYGACC --motif TYRGGAAG \
        --motif 20-TGAAGYGACTTAGCA --motif 37-CAGGAGATGYRG --motif 62-GCAACTAAGMC --motif YWWCACTTTCACTTT \
        --motif 61-CTGCATTGGCA --motif 15-GTTGCAA --motif 64-AGAATCCCCATGGA --motif 16-ACACACACACACACA \
        --motif 88-CAACTACTGAGCC --motif 28-AACAACAACAAM --motif 89-AAGTGAAGTG --motif 18-TAAAAAAAAAAAAAA \
        --motif 79-AYTGAGCACYTACT --motif 91-AGTCCTAACCACTG --motif 21-TATTTTWAAAATA --motif 85-AGCACAGCACAGC \
        --motif 73-TCATTCATTCATTCA --motif 22-AAGATCCCACATGCY --motif 24-CCCCACCCCCACCM --motif 26-TAAATAAATAW \
        --motif 63-ACTCATTGGAAAAGA --motif 42-WTTTTAAAAW --motif 48-GCAGCAGCAGCAGCA --motif 72-AGTTCCCTGACCAGG \
        --motif 29-CCTCAGTTTCCTCAT --motif 59-AAACTAACACAACAT --motif 84-GATCTTAGTTCCCY --motif 38-ATGGACATGAGTTTG \
        --motif 41-ACGCCAGGCCTCCCT --motif 60-AAGACCCAGCACAGC --motif 107-ACTGAGCTATGA --motif 97-HGCATGCATGCD \
        --motif 96-TWAAATGARATA --motif 45-CCAAGGTCACACA --motif 87-GAACTGAACTGAA --motif 44-ATGAGATGGTTGGAT \
        --motif 39-CATCACCRACTCM --motif 51-CCCATTTTACAG --motif 68-AGGCTCAGAGAGGT --motif 81-CGGCCGCCGGCCC \
        --motif 76-GTGCCAGGCACTGT --motif 75-CSCGGCCCGGCCSG --motif 86-AAGGACTGATGCT --motif 52-TCTTTCCCAGCATCA \
        --motif 70-AGAGCACAGGCTC --motif 67-CAACTAGAGAAARSC --motif 71-ATGATGATGATG --motif 83-GACTCTCAAGAGTC \
        --motif 106-GTAAAKCAACTAT --motif 109-CCTAAAGGAAATCA --motif 55-CAGGAGGAGAAGGGG --motif 100-CAATGAATATTCAGG \
        --motif 90-TAGCACAGGGAACT --motif 98-GCAAGGAGATCMAA \
        /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/xstreme_out/streme_out/streme.xml \
        {input_fasta}"
    subprocess.run(fimo_denovo_cmd, shell=True, check=True)
    subprocess.run(
        f"conda run -n meme fimo --verbosity 1 --oc {output_dir} {motif_db} {input_fasta}",
        check=True,
        shell=True,
    )
    # combine fimo results
    try:
        denovo_fimo_out = pd.read_csv(
            denovo_out_dir + "/fimo.tsv", sep="\t", comment="#"
        )
    except (FileNotFoundError, pd.errors.EmptyDataError):
        denovo_fimo_out = pd.DataFrame()
    try:
        fimo_out = pd.read_csv(output_dir + "/fimo.tsv", sep="\t", comment="#")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        fimo_out = pd.DataFrame()
        # 合并fimo结果
    combined_fimo_out = pd.concat([denovo_fimo_out, fimo_out], ignore_index=True)
    combined_fimo_out = (
        combined_fimo_out[combined_fimo_out["q-value"] < 0.05]
        if not combined_fimo_out.empty
        else combined_fimo_out
    )
    return combined_fimo_out


@logger.catch
def parse_fimo_results(fimo_output):
    """解析fimo结果并返回motif计数"""
    try:
        df = pd.read_csv(fimo_output, sep="\t", comment="#")
        return len(df)
    except FileNotFoundError:
        return 0


@logger.catch
def generate_random_spacer(length, seq_seed):
    """生成随机spacer序列"""
    import random

    # while True:
    # 生成随机序列
    seq = "".join(random.choices(["A", "T", "C", "G"], k=length))
    # 保存为fasta格式，临时文件，调整命名避免重复
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_fasta:
        tmp_fasta.write(f">random_spacer\n{seq}\n")
    # 运行fimo，并且解析结果，如果motif数量小于阈值，返回seq
    # logger.info(f"Running fimo for sequence: {seq_seed}")
    fimo_res = run_fimo(tmp_fasta.name, MOTIF_DB, seq_seed)
    motif_count = fimo_res.shape[0]
    logger.info(f"Motif count for sequence {seq_seed}: {motif_count}")
    # if motif_count < MOTIF_THRESHOLD:
    return seq, motif_count
    os.remove(tmp_fasta.name)
    os.remove(fimo_res)


@logger.catch
def blast_denovo_seq(seq, seq_seed, seq_name="group2"):
    """和已知序列去比较，blast或搜索，没有重合（identity小于90）就保存到fasta文件中"""
    blast_db = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/all_cre"
    # save seq as a temp fasta file, adjust name to avoid conflict
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_fasta:
        tmp_fasta.write(f">random_seq\n{seq}\n")
        new_fa = tmp_fasta.name
    blast_out = f"/media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/exp_seq/blast_tmp/{seq_name}_{seq_seed}_blast_out"
    logger.info(f"Running blast for sequence: {seq_name}_{seq_seed}")
    subprocess.run(
        f"conda run -n R blastn -query {new_fa} -db {blast_db} -evalue 1e-5 -strand both -perc_identity 99 -task dc-megablast -outfmt 6 -out {blast_out} -num_threads 30",
        shell=True,
        check=True,
    )
    # return line number of blast_out counted by shell command
    hit_num = subprocess.run(
        f"wc -l {blast_out}",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    return int(hit_num.stdout.strip().split()[0])
    os.remove(new_fa)
    os.remove(blast_out)


@logger.catch
def generate_valid_sequence(
    processed_data, combined_motifs, seq_seed=1, seq_name="group2"
):
    # TODO: 这个作为主流程
    # TODO: 写一个generate_motif_seq和generate_spacer_seq的函数，分别生成motif和spacer的序列，然后拼接成一个序列，然后和已知序列去比较，blast或搜索，没有重合就保存到fasta文件中
    """根据处理过的数据生成符合要求的序列，使用generate_motif_sequences和generate_random_spacer分别生成motif和spacer，连接在一起，保存到fasta文件中"""

    # output_fasta = f"/media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/exp_seq/{seq_name}_seq.fasta"
    # 生成motif和spacer的序列，然后拼接成一个序列，然后和已知序列去比较，blast或搜索，没有重合就保存到fasta文件中
    blast_hit = 1
    while True:
        seq_temp = []
        for seq_type, seq_dist in processed_data.items():
            if seq_type.startswith("motif_"):
                # 生成motif序列
                motif = random.choices(seq_dist)[0]
                seq = generate_motif_sequences(combined_motifs[motif])
                logger.info(
                    f"generated sequence for {seq_type}: {motif}, the seq is {seq}"
                )
            elif seq_type.startswith("spacer_"):
                # 生成spacer序列
                spacer_length = 100
                while spacer_length > 50:
                    spacer_length = random.choices(seq_dist)[0]
                spacer_motif_count = 1000
                while spacer_motif_count > MOTIF_THRESHOLD:
                    # 并行生成spacer序列，保留spacer如果spacer_motif_count小于MOTIF_THRESHOLD
                    with concurrent.futures.ProcessPoolExecutor(
                        max_workers=10
                    ) as executor:
                        futures = {
                            executor.submit(
                                generate_random_spacer,
                                int(spacer_length),
                                seq_seed2,
                            ): seq_seed2
                            for seq_seed2 in range(10)
                        }
                        for future in concurrent.futures.as_completed(futures):
                            seq_seed3 = futures[future]
                            try:
                                seq, spacer_motif_count = future.result()
                                if spacer_motif_count < MOTIF_THRESHOLD:
                                    break
                            except Exception as e:
                                logger.error(
                                    f"Error generating spacer sequence {seq_seed}: {e}"
                                )
                # seq = generate_random_spacer(int(spacer_length), seq_seed)
                logger.info(f"generated sequence for {seq_type}: {spacer_length}")
            seq_temp.append(seq)
        # save to output_fasta
        seq = "".join(seq_temp)
        if len(seq) > 500:
            continue
        blast_hit = blast_denovo_seq(seq, seq_seed)
        if blast_hit == 0:
            break
    record = SeqRecord(Seq(seq), id=f"{seq_name}_seq_{seq_seed}", description="")
    return record


@logger.catch
def generate_exp_seq(processed_data, seq_num=15, seq_name="group2"):
    """并行生成序列，保存到fasta文件中"""

    combined_meme = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/xstreme_out/combined.meme"
    combined_motifs = merge_meme_files(
        MOTIF_DB,
        combined_meme,
    )
    with concurrent.futures.ProcessPoolExecutor(max_workers=seq_num) as executor:
        futures = {
            executor.submit(
                generate_valid_sequence,
                processed_data,
                combined_motifs,
                seq_seed,
                seq_name,
            ): seq_seed
            for seq_seed in range(seq_num)
        }
        for future in concurrent.futures.as_completed(futures):
            seq_seed = futures[future]
            try:
                record = future.result()
                # 保存到fasta文件中
                output_fasta = f"/media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/exp_seq/{seq_name}_seq_{seq_seed}.fasta"
                SeqIO.write(record, output_fasta, "fasta")
                logger.info(f"Generated sequence {seq_seed} saved to {output_fasta}")
            except Exception as e:
                logger.error(f"Error generating sequence {seq_seed}: {e}")


if __name__ == "__main__":
    # read checkpoint file of alignment as test input
    # test_checkpoint = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/cache_checkpoint/checkpoint_eada1db5.pkl"
    # 从命令行读取name参数
    parser = argparse.ArgumentParser(
        description="Draw sequences from alignment results."
    )
    parser.add_argument(
        "--name", type=str, help="The name of the file to process", required=True
    )
    args = parser.parse_args()
    name = args.name
    test_checkpoint = f"/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/top/{name}/all_msa_result.pkl"
    # test_checkpoint = name
    with open(test_checkpoint, "rb") as f:
        test_checkpoint = pickle.load(f)
    # alignment = test_checkpoint["alignments"]
    alignment = test_checkpoint
    alignment_df = dict_to_dataframe(
        alignment["alignments"][list(alignment["alignments"].keys())[-1]]
    )
    cluster_df = dict_to_dataframe(
        alignment["cluster_alignments"][
            list(alignment["cluster_alignments"].keys())[-1]
        ]
    )
    cluster_stat = calculate_cluster_scores(cluster_df)
    (cluster_stat["ValidCount"] > 2218).sum()
    # 3
    (cluster_stat["TopFrequency"] > 0.5).sum()
    # 20
    # cluster_df中TopFrequency最大的列
    cluster_df[cluster_stat["TopFrequency"].idxmax()]
    cluster_stat.loc[cluster_stat[cluster_stat["TopFrequency"] > 0.5].index.tolist()]
    # cluster_df中TopFrequency大于0.5的列
    cluster_df[cluster_stat[cluster_stat["TopFrequency"] > 0.5].index.tolist()]
    alignment_df[cluster_stat[cluster_stat["TopFrequency"] > 0.5].index.tolist()]
    # cluster_df中Valid_count最大的列
    cluster_df[cluster_stat["ValidCount"].idxmax()]
    # cluster_df中UniqueCount最大的列
    cluster_df[cluster_stat["UniqueCount"].idxmax()]
    cluster_df[cluster_stat[cluster_stat["ValidCount"] > 1000].index.tolist()]
    # 调整存储结构，方便抽取序列
    validCountThreshold = (
        cluster_stat.sort_values(by="ValidCount", ascending=False).iloc[15][
            "ValidCount"
        ]
        + 1
    )
    valid_stat = cluster_stat[cluster_stat["ValidCount"] > validCountThreshold]
    with open(
        f"/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/top/{name}/valid_stat.pkl",
        "wb",
    ) as f:
        pickle.dump(valid_stat, f)
    processed_data = process_alignment_dict(
        alignment_df, cluster_stat, threshold=validCountThreshold
    )
    # save processed_data
    with open(
        f"/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/top/{name}/processed_data.pkl",
        "wb",
    ) as f:
        pickle.dump(processed_data, f)
    generate_exp_seq(processed_data, seq_num=15, seq_name=name)
    subprocess.run(
        "for dd in `seq 0 14`;do cat /media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/exp_seq/group1_seq_${dd}.fasta >> /media/Data/zhangz/chip/analysis/summary2/blast_all/dual_luci/exp_seq/group1_seq.fasta;done",
        shell=True,
        check=True,
    )


def read_res():
    name = "group1"
    with open(
        f"/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/top/{name}/valid_stat.pkl",
        "rb",
    ) as f:
        valid_stat = pickle.load(f)
    with open(
        f"/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/top/{name}/processed_data.pkl",
        "rb",
    ) as f:
        processed_data = pickle.load(f)
    motif_cluster = pd.read_csv(
        "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/meme/JASPAR_denovo_motifs_clusters.tab",
        sep="\t",
        header=0,
    )
    from collections import Counter

    motif_list = processed_data["motif_53"]
    motif_counter = Counter(motif_list)
    motif_freq_df = pd.DataFrame.from_dict(
        motif_counter, orient="index", columns=["count"]
    )
    motif_freq_df["motif_id_base"] = motif_freq_df.index
    motif_freq_df["motif_cluster"] = motif_freq_df["motif_id_base"].apply(
        lambda x: (
            motif_cluster[motif_cluster["id"].str.contains(x.split(".")[0])][
                "cluster"
            ].values[0]
            if not motif_cluster[
                motif_cluster["id"].str.contains(x.split(".")[0])
            ].empty
            else None
        )
    )

    set(processed_data[valid_stat["ValidCount"].idxmax()])
    name = "group2"
    with open(
        f"/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/top/{name}/valid_stat.pkl",
        "rb",
    ) as f:
        valid_stat = pickle.load(f)
