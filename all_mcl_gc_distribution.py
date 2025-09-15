#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Author: zhangz
@Date: 2025-05-19
@Description: 计算参与mcl聚类的所有序列的GC含量分布
"""
import os
import sys
from Bio import SeqIO
import csv
from pathlib import Path
from loguru import logger
import concurrent.futures
import numpy as np
import pandas as pd

logger.add(
    sys.stdout,
    level="INFO",
    format="{time} {level} {message}",
    backtrace=True,
    diagnose=True,
)


@logger.catch
def calculate_gc(seq):
    """
    计算序列的GC含量
    """
    seq = seq.upper()
    gc_content = (
        ((seq.count("G") + seq.count("C")) / len(seq)) * 100 if len(seq) > 0 else 0
    )
    return gc_content


@logger.catch
def calculate_gc_distribution(fa_file):
    """
    并行计算用于mcl聚类的所有序列的GC含量，形成一个分布，用于判断产生的随机序列是否在这个分布的95%区间内
    """
    gc_contents = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=90) as executor:
        futures = []
        for record in SeqIO.parse(fa_file, "fasta"):
            future = executor.submit(calculate_gc, str(record.seq))
            futures.append(future)

        for future in concurrent.futures.as_completed(futures):
            gc_contents.append(future.result())
    return gc_contents


@logger.catch
def main():
    """
    主函数
    """
    # 输入文件路径
    fa_file = (
        "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/all_cre.fa"
    )

    # 计算GC含量分布
    gc_contents = calculate_gc_distribution(fa_file)

    # save gc_contents
    gc_contents_df = pd.DataFrame(gc_contents, columns=["GC_Content"])
    gc_contents_df.to_csv(
        "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/gc_contents.csv",
        index=False,
    )
    # gc_contents_df = pd.read_csv(
    #     "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/gc_contents.csv"
    # )
    # 计算95%分位数
    gc_975th_percentile = np.percentile(gc_contents, 97.5)

    # 输出结果
    logger.info(f"GC含量97.5%分位数: {gc_975th_percentile}")

    # 计算2.5%分位数
    gc_025th_percentile = np.percentile(gc_contents, 2.5)
    logger.info(f"GC含量2.5%分位数: {gc_025th_percentile}")


if __name__ == "__main__":
    main()
