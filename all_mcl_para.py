#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from multiprocessing import Manager
from loguru import logger
import subprocess
import sys
import concurrent.futures
from Bio import SeqIO
from functools import partial
import os
import itertools
import re

# from blast_all_mcl import make_similarity_matrix

## the all-against-all blast task has been completed, but the mcl clustering has not been done yet
## so we need to do the mcl clustering for each species
## this script only do mcl clustering

species = [
    "Ovis_aries",
    "Bos_taurus",
    "Neophocaena_asiaeorientalis",
    "Sus_scrofa",
    "Lama_glama",
    "Mustela_putorius",
    "Canis_lupus",
    "Felis_catus",
    "Equus_asinus",
    "Equus_caballus",
    "Rhinolophus_pusillus",
    "Myotis_chinensis",
    "Atelerix_albiventris",
    "Mus_musculus",
    "Rattus_norvegicus",
    "Cavia_porcellus",
    "Oryctolagus_cuniculus",
    "Macaca_mulatta",
    "Rhinopithecus_roxellana",
    "Tupaia_belangeri",
    "Procavia_capensis",
    "Petaurus_breviceps",
]


@logger.catch
def process_blast_chunk(chunk_file, mcl_input_file, cre_ids, processed_pairs):
    cre_id_set = set(
        cre_ids["query_id"].values
    )  # 将 cre_ids["query_id"] 转换为集合以加速查找
    lines_to_write = []

    with open(chunk_file, "r") as f:
        logger.info(f"processing chunk {chunk_file}")
        for line in itertools.islice(f, 0, None):
            line = line.strip().split("\t")
            query = line[0]
            query_name = "_".join(query.split("_")[0:5])
            subject = line[1]
            subject_name = "_".join(subject.split("_")[0:5])
            similarity = float(line[2])  # 提前转换为 float 类型

            if (
                (query_name, subject_name) not in processed_pairs
                and similarity >= 50
                and query_name in cre_id_set
                and query_name != subject_name
            ):
                lines_to_write.append(f"{query}\t{subject}\t{similarity}\n")
                processed_pairs.append((query_name, subject_name))
                processed_pairs = list(set(processed_pairs))

    with open(mcl_input_file, "a") as out:
        out.writelines(lines_to_write)
    logger.info(f"{chunk_file} done, remove")

    # remove chunk file
    os.remove(chunk_file)


@logger.catch
def main(blast_out, mcl_input_file, cre_ids):
    # 分割 BLAST 结果文件
    # chunk_size = 100000  # 每个小文件包含的行数
    # chunk_files = []
    # with open(blast_out, "r") as f:
    #     chunk_file = None
    #     chunk = []
    #     # remove former chunk files
    #     subprocess.run(
    #         "find /media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/ -name 'all_cre_blast.out_chunk_*.txt' -exec rm -f {} +",
    #         shell=True,
    #     )
    #     for i, line in enumerate(f):
    #         if i % chunk_size == 0:
    #             logger.info(f"processing chunk {i//chunk_size}")
    #             if chunk_file:
    #                 chunk_file.writelines(chunk)
    #                 chunk_file.close()
    #             chunk_filename = f"/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/chunk_{i//chunk_size}.txt"
    #             chunk_files.append(chunk_filename)
    #             chunk_file = open(chunk_filename, "w")
    #             chunk = []
    #         chunk.append(line)
    #     if chunk_file:
    #         chunk_file.writelines(chunk)
    #         chunk_file.close()

    chunk_dir = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/"
    chunk_files = [
        os.path.join(chunk_dir, f)
        for f in os.listdir(chunk_dir)
        if f.startswith("chunk_")
    ]

    # # 注册自定义集合类型
    # @logger.catch
    # def Manager():
    #     m = mp.Manager()
    #     m.register("custom_set", set)
    #     return m

    # manager = Manager()
    # @logger.catch
    # def create_set():
    #     return set()

    # manager = Manager()
    # manager.register("custom_set", create_set)
    # manager.register("custom_set", create_set)
    # processed_pairs = manager.custom_set()
    # manager = mp.Manager()
    # # 注册自定义类型
    # manager.register("custom_set", create_set)
    # # 创建共享的 set 对象
    # processed_pairs = manager.custom_set()
    # processed_pairs = manager.set()
    # AttributeError: 'SyncManager' object has no attribute 'set'

    # 使用 Manager 来管理进程安全的列表
    manager = Manager()
    processed_pairs = manager.list()
    logger.info("start transform blast results")

    # 使用多线程处理每个小文件
    with concurrent.futures.ProcessPoolExecutor(max_workers=50) as executor:
        futures = []
        for chunk_file in chunk_files:
            futures.append(
                executor.submit(
                    process_blast_chunk,
                    chunk_file,
                    mcl_input_file,
                    cre_ids,
                    processed_pairs,
                )
            )
        concurrent.futures.wait(futures)

    # 保存 processed_pairs
    processed_file = blast_out.replace(".out", "_processed_pairs.txt")
    with open(processed_file, "w") as f:
        for pair in processed_pairs:
            f.write(f"{pair[0]}\t{pair[1]}\n")
    logger.info(f"processed_pairs saved as {processed_file}")

    # # 合并所有小文件的处理结果
    # with open(mcl_input_file, "w") as out:
    #     for chunk_file in chunk_files:
    #         with open(chunk_file, "r") as chunk:
    #             out.write(chunk.read())
    #         os.remove(chunk_file)  # 删除临时小文件

    logger.info(f"make similarity matrix for {blast_out} done")
    mci_file = mcl_input_file.replace(".txt", ".mci")
    tab_file = mcl_input_file.replace(".txt", ".tab")
    mcl_file = mcl_input_file.replace("_input.txt", ".mcl")
    logger.info("start mcl input convert by mxcload")
    subprocess.run(
        f"mcxload -abc {mcl_input_file} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {mci_file} -write-tab {tab_file}",
        shell=True,
    )
    logger.info("start mcl clustering")
    subprocess.run(
        f"mcl {mci_file} -I 2 -use-tab {tab_file} -o {mcl_file} -V pruning -te 60",
        shell=True,
    )


@logger.catch
def make_cre_ids(species):
    def process_file(file_path, sp):
        for chunk in pd.read_csv(file_path, sep=",", header=0, chunksize=1000):
            chunk["unique_id"] = chunk.apply(
                lambda row: f"{row['unique_id']}_{row['peak']}", axis=1
            )
            chunk["query_id"] = chunk.apply(lambda row: f"{sp}_{row['id']}", axis=1)
            filtered_chunk = chunk[chunk["length"] < 500][["unique_id", "query_id"]]
            yield filtered_chunk

    cre_ids = pd.DataFrame()
    for sp in species:
        file_path = (
            f"/media/Data/zhangz/chip/analysis/summary2/anno_sum/{sp}/{sp}_all.csv"
        )
        for filtered_chunk in process_file(file_path, sp):
            cre_ids = pd.concat([cre_ids, filtered_chunk], ignore_index=True)

    return cre_ids


if __name__ == "__main__":
    blast_out = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/all_cre_blast.out"
    mcl_input_file = "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/chunk/mcl/all_cre_blast_para_mcl_input.txt"
    # cre_ids = make_cre_ids(species)
    # cre_ids.to_csv("/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/all_cre_ids.txt", index=False)
    cre_ids = pd.read_csv(
        "/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/all_cre_ids.txt",
        sep=",",
        header=0,
    )
    main(blast_out, mcl_input_file, cre_ids)
