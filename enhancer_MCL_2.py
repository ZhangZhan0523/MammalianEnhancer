import pandas as pd
import numpy as np
import multiprocessing as mp
from loguru import logger
import subprocess
import sys


@logger.catch
def make_similarity_matrix(sp):
    # 为物种sp创建相似性矩阵
    # 读取BLAST结果
    blast_res = f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}_blast.out"
    ## query id, subject id, identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    mcl_input_file = (
        f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}_mcl_input.txt"
    )
    cre_info = pd.read_csv(
        f"/media/Data/zhangz/chip/analysis/summary2/anno_sum/{sp}/{sp}_all.csv",
        sep=",",
        header=0,
    )
    cre_info["unique_id"] = cre_info.apply(
        lambda row: f"{row['unique_id']}_{row['peak']}", axis=1
    )
    ## filter cres by length and
    processed_pairs = set()

    try:
        with open(blast_res, "r") as f, open(mcl_input_file, "w") as out:
            for line in f:
                line = line.strip().split("\t")
                query = line[0]
                query_name = "_".join(query.split("_")[0:3])
                subject = line[1]
                subject_name = "_".join(subject.split("_")[0:3])
                similarity = line[2]
                if (query_name, subject_name) in processed_pairs or (
                    subject_name,
                    query_name,
                ) in processed_pairs:
                    continue
                if float(line[2]) < 50:
                    continue
                if (
                    float(cre_info["length"][cre_info["unique_id"] == query].values[0])
                    > 500
                ):
                    continue
                if (
                    float(
                        cre_info["length"][cre_info["unique_id"] == subject].values[0]
                    )
                    > 500
                ):
                    continue
                if query_name != subject_name:
                    out.write(f"{query}\t{subject}\t{similarity}\n")
                    processed_pairs.add((query_name, subject_name))
        logger.info(
            f"Finished processing {sp}.成功将 BLAST 结果转换为 MCL 输入格式，输出文件为 {mcl_input_file}"
        )
        subprocess.run(
            f"mcxload -abc {mcl_input_file} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o /media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.mci -write-tab /media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.tab",
            shell=True,
        )
        subprocess.run(
            f"mcl -te 4 /media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.mci -I 2 -use-tab /media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.tab -o /media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.mcl -V pruning",
            shell=True,
        )
    except FileNotFoundError:
        logger.error(
            f"FileNotFoundError: 未找到指定的 BLAST 结果文件: {blast_file_path}"
        )


@logger.catch
def find_popular(sp):
    # 为物种sp找到最受欢迎的聚类
    mcl_res = f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.mcl"
    mcl_res_df = pd.read_csv(mcl_res, sep="\t", header=None)
    mcl_res_df.columns = ["gene", "cluster"]
    cluster_count = mcl_res_df["cluster"].value_counts()
    popular_cluster = cluster_count.idxmax()
    popular_cluster_size = cluster_count.max()
    logger.info(
        f"物种 {sp} 中最受欢迎的聚类为 {popular_cluster}，包含基因数量为 {popular_cluster_size}"
    )
    ## 保存最受欢迎的聚类
    popular_cluster_genes = mcl_res_df[mcl_res_df["cluster"] == popular_cluster]
    popular_cluster_genes.to_csv(
        f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}_popular_cluster_genes.txt",
        sep="\t",
        index=False,
    )


@logger.catch
def process_sp(sp):
    logger.info(f"开始处理物种 {sp}")
    make_similarity_matrix(sp)
    # find_popular(sp)


if __name__ == "__main__":
    # 为每个物种创建相似性矩阵
    species = [
        "Ovis_aries",
        "Bos_taurus",
        # "Neophocaena_asiaeorientalis",
        "Sus_scrofa",
        # "Lama_glama",
        "Mustela_putorius",
        "Canis_lupus",
        "Felis_catus",
        # "Equus_asinus",
        "Equus_caballus",
        # "Rhinolophus_pusillus",
        "Myotis_chinensis",
        "Atelerix_albiventris",
        "Mus_musculus",
        "Rattus_norvegicus",
        # "Cavia_porcellus",
        "Oryctolagus_cuniculus",
        "Macaca_mulatta",
        "Rhinopithecus_roxellana",
        "Tupaia_belangeri",
        # "Procavia_capensis",
        "Petaurus_breviceps",
    ]
    logger.add(
        sink=sys.stdout,
        level="INFO",
        format="{time} {level} {message}",
        backtrace=True,
        diagnose=True,
    )
    with mp.Pool(processes=len(species)) as pool:
        pool.map(process_sp, species)

# import pytest
# import pandas as pd
# from enhancer_MCL import make_similarity_matrix

# # 测试正常情况
# def test_make_similarity_matrix_normal():
#     # 模拟数据
#     sp = "test_species"
#     cre_info = pd.DataFrame({
#         'unique_id': ['query1', 'query2', 'ubject1', 'ubject2'],
#         'length': [100, 200, 300, 400]
#     })

#     # 调用方法
#     make_similarity_matrix(sp)

#     # 断言文件是否生成
#     assert os.path.exists(f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}_mcl_input.txt")
#     assert os.path.exists(f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.mci")
#     assert os.path.exists(f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.tab")
#     assert os.path.exists(f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/{sp}.mcl")

# # 测试文件不存在情况
# def test_make_similarity_matrix_file_not_found():
#     sp = "test_species"
#     with pytest.raises(FileNotFoundError):
#         make_similarity_matrix(sp)

# # 测试其他异常情况
# def test_make_similarity_matrix_other_exception():
#     sp = "test_species"
#     # 模拟引发其他异常的情况
#     cre_info = pd.DataFrame({
#         'unique_id': ['query1', 'query2', 'ubject1', 'ubject2'],
#         'length': [100, 200, 300, 400]
#     })
#     with pytest.raises(Exception) as e:
#         make_similarity_matrix(sp)
#     assert str(e.value) == "模拟的其他异常"
