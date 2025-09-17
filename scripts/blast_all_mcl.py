import pandas as pd
import numpy as np
import multiprocessing as mp
from loguru import logger
import subprocess
import sys
from Bio import SeqIO
from functools import partial

"""
add species name to fasta record id and combine all fasta from all species to a single fasta 
then all against all blast, then mcl clustering
"""


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
def add_prefix_to_fasta(sp, output):
    """
    为fasta文件中的序列名称添加前缀 "sp_"。
    example: Mus_musculus_enh_1_B_Peak_1

    :param sp: 物种名称
    """
    # 读取fasta文件
    fasta_file = (
        f"/media/Data/zhangz/chip/analysis/{sp}/anno/blast2/fastadb/{sp}_all_cre.fa"
    )
    awk_command = f"awk '/^>/ {{print \">{sp}_\" substr($0, 2); next}} {{print}}' {fasta_file} >> {output}"
    subprocess.run(awk_command, shell=True, check=True)
    logger.info(f"rename and combine fasta file for {sp} done")


@logger.catch
def combine_fasta_files(
    species=species,
    output="/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre.fa",
):
    """
    将所有物种的fasta文件合并到一个文件中。

    :param species: 物种列表
    :param output: 输出文件路径
    """
    # remove old output file
    subprocess.run(f"rm -f {output}", shell=True, check=True)
    for sp in species:
        add_prefix_to_fasta(sp, output)


@logger.catch
def make_fastadb_and_blast(
    fasta="/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre.fa",
):
    """
    为fasta文件创建blast数据库并进行blast。

    :param fasta: fasta文件路径
    """
    # make blast db
    db_dir = fasta.replace(".fa", "_db")
    new_fa = f"{db_dir}/all_cre.fa"
    subprocess.run(f"mkdir -p {db_dir}", shell=True, check=True)
    subprocess.run(f"cd {db_dir}", shell=True, check=True)
    subprocess.run(f"mv {fasta} {db_dir}", shell=True, check=True)
    subprocess.run(
        f"makeblastdb -in {new_fa} -dbtype nucl -out all_cre",
        shell=True,
        check=True,
    )
    # blast
    subprocess.run(
        f"blastn -query {db_dir}/all_cre.fa -db {new_fa.replace('.fa', '')} -evalue 1e-5 -strand both -perc_identity 50 -task dc-megablast -outfmt 6 -out {new_fa.replace('.fa', '_blast.out')} -num_threads 30",
        shell=True,
        check=True,
    )
    logger.info(f"blast for {fasta} done")


@logger.catch
def make_similarity_matrix(
    species=species,
    blast_out="/media/Data/zhangz/chip/analysis/summary2/blast_all/all_cre_db/all_cre_blast.out",
):
    """
    为blast结果创建相似性矩阵。

    :param blast_out: blast结果文件路径
    """

    ## a dataframe contain unique_id from all species
    cre_ids = pd.DataFrame()
    for sp in species:
        cre_info = pd.read_csv(
            f"/media/Data/zhangz/chip/analysis/summary2/anno_sum/{sp}/{sp}_all.csv",
            sep=",",
            header=0,
        )
        cre_info["unique_id"] = cre_info.apply(
            lambda row: f"{row['unique_id']}_{row['peak']}", axis=1
        )
        cre_info["query_id"] = cre_info.apply(lambda row: f"{sp}_{row['id']}", axis=1)
        cre_ids = pd.concat(
            [cre_ids, cre_info[cre_info["length"] < 500][["unique_id", "query_id"]]]
        )

    ## filter cres by length and
    processed_pairs = set()
    mcl_input_file = blast_out.replace(".out", "_mcl_input.txt")

    try:
        with open(blast_out, "r") as f, open(mcl_input_file, "w") as out:
            for line in f:
                line = line.strip().split("\t")
                query = line[0]
                query_name = "_".join(query.split("_")[0:5])
                subject = line[1]
                subject_name = "_".join(subject.split("_")[0:5])
                similarity = line[2]
                if (query_name, subject_name) in processed_pairs:
                    continue
                if float(similarity) < 50:
                    continue
                if query_name not in cre_ids["query_id"].values:
                    continue
                if query_name != subject_name:
                    out.write(f"{query}\t{subject}\t{similarity}\n")
                    processed_pairs.add((query_name, subject_name))
        logger.info(f"make similarity matrix for {blast_out} done")
        mci_file = mcl_input_file.replace(".txt", ".mci")
        tab_file = mcl_input_file.replace(".txt", ".tab")
        mcl_file = mcl_input_file.replace("_input.txt", ".mcl")
        subprocess.run(
            f"mcxload -abc {mcl_input_file} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {mci_file} -write-tab {tab_file}",
            shell=True,
        )
        subprocess.run(
            f"mcl {mci_file} -I 2 -use-tab {tab_file} -o {mcl_file} -V pruning -te 30",
            shell=True,
        )
        # save processed_pairs
        with open(blast_out.replace(".out", "_processed_pairs.txt"), "w") as f:
            for pair in processed_pairs:
                f.write(f"{pair[0]}\t{pair[1]}\n")
    except FileNotFoundError:
        logger.error(f"FileNotFoundError: 未找到指定的 BLAST 结果文件: {blast_out}")


@logger.catch
def main():
    combine_fasta_files(species)
    make_fastadb_and_blast()
    make_similarity_matrix()


if __name__ == "__main__":
    main()
