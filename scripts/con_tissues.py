#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
we have mapped cres of each tissue and species to every another species genome,
and intersection has been made between orthologous sequence and cres of the same tissue
in target species. now we want to know how many orthologous sequence can be found
active in other tissues of target species.
we will still use intersectBed -f 0.5 -e to get the intersection between orthologous sequence
use concurrent.futures to run the intersectBed command in parallel
"""

# packages
import os
import concurrent.futures
import subprocess
import itertools
import pandas as pd
from loguru import logger


@logger.catch
def intersect(comb):
    index, row = comb
    sp1, sp2, tis1, tis2, ele = row
    halper_file = f"/media/Data/zhangz/chip/analysis/{sp1}/compare2/{tis1}/{ele}/{sp1}2{sp2}_{tis1}_{ele}_halper.bed4"
    target_file = f"/media/Data/zhangz/chip/analysis/{sp2}/anno/peaks/{sp2}_{tis2}_{ele}_optimal_abbr_sorted.narrowPeak"
    out_file = f"/media/Data/zhangz/chip/analysis/{sp1}/compare2/{tis1}/{ele}/{sp1}2{sp2}_{tis1}_{tis2}_{ele}_intersect.bed"
    cmd = (
        f"intersectBed -a {halper_file} -b {target_file} -f 0.5 -e -wa -u > {out_file}"
    )
    # run the command
    subprocess.run(cmd, shell=True)
    logger.info(f"{sp1}2{sp2}_{tis1}_{tis2}_{ele} done!")
    return out_file


@logger.catch
def count_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()
    return len(lines)


@logger.catch
def count_intersect(comb):
    index, row = comb
    sp1, sp2, tis1, tis2, ele = row
    in_file = f"/media/Data/zhangz/chip/analysis/{sp1}/compare2/{tis1}/{ele}/{sp1}2{sp2}_{tis1}_{tis2}_{ele}_intersect.bed"
    count = count_lines(in_file)
    return count


@logger.catch
def test(comb):
    index, row = comb
    # sp1 = comb['sp1']
    # sp2 = comb['sp2']
    # tis1 = comb['tis1']
    # tis2 = comb['tis2']
    # ele = comb['ele']
    sp1, sp2, tis1, tis2, ele = row
    print(sp1, sp2, tis1, tis2, ele)


@logger.catch
def tissue_intersect(comb):
    idx, row = comb
    sp1, sp2, tis1, ele = row
    num_bkl = None
    num_bk = None
    num_bl = None
    num_kl = None
    b_file = f"/media/Data/zhangz/chip/analysis/{sp1}/compare2/{tis1}/{ele}/{sp1}2{sp2}_{tis1}_Brain_{ele}_intersect.bed"
    k_file = f"/media/Data/zhangz/chip/analysis/{sp1}/compare2/{tis1}/{ele}/{sp1}2{sp2}_{tis1}_Kidney_{ele}_intersect.bed"
    l_file = f"/media/Data/zhangz/chip/analysis/{sp1}/compare2/{tis1}/{ele}/{sp1}2{sp2}_{tis1}_Liver_{ele}_intersect.bed"
    try:
        b_file = pd.read_csv(b_file, sep="\t", header=None)
        b_file.columns = ["chr", "start", "end", "name"]
    except pd.errors.EmptyDataError:
        b_file = pd.DataFrame()
        num_bkl = 0
        num_bk = 0
        num_bl = 0
    try:
        k_file = pd.read_csv(k_file, sep="\t", header=None)
        k_file.columns = ["chr", "start", "end", "name"]
    except pd.errors.EmptyDataError:
        k_file = pd.DataFrame()
        num_bkl = 0
        num_bk = 0
        num_kl = 0
    try:
        l_file = pd.read_csv(l_file, sep="\t", header=None)
        l_file.columns = ["chr", "start", "end", "name"]
    except pd.errors.EmptyDataError:
        l_file = pd.DataFrame()
        num_bkl = 0
        num_bl = 0
        num_kl = 0
    # count the number of intersection
    if num_bkl != 0:
        num_bkl = len(set(b_file["name"]) & set(k_file["name"]) & set(l_file["name"]))
    if num_bk != 0:
        num_bk = len(set(b_file["name"]) & set(k_file["name"]))
    if num_bl != 0:
        num_bl = len(set(b_file["name"]) & set(l_file["name"]))
    if num_kl != 0:
        num_kl = len(set(k_file["name"]) & set(l_file["name"]))
    res = pd.Series([num_bkl, num_bk, num_bl, num_kl], index=["bkl", "bk", "bl", "kl"])
    return res


@logger.catch
def main():
    species = [
        "Ovis_aries",
        "Bos_taurus",
        # "Neophocaena_asiaeorientalis",
        "Sus_scrofa",
        "Lama_glama",
        "Mustela_putorius",
        "Canis_lupus",
        "Felis_catus",
        "Equus_asinus",
        "Equus_caballus",
        "Rhinolophus_pusillus",
        "Rhinolophus_ferrumequinum",
        "Hipposideros_larvatus",
        "Myotis_ricketti",
        "Myotis_chinensis",
        "Atelerix_albiventris",
        "Mus_musculus",
        "Rattus_norvegicus",
        "Cavia_porcellus",
        "Oryctolagus_cuniculus",
        "Macaca_mulatta",
        "Procavia_capensis",
        "Tupaia_belangeri",
        # "Rhinopithecus_roxellana",
        "Petaurus_breviceps",
    ]
    tissues = ["Brain", "Kidney", "Liver"]
    elements = ["enhancer", "promoter"]
    combins = list(itertools.permutations(species, 2))
    res = list(itertools.product(combins, tissues, tissues, elements))
    res_df = pd.DataFrame(res, columns=["sps", "tis1", "tis2", "ele"])
    res_df[["sp1", "sp2"]] = res_df["sps"].apply(pd.Series)
    res_df = res_df.drop("sps", axis=1)
    res_df = res_df[["sp1", "sp2", "tis1", "tis2", "ele"]]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # 在进程池中运行函数
        results = executor.map(intersect, res_df.iterrows())

    with concurrent.futures.ProcessPoolExecutor() as executor:
        # combins = list(itertools.combinations(species, 2))
        results2 = executor.map(
            count_intersect, res_df[["sp1", "sp2", "tis1", "tis2", "ele"]].iterrows()
        )
    re_list = list(results2)
    # 将results2转换为一个Series
    results2_series = pd.Series(re_list)

    # 将这个Series作为一个新的列添加到res_df中
    res_df["count"] = results2_series
    # save the results
    res_df.to_csv(
        "/media/Data/zhangz/chip/analysis/conservation3/intersect_count.csv",
        index=False,
    )

    df2 = res_df.drop(["tis2", "count"], axis=1)
    df2 = df2.drop_duplicates()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results3 = executor.map(tissue_intersect, df2.iterrows())
    # for test
    # with concurrent.futures.ProcessPoolExecutor() as executor::
    #     results = executor.map(test, res_df.iloc[0:5].iterrows())
    # combins_df = pd.DataFrame(itertools.combinations(species, 2), columns = ['sp1', 'sp2'])
    # convert results3 into dataframe and merge with df2
    res3_list = list(results3)
    res3_df = pd.DataFrame(res3_list)
    res3_df = pd.concat([df2, res3_df], axis=1)
    res3_df.to_csv(
        "/media/Data/zhangz/chip/analysis/conservation3/tissue_intersect.csv",
        index=False,
    )


@logger.catch
def main2():
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
        "Rhinolophus_ferrumequinum",
        "Hipposideros_larvatus",
        "Myotis_ricketti",
        "Myotis_chinensis",
        "Atelerix_albiventris",
        "Mus_musculus",
        "Rattus_norvegicus",
        "Cavia_porcellus",
        "Oryctolagus_cuniculus",
        "Macaca_mulatta",
        "Procavia_capensis",
        "Tupaia_belangeri",
        "Rhinopithecus_roxellana",
        "Petaurus_breviceps",
    ]
    tissues = ["Brain", "Kidney", "Liver"]
    elements = ["enhancer", "promoter"]
    combins = list(itertools.permutations(species, 2))
    # jt_combins
    res = list(itertools.product(combins, tissues, tissues, elements))
    res_df = pd.DataFrame(res, columns=["sps", "tis1", "tis2", "ele"])
    res_df[["sp1", "sp2"]] = res_df["sps"].apply(pd.Series)
    res_df = res_df.drop("sps", axis=1)
    res_df = res_df[["sp1", "sp2", "tis1", "tis2", "ele"]]
    res_df = res_df[
        ~(
            (
                (res_df["sp1"] == "Neophocaena_asiaeorientalis")
                & (res_df["tis1"] != "Brain")
            )
            | (
                (res_df["sp2"] == "Neophocaena_asiaeorientalis")
                & (res_df["tis2"] != "Brain")
            )
            | (
                (res_df["sp1"] == "Rhinopithecus_roxellana")
                & (res_df["tis1"] == "Brain")
            )
            | (
                (res_df["sp2"] == "Rhinopithecus_roxellana")
                & (res_df["tis2"] == "Brain")
            )
        )
    ]
    res_df = res_df[
        (res_df["sp1"] == "Neophocaena_asiaeorientalis")
        | (res_df["sp2"] == "Neophocaena_asiaeorientalis")
        | (res_df["sp1"] == "Rhinopithecus_roxellana")
        | (res_df["sp2"] == "Rhinopithecus_roxellana")
    ]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # 在进程池中运行函数
        results = executor.map(intersect, res_df.iterrows())

    with concurrent.futures.ProcessPoolExecutor() as executor:
        # combins = list(itertools.combinations(species, 2))
        results2 = executor.map(
            count_intersect, res_df[["sp1", "sp2", "tis1", "tis2", "ele"]].iterrows()
        )
    re_list = list(results2)
    # 将results2转换为一个Series
    results2_series = pd.Series(re_list)

    # 将这个Series作为一个新的列添加到res_df中
    res_df["count"] = results2_series
    # save the results
    res_df.to_csv(
        "/media/Data/zhangz/chip/analysis/conservation3/intersect_count2.csv",
        index=False,
    )

    df2 = res_df.drop(["tis2", "count"], axis=1)
    df2 = df2.drop_duplicates()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results3 = executor.map(tissue_intersect, df2.iterrows())
    # for test
    # with concurrent.futures.ProcessPoolExecutor() as executor::
    #     results = executor.map(test, res_df.iloc[0:5].iterrows())
    # combins_df = pd.DataFrame(itertools.combinations(species, 2), columns = ['sp1', 'sp2'])
    # convert results3 into dataframe and merge with df2
    res3_list = list(results3)
    res3_df = pd.DataFrame(res3_list)
    res3_df = pd.concat([df2, res3_df], axis=1)
    res3_df.to_csv(
        "/media/Data/zhangz/chip/analysis/conservation3/tissue_intersect2.csv",
        index=False,
    )


if __name__ == "__main__":
    # main()
    main2()
