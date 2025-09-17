#! /usr/bin/env python
# coding = utf-8

import pandas as pd
import sys
import os
import numpy as np
import subprocess
from loguru import logger
import itertools
import sh
from remind import remind

"""
统计每个peak的比对情况以及在target基因组中的功能保守性，
输入信息为物种标号，每次执行一个物种
循环套层为：query物种-组织-元件-target物种
align统计标准：halper结果中出现元件名即为可比对
功能保守性标准：bedtools -e -r 0.5
halper结果需要转化为bed文件，取4列即可
"""


class QuerySp:
    # query species

    def __init__(self, SP_NUM):
        self.SPECIES = (
            "Canis_lupus",
            "Mustela_putorius",
            "Felis_catus",
            "Equus_asinus",
            "Equus_caballus",
            "Bos_taurus",
            "Ovis_aries",
            "Sus_scrofa",
            "Lama_glama",
            "Rhinolophus_pusillus",
            "Rhinolophus_ferrumequinum",
            "Hipposideros_larvatus",
            "Myotis_chinensis",
            "Myotis_ricketti",
            "Atelerix_albiventris",
            "Rattus_norvegicus",
            "Mus_musculus",
            "Cavia_porcellus",
            "Oryctolagus_cuniculus",
            "Macaca_mulatta",
            "Tupaia_belangeri",
            "Procavia_capensis",
            "Petaurus_breviceps",
            "Rhinopithecus_roxellana",
            "Neophocaena_asiaeorientalis",
        )
        self.chosen_sp = self.SPECIES[SP_NUM]
        logger.add(
            f"/media/Data/zhangz/chip/scripts2/log/align_con/align_con_{self.chosen_sp}.log",
            backtrace=True,
        )
        self.other_sp = list(self.SPECIES)
        self.other_sp.remove(self.chosen_sp)
        self.ELES = ("enhancer", "promoter")
        if self.chosen_sp == "Rhinopithecus_roxellana":
            self.TISSUES = ("Kidney", "Liver")
        elif self.chosen_sp == "Neophocaena_asiaeorientalis":
            self.TISSUES = ("Brain",)
        else:
            self.TISSUES = ("Kidney", "Liver", "Brain")

    @logger.catch
    def _align_consrv(self):
        @logger.catch
        def count_align(tis, ele):
            out_file = f"/media/Data/zhangz/chip/analysis/{self.chosen_sp}/compare2/{tis}/{ele}/{self.chosen_sp}_{tis}_{ele}_stats.tsv"
            query_bed = f"/media/Data/zhangz/chip/analysis/{self.chosen_sp}/anno/peaks/{self.chosen_sp}_{tis}_{ele}_optimal_abbr_sorted.narrowPeak"
            query_df = pd.read_csv(query_bed, sep="\t", header=None)
            peak_list = query_df[3].tolist()
            count_df = pd.DataFrame(
                np.zeros((len(peak_list), 24), np.int8),
                index=peak_list,
                columns=self.other_sp,
            )
            count_df.index.name = "peak"
            count_df["length"] = 0
            count_df["align"] = 0
            for sp in self.other_sp:
                target_file = f"/media/Data/zhangz/chip/analysis/{self.chosen_sp}/compare2/{tis}/{ele}/{self.chosen_sp}2{sp}_{tis}_{ele}_halper.bed"
                with open(target_file) as t:
                    target_peak_list = [line.strip().split("\t")[4] for line in t]
                    for peak in peak_list:
                        query_length = int(
                            query_df.loc[query_df[3] == peak, 2]
                            - query_df.loc[query_df[3] == peak, 1]
                            + 1
                        )
                        count_df.loc[peak, "length"] = query_length
                        if peak in target_peak_list:
                            count_df.loc[peak, sp] = 1
                            count_df.loc[peak, "align"] += 1
            count_df.to_csv(out_file, sep="\t", index=True)

        @logger.catch
        def count_consrv(tis, ele):
            out_file = f"/media/Data/zhangz/chip/analysis/{self.chosen_sp}/compare2/{tis}/{ele}/{self.chosen_sp}_{tis}_{ele}_stats.tsv"
            count_df = pd.read_csv(out_file, sep="\t", index_col=0)
            count_df["consrv"] = 0
            for sp in self.other_sp:
                if sp == "Neophocaena_asiaeorientalis" and tis in [
                    "Liver",
                    "Kidney",
                ]:
                    continue
                elif sp == "Rhinopithecus_roxellana" and (tis == "Brain"):
                    continue
                target_file = f"/media/Data/zhangz/chip/analysis/{self.chosen_sp}/compare2/{tis}/{ele}/{self.chosen_sp}2{sp}_{tis}_{ele}_halper.bed"
                target_bed = f"/media/Data/zhangz/chip/analysis/{self.chosen_sp}/compare2/{tis}/{ele}/{self.chosen_sp}2{sp}_{tis}_{ele}_halper.bed4"
                target_peak = f"/media/Data/zhangz/chip/analysis/{sp}/anno/peaks/{sp}_{tis}_{ele}_optimal_abbr_sorted.narrowPeak"
                subprocess.run(
                    f"cut -f 1,2,3,5 {target_file} > {target_bed}", shell=True
                )
                subprocess.run(
                    f"bedtools intersect -a {target_bed} -b {target_peak} -e -f 0.5 > {target_bed}.func",
                    shell=True,
                )
                count_df[sp + "_con"] = 0
                try:
                    func_df = pd.read_csv(target_bed + ".func", sep="\t", header=None)
                    for peak in count_df.index.tolist():
                        if peak in func_df[3].tolist():
                            count_df.loc[peak, sp + "_con"] = 1
                            count_df.loc[peak, "consrv"] += 1
                except pd.errors.EmptyDataError:
                    continue
            count_df.to_csv(out_file + ".con", sep="\t", index=True)

        for tis in self.TISSUES:
            for ele in self.ELES:
                logger.info(f"calculating alignable for {tis} {ele}")
                count_align(tis, ele)
                logger.info(f"calculating conservation for {tis} {ele}")
                count_consrv(tis, ele)


if __name__ == "__main__":
    SP_NUM = int(sys.argv[1])
    sp = QuerySp(SP_NUM)
    sp._align_consrv()
    logger.info(f"finished: align peaks of {sp.chosen_sp}")
    remind(
        f"finished: align peaks of {sp.chosen_sp}",
        f"finished: align peaks of {sp.chosen_sp}",
    )
