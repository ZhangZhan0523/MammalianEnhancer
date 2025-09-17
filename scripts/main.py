#!/usr/bin/env python
# coding=utf-8


import pandas as pd
import numpy as np
import sh
import os
import json
import sys
import subprocess
from loguru import logger


"""
把每个物种的每个组织定义成类，抗体、重复等作为对象的属性，
把每一个分析步骤定义成类的方法
再写一个实例化类的方法，重命名类，例如TEMP.__name__ = f'Pa_{name}'
然后根据数据表迭代调用该方法
"""


class SpTi:
    "每个物种每个组织的基类"
    # 初始路径设置可以使用类变量，但是其他不应当使用类变量
    # 把basicD改成实例变量 并且通过表格赋值
    # py_cwd = os.getcwd()
    # py_cwd = sys.path[0]
    # basicD = os.path.abspath(os.path.join(py_cwd, os.pardir)) #去掉路径里的..
    # clean_reads_dir = ''
    # chip_ac_dir = ''
    # chip_me_dir = ''
    # anno_dir = ''
    # ac_json = ""#directory
    # me_json = ""
    # ac_rep1_read1_ln = ''
    # ac_rep1_read2_ln = ''
    # ac_rep2_read1_ln = ''
    # ac_rep2_read2_ln = ''
    # me_rep1_read1_ln = ''
    # me_rep1_read2_ln = ''
    # me_rep2_read1_ln = ''
    # me_rep2_read2_ln = ''
    # inp_rep1_read1_ln = ''
    # inp_rep1_read2_ln = ''
    # inp_rep2_read1_ln = ''
    # inp_rep2_read2_ln = ''
    # """ac_rep1_read1 = 0, ac_rep1_read2 = 0, ac_rep2_read1 = 0, ac_rep2_read2 = 0,\
    #         me_rep1_read1 = 0, me_rep1_read2 = 0, me_rep2_read1 = 0, me_rep2_read2 = 0,\
    #             inp_rep1_read1 = 0, inp_rep1_read2 = 0, inp_rep2_read1 = 0, inp_rep1_read2 = 0,\
    # """

    @logger.catch
    def __init__(
        self,
        species,
        tissue,
        dd,
        ac=0,
        me=0,
        inp=0,
        ac_nm1=0,
        ac_nm2=0,
        me_nm1=0,
        me_nm2=0,
        inp_nm1=0,
        inp_nm2=0,
        ac_rep1_trimmed1=0,
        ac_rep1_trimmed2=0,
        ac_rep2_trimmed1=0,
        ac_rep2_trimmed2=0,
        me_rep1_trimmed1=0,
        me_rep1_trimmed2=0,
        me_rep2_trimmed1=0,
        me_rep2_trimmed2=0,
        inp_rep1_trimmed1=0,
        inp_rep1_trimmed2=0,
        inp_rep2_trimmed1=0,
        inp_rep2_trimmed2=0,
        json_ac=0,
        json_me=0,
        ac_peak=0,
        me_peak=0,
        genome_fa=0,
        genome_tsv=0,
        genome_gff=0,
        ana_dir=0,
    ):
        """
        初始化，传入参数，可以传入初始参数进行完整分析，也可以改变默认参数值，
        传入中间属性从而只执行部分方法。
        """
        self.species = species
        self.tissue = tissue
        self.ac = ac
        self.me = me
        self.inp = inp
        self.dd = dd
        self.ac_nm1 = ac_nm1
        self.ac_nm2 = ac_nm2
        self.me_nm1 = me_nm1
        self.me_nm2 = me_nm2
        self.inp_nm1 = inp_nm1
        self.inp_nm2 = inp_nm2
        # self.ac_rep1_read1 = ac_rep1_read1
        # self.ac_rep1_read2 = ac_rep1_read2
        # self.ac_rep2_read1 = ac_rep2_read1
        # self.ac_rep2_read2 = ac_rep2_read2
        # self.me_rep1_read1 = me_rep1_read1
        # self.me_rep1_read2 = me_rep1_read2
        # self.me_rep2_read1 = me_rep2_read1
        # self.me_rep2_read2 = me_rep2_read2
        # self.inp_rep1_read1 = inp_rep1_read1
        # self.inp_rep1_read2 = inp_rep1_read2
        # self.inp_rep2_read1 = inp_rep2_read1
        # self.inp_rep2_read2 = inp_rep2_read2
        self.ac_rep1_trimmed1 = ac_rep1_trimmed1
        self.ac_rep1_trimmed2 = ac_rep1_trimmed2
        self.ac_rep2_trimmed1 = ac_rep2_trimmed1
        self.ac_rep2_trimmed2 = ac_rep2_trimmed2
        self.me_rep1_trimmed1 = me_rep1_trimmed1
        self.me_rep1_trimmed2 = me_rep1_trimmed2
        self.me_rep2_trimmed1 = me_rep2_trimmed1
        self.me_rep2_trimmed2 = me_rep2_trimmed2
        self.inp_rep1_trimmed1 = inp_rep1_trimmed1
        self.inp_rep1_trimmed2 = inp_rep1_trimmed2
        self.inp_rep2_trimmed1 = inp_rep2_trimmed1
        self.inp_rep2_trimmed2 = inp_rep2_trimmed2
        self.json_ac = json_ac
        self.json_me = json_me
        self.ac_peak = ac_peak
        self.me_peak = me_peak
        self.genome_fa = genome_fa
        self.genome_tsv = genome_tsv
        self.genome_gff = genome_gff
        self.ana_dir = ana_dir
        self.clean_reads_dir = os.path.join(self.ana_dir, self.species, "reads")
        self.log_file = logger.add(
            os.path.join(
                self.ana_dir,
                "analysis",
                self.species,
                "log",
                self.species + "_" + self.tissue + "{time}.log",
            ),
            backtrace=True,
            diagnose=True,
            rotation="10 MB",
            encoding="utf-8",
            enqueue=True,
            level="INFO",
        )
        logger.info("初始化{0} {1}的分析".format(self.species, self.tissue))
        logger.info(
            "log文件位置是{0}".format(
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "log",
                    self.species + "_" + self.tissue + "{time}.log",
                )
            )
        )
        # print(self.species, self.tissue, "clean_reads_dir: ", self.clean_reads_dir)
        if self.ac != 0:
            logger.info("建立{0} {1}的H3K27ac分析分析文件夹".format(self.species, self.tissue))
            self.chip_ac_dir = os.path.join(
                self.ana_dir, "analysis", self.species, "chip", self.tissue, "H3K27ac"
            )
            # print(self.species, self.tissue, "chip_ac_dir: ", self.chip_ac_dir)
            logger.info(
                "{0} {1}的H3K27ac分析分析文件夹位置是{2}".format(
                    self.species, self.tissue, self.chip_ac_dir
                )
            )
        if self.me != 0:
            logger.info("建立{0} {1}的H3K4me3分析分析文件夹".format(self.species, self.tissue))
            self.chip_me_dir = os.path.join(
                self.ana_dir, "analysis", self.species, "chip", self.tissue, "H3K4me3"
            )
            # print(self.species, self.tissue, "chip_me_dir: ", self.chip_me_dir)
            logger.info(
                "{0} {1}的H3K4me3分析分析文件夹位置是{2}".format(
                    self.species, self.tissue, self.chip_me_dir
                )
            )
        self.anno_dir = os.path.join(self.ana_dir, "analysis", self.species, "anno")
        # print(self.species, self.tissue, "anno_dir: ", self.anno_dir)
        logger.info(
            "{0} {1}的注释文件夹位置是{2}".format(self.species, self.tissue, self.anno_dir)
        )
        self.ac_json = os.path.join(
            self.chip_ac_dir, self.species + "_" + self.tissue + "_H3K27ac.json"
        )
        self.me_json = os.path.join(
            self.chip_me_dir, self.species + "_" + self.tissue + "_H3K4me3.json"
        )
        self.ac_file = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [
                    self.species,
                    self.tissue,
                    "H3K27ac_overlap.optimal_peak.narrowPeak.gz",
                ]
            ),
        )
        self.me_file = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [
                    self.species,
                    self.tissue,
                    "H3K4me3_overlap.optimal_peak.narrowPeak.gz",
                ]
            ),
        )
        self.ac_gzip = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [
                    self.species,
                    self.tissue,
                    "H3K27ac_overlap.optimal_peak.narrowPeak.gz",
                ]
            ),
        )
        self.me_gzip = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [
                    self.species,
                    self.tissue,
                    "H3K4me3_overlap.optimal_peak.narrowPeak.gz",
                ]
            ),
        )
        self.ac_peak = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [self.species, self.tissue, "H3K27ac_overlap.optimal_peak.narrowPeak"]
            ),
        )
        self.me_peak = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [self.species, self.tissue, "H3K4me3_overlap.optimal_peak.narrowPeak"]
            ),
        )
        self.ac_qc_file = os.path.join(
            self.chip_ac_dir,
            "chipRes",
            "_".join([self.species, self.tissue, "H3K27ac_qc.json"]),
        )
        self.me_qc_file = os.path.join(
            self.chip_me_dir,
            "chipRes",
            "_".join([self.species, self.tissue, "H3K4me3_qc.json"]),
        )
        self.base_dir = os.path.join(self.ana_dir, "analysis", self.species)
        self.store_dir = "/media/Data/zhangz/usb1/chip/analysis"

    @logger.catch
    def data_prep(self):
        """
        在对应的文件夹中建立原始测序数据clean data的软链接
        """
        if not os.path.exists(self.clean_reads_dir):
            sh.mkdir(self.clean_reads_dir, p=True)
        for chip in [
            self.ac_nm1,
            self.ac_nm2,
            self.me_nm1,
            self.me_nm2,
            self.inp_nm1,
            self.inp_nm2,
        ]:
            if chip != "0":
                ln_dir = os.path.join(
                    self.ana_dir, "analysis", self.species, "reads", chip
                )
                if not os.path.exists(ln_dir):
                    sh.mkdir(ln_dir, p=True)  # 可以先sh.mkdir(ln_dir, p = True)
                try:
                    os.readlink(
                        os.path.join(
                            self.ana_dir,
                            "analysis",
                            self.species,
                            "reads",
                            chip,
                            chip + "_1.fq.gz",
                        )
                    )
                except FileNotFoundError:
                    sh.ln(
                        "-s",
                        os.path.join(self.dd, chip, chip + "_1.fq.gz"),
                        os.path.join(
                            self.ana_dir,
                            "analysis",
                            self.species,
                            "reads",
                            chip,
                            chip + "_1.fq.gz",
                        ),
                    )
                try:
                    os.readlink(
                        os.path.join(
                            self.ana_dir,
                            "analysis",
                            self.species,
                            "reads",
                            chip,
                            chip + "_2.fq.gz",
                        )
                    )
                except FileNotFoundError:
                    sh.ln(
                        "-s",
                        os.path.join(self.dd, chip, chip + "_2.fq.gz"),
                        os.path.join(
                            self.ana_dir,
                            "analysis",
                            self.species,
                            "reads",
                            chip,
                            chip + "_2.fq.gz",
                        ),
                    )
                    # print("make soft links for", chip, "paired reads")
                    logger.info("make soft links for {0} paired reads".format(chip))
        for dire in [
            self.chip_ac_dir,
            self.chip_me_dir,
            os.path.join(self.anno_dir, "peaks"),
        ]:
            if not os.path.exists(dire):
                sh.mkdir(dire, p=True)
                # print("Creating directories for chip ac, me, and annotations...")
                logger.info(
                    "Creating directories for chip H3K27ac, H3K4me3, and annotations..."
                )
                logger.info("now:\n{0}".format(dire))
        return f'data preparation for {self.species + "_" + self.tissue} done'
        # for ln_org, ln_des in zip([self.ac_rep1_read1, self.ac_rep1_read2, self.ac_rep2_read1, self.ac_rep2_read2, \
        #     self.me_rep1_read1, self.me_rep1_read2, self.me_rep2_read1, self.me_rep2_read2, \
        #         self.inp_rep1_read1, self.inp_rep1_read2, self.inp_rep2_read1, self.inp_rep2_read2], \
        #             [SpTi.ac_rep1_read1, SpTi.ac_rep1_read2, SpTi.ac_rep2_read1, SpTi.ac_rep2_read2, \
        #                 SpTi.me_rep1_read1, SpTi.me_rep1_read2, SpTi.me_rep2_read1, SpTi.me_rep2_read2, \
        #                     SpTi.inp_rep1_read1, SpTi.inp_rep1_read2, SpTi.inp_rep2_read1, SpTi.inp_rep2_read2]):
        #     ln_des = os.path.join()

    @logger.catch
    def trim(self):
        """
        trim_galore
        """
        for chip in [
            self.ac_nm1,
            self.ac_nm2,
            self.me_nm1,
            self.me_nm2,
            self.inp_nm1,
            self.inp_nm2,
        ]:
            trim_code = 0
            if chip != "0":
                trim_dir = os.path.join(
                    self.ana_dir, "analysis", self.species, "reads", chip, "trim"
                )
                if not os.path.exists(trim_dir):
                    sh.mkdir(trim_dir)
                else:
                    # print(self.species, chip, "trim文件夹已经存在")
                    logger.info("{0} {1} trim文件夹已经存在".format(self.species, chip))
                ln1 = os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    chip,
                    chip + "_1.fq.gz",
                )
                ln2 = os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    chip,
                    chip + "_2.fq.gz",
                )
                trim_cmd = " ".join(
                    [
                        "trim_galore",
                        "-q 20 -length 20 --stringency 3 -e 0.1 --gzip -j 4 -o",
                        trim_dir,
                        "-paired",
                        ln1,
                        ln2,
                        "--fastqc",
                    ]
                )
                logger.info("start trim_galore for {0} {1}".format(self.species, chip))
                logger.info("trim_galore command: {0}".format(trim_cmd))
                logger.info(
                    "trim_galore log file: {0}".format(
                        os.path.join(
                            self.ana_dir,
                            "analysis",
                            self.species,
                            "reads",
                            chip,
                            "trim",
                            "trim.log",
                        )
                    )
                )
                with open(
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        chip,
                        "trim",
                        "trim.log",
                    ),
                    "w",
                ) as trim_log:
                    trim_g = subprocess.run(
                        trim_cmd,
                        shell=True,
                        stderr=subprocess.STDOUT,
                        stdout=trim_log,
                        encoding="utf-8",
                    )
                logger.info(
                    "trim_galore for {0} {1} done, returncode is {2}".format(
                        self.species, chip, trim_g.returncode
                    )
                )
                if trim_g.returncode == 0:
                    trim_mes = " ".join(
                        [
                            self.species,
                            chip,
                            "Trim success, log file: ",
                            os.path.join(
                                self.ana_dir,
                                "analysis",
                                self.species,
                                "reads",
                                chip,
                                "trim",
                                "trim.log",
                            ),
                        ]
                    )
                else:
                    trim_mes = " ".join(
                        [
                            self.species,
                            chip,
                            "Trim failure, log file: ",
                            os.path.join(
                                self.ana_dir,
                                "analysis",
                                self.species,
                                "reads",
                                chip,
                                "trim",
                                "trim.log",
                            ),
                        ]
                    )
                trim_code += trim_g.returncode
        return (
            "all trimming has been successfully done..."
            if trim_code == 0
            else "something wrong with the trimming, please check your log file"
        )

    @logger.catch
    def make_json(self):
        """
        生成chip流程的json
        """
        if self.ac != 0:
            # ac_j_cont = {}
            # ac_j_cont["chip.title"] = " ".join([self.species, self.tissue, "H3K27ac"])
            ac_j_cont = {"chip.title": " ".join([self.species, self.tissue, "H3K27ac"])}
            ac_j_cont["chip.description"] = " ".join(
                ["call peaks for H3K27ac of", self.species, self.tissue]
            )
            ac_j_cont["chip.pipeline_type"] = "histone"
            ac_j_cont["chip.genome_tsv"] = self.genome_tsv
            ac_j_cont["chip.paired_end"] = True
            ac_j_cont["chip.ctl_paired_end"] = True
            ac_j_cont["chip.always_use_pooled_ctl"] = True
            ac_j_cont["chip.fastqs_rep1_R1"] = [
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    self.ac_nm1,
                    "trim",
                    self.ac_nm1 + "_1_val_1.fq.gz",
                )
            ]
            ac_j_cont["chip.fastqs_rep1_R2"] = [
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    self.ac_nm1,
                    "trim",
                    self.ac_nm1 + "_2_val_2.fq.gz",
                )
            ]
            ac_j_cont["chip.ctl_fastqs_rep1_R1"] = [
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    self.inp_nm1,
                    "trim",
                    self.inp_nm1 + "_1_val_1.fq.gz",
                )
            ]
            ac_j_cont["chip.ctl_fastqs_rep1_R2"] = [
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    self.inp_nm1,
                    "trim",
                    self.inp_nm1 + "_2_val_2.fq.gz",
                )
            ]
            if self.ac == 2:
                ac_j_cont["chip.fastqs_rep2_R1"] = [
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        self.ac_nm2,
                        "trim",
                        self.ac_nm2 + "_1_val_1.fq.gz",
                    )
                ]
                ac_j_cont["chip.fastqs_rep2_R2"] = [
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        self.ac_nm2,
                        "trim",
                        self.ac_nm2 + "_2_val_2.fq.gz",
                    )
                ]
            if self.inp == 2:
                ac_j_cont["chip.ctl_fastqs_rep2_R1"] = [
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        self.inp_nm2,
                        "trim",
                        self.inp_nm2 + "_1_val_1.fq.gz",
                    )
                ]
                ac_j_cont["chip.ctl_fastqs_rep2_R2"] = [
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        self.inp_nm2,
                        "trim",
                        self.inp_nm2 + "_2_val_2.fq.gz",
                    )
                ]
            # print(self.ac_json)
            logger.info(
                "start to write json file for {0} {1} H3K27ac: {2}".format(
                    self.species, self.tissue, self.ac_json
                )
            )
            # TODO: if self.chip_ac_json == "",exit code=2,error=""
            with open(self.ac_json, "w") as f:
                f.write(json.dumps(ac_j_cont, indent=4, separators=(",", " : ")))
        if self.me != 0:
            me_j_cont = {}
            me_j_cont["chip.title"] = " ".join([self.species, self.tissue, "H3K4me3"])
            me_j_cont["chip.description"] = " ".join(
                ["call peaks for H3K4me3 of ", self.species, self.tissue]
            )
            me_j_cont["chip.pipeline_type"] = "histone"
            me_j_cont["chip.genome_tsv"] = self.genome_tsv
            me_j_cont["chip.paired_end"] = True
            me_j_cont["chip.ctl_paired_end"] = True
            me_j_cont["chip.always_use_pooled_ctl"] = True
            me_j_cont["chip.fastqs_rep1_R1"] = [
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    self.me_nm1,
                    "trim",
                    self.me_nm1 + "_1_val_1.fq.gz",
                )
            ]
            me_j_cont["chip.fastqs_rep1_R2"] = [
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    self.me_nm1,
                    "trim",
                    self.me_nm1 + "_2_val_2.fq.gz",
                )
            ]
            me_j_cont["chip.ctl_fastqs_rep1_R1"] = [
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    self.inp_nm1,
                    "trim",
                    self.inp_nm1 + "_1_val_1.fq.gz",
                )
            ]
            me_j_cont["chip.ctl_fastqs_rep1_R2"] = [
                os.path.join(
                    self.ana_dir,
                    "analysis",
                    self.species,
                    "reads",
                    self.inp_nm1,
                    "trim",
                    self.inp_nm1 + "_2_val_2.fq.gz",
                )
            ]
            if self.me == 2:
                me_j_cont["chip.fastqs_rep2_R1"] = [
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        self.me_nm2,
                        "trim",
                        self.me_nm2 + "_1_val_1.fq.gz",
                    )
                ]
                me_j_cont["chip.fastqs_rep2_R2"] = [
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        self.me_nm2,
                        "trim",
                        self.me_nm2 + "_2_val_2.fq.gz",
                    )
                ]
            if self.inp == 2:
                me_j_cont["chip.ctl_fastqs_rep2_R1"] = [
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        self.inp_nm2,
                        "trim",
                        self.inp_nm2 + "_1_val_1.fq.gz",
                    )
                ]
                me_j_cont["chip.ctl_fastqs_rep2_R2"] = [
                    os.path.join(
                        self.ana_dir,
                        "analysis",
                        self.species,
                        "reads",
                        self.inp_nm2,
                        "trim",
                        self.inp_nm2 + "_2_val_2.fq.gz",
                    )
                ]
            # print(self.me_json)
            logger.info(
                "start to write json file for {0} {1} H3K4me3: {2}".format(
                    self.species, self.tissue, self.me_json
                )
            )
            with open(self.me_json, "w") as f:
                f.write(json.dumps(me_j_cont, indent=4, separators=(",", " : ")))
        return self.ac_json, self.me_json

    @logger.catch
    def call_peak(self):
        """
        执行encode的chip流程
        """
        if self.ac != 0:
            os.chdir(self.chip_ac_dir)
            ac_chip_cmd = " ".join(
                [
                    "caper run /media/tyloo/zhangz/chip-seq-pipeline2/chip.wdl -i",
                    self.ac_json,
                    "--singularity",
                ]
            )
            logger.info(
                "start to call peaks for {0} {1} H3K27ac: {2}".format(
                    self.species,
                    self.tissue,
                    os.path.join(self.chip_ac_dir, "chip.log"),
                )
            )
            with open(os.path.join(self.chip_ac_dir, "chip.log"), "w") as ac_chip_log:
                ac_chip = subprocess.run(
                    ac_chip_cmd,
                    shell=True,
                    stderr=subprocess.STDOUT,
                    stdout=ac_chip_log,
                    encoding="utf-8",
                )
            if ac_chip.returncode == 0:
                # print(self.species, self.tissue, "H3K27ac chip call peaks success, log file: ", os.path.join(self.chip_ac_dir, "chip.log"))
                logger.info(
                    "{0} {1} H3K27ac chip call peaks success, log file: {2}".format(
                        self.species,
                        self.tissue,
                        os.path.join(self.chip_ac_dir, "chip.log"),
                    )
                )
            else:
                # print(self.species, self.tissue, "H3K27ac chip call peaks failure, log file: ", os.path.join(self.chip_ac_dir, "chip.log"))
                logger.error(
                    "{0} {1} H3K27ac chip call peaks failure, log file: {2}".format(
                        self.species,
                        self.tissue,
                        os.path.join(self.chip_ac_dir, "chip.log"),
                    )
                )
        if self.me != 0:
            os.chdir(self.chip_me_dir)
            logger.info(
                "start to call peaks for {0} {1} H3K4me3: {2}".format(
                    self.species,
                    self.tissue,
                    os.path.join(self.chip_me_dir, "chip.log"),
                )
            )
            me_chip_cmd = " ".join(
                [
                    "caper run /media/tyloo/zhangz/chip-seq-pipeline2/chip.wdl -i",
                    self.me_json,
                    "--singularity",
                ]
            )
            with open(os.path.join(self.chip_me_dir, "chip.log"), "w") as me_chip_log:
                me_chip = subprocess.run(
                    me_chip_cmd,
                    shell=True,
                    stderr=subprocess.STDOUT,
                    stdout=me_chip_log,
                    encoding="utf-8",
                )
            if me_chip.returncode == 0:
                # print(self.species, self.tissue, "H3K4me3 chip call peaks success, log file: ", os.path.join(self.chip_me_dir, "chip.log"))
                logger.info(
                    "{0} {1} H3K4me3 chip call peaks success, log file: {2}".format(
                        self.species,
                        self.tissue,
                        os.path.join(self.chip_me_dir, "chip.log"),
                    )
                )
            else:
                # print(self.species, self.tissue, "H3K4me3 chip call peaks failure, log file: ", os.path.join(self.chip_me_dir, "chip.log"))
                logger.error(
                    "{0} {1} H3K4me3 chip call peaks failure, log file: {2}".format(
                        self.species,
                        self.tissue,
                        os.path.join(self.chip_me_dir, "chip.log"),
                    )
                )
        if ac_chip.returncode + me_chip.returncode == 0:
            call_mes = (
                "call peaks for "
                + self.species
                + " "
                + self.tissue
                + " successfully done!"
            )
        else:
            call_mes = (
                "call peaks for "
                + self.species
                + " "
                + self.tissue
                + " went wrong! Please check settings and log file."
            )
        return call_mes

    @logger.catch
    def org(self):
        """
        执行croo组织chip流程的结果
        """
        os.chdir(os.path.join(self.chip_ac_dir, "chip/*/"))
        sh.croo("metadata.json")
        os.chdir(os.path.join(self.chip_me_dir, "chip/*/"))
        sh.croo("metadata.json")

    @logger.catch
    def get_peak(self):
        """
        整理chip流程产生的peak
        """
        # TODO: 与volunme compression合并
        if self.ac != 0:
            os.chdir(os.path.join(self.anno_dir, "peaks"))
            ac_get_peak_cmd = " ".join(
                [
                    "find",
                    self.chip_ac_dir,
                    "-mindepth 6 -name overlap.optimal_peak.narrowPeak.gz ! -type l",
                ]
            )
            # ac_peak = sh.find(self.chip_ac_dir, mindepth = 6, name = "overlap.optimal_peak.narrowPeak.gz")
            ac_peak = subprocess.run(
                ac_get_peak_cmd, shell=True, stdout=subprocess.PIPE
            )
            # ac_file = os.path.join(self.anno_dir, "_".join([self.species, self.tissue, "H3K27ac_overlap.optimal_peak.narrowPeak.gz"]))
            ac_peak_out = ac_peak.stdout.decode()
            """
            TypeError: sequence item 1: expected str instance, bytes found
            """
            ac_cp_cmd = " ".join(["cp -f", ac_peak_out.rstrip(), self.ac_file])
            # overwrite the old output
            ac_cp = subprocess.run(ac_cp_cmd, shell=True, stdout=subprocess.PIPE)
            # if ac_cp.returncode == 0:
            #     print(" ".join(["H3K27ac peaks for", self.species, self.tissue, "successfully copied! The location is", \
            #         os.path.join(self.anno_dir, self.species + self.tissue + "_H3K27ac_overlap.optimal_peak.narrowPeak.gz")]))
            ac_gzip_cmd = " ".join(["gzip -d -f", self.ac_file])
            # overwrite the old output
            ac_gzip = subprocess.run(ac_gzip_cmd, shell=True, stdout=subprocess.PIPE)
            # if ac_gzip.returncode == 0:
            #     print("Decompressing successfully done.")
        if self.me != 0:
            os.chdir(os.path.join(self.anno_dir, "peaks"))
            me_get_peak_cmd = " ".join(
                [
                    "find",
                    self.chip_me_dir,
                    "-mindepth 6 -name overlap.optimal_peak.narrowPeak.gz ! -type l",
                ]
            )
            # me_peak = sh.find(self.chip_me_dir, mindepth = 6, name = "overlap.optimal_peak.narrowPeak.gz")
            me_peak = subprocess.run(
                me_get_peak_cmd, shell=True, stdout=subprocess.PIPE
            )
            # me_file = os.path.join(self.anno_dir, "_".join([self.species, self.tissue, "H3K4me3_overlap.optimal_peak.narrowPeak.gz"]))
            me_peak_out = me_peak.stdout.decode()
            me_cp_cmd = " ".join(["cp -f", me_peak_out.rstrip(), self.me_file])
            # overwrite the old output
            me_cp = subprocess.run(me_cp_cmd, shell=True, stdout=subprocess.PIPE)
            # if me_cp.returncode == 0:
            #     print(" ".join(["H3K4me3 peaks for", self.species, self.tissue, "successfully copied! The location is", \
            #         os.path.join(self.anno_dir, self.species + self.tissue + "_H3K4me3_overlap.optimal_peak.narrowPeak.gz")]))
            me_gzip_cmd = " ".join(["gzip -d -f", self.me_file])
            # overwrite
            subprocess.run(me_gzip_cmd, shell=True, stdout=subprocess.PIPE)
        #     if me_gzip.returncode == 0:
        #         print("Decompressing successfully done.")
        # returncode = ac_cp.returncode + ac_gzip.returncode + me_cp.returncode + me_gzip.returncode
        # if returncode == 0:
        #     get_peak_mes = "peaks successfully found."
        # else:
        get_peak_mes = "peaks successfully found."
        return get_peak_mes
        """
        例如：/media/usb1/chip/analysis/Myotis_chinensis/chip/
        Kidney/H3K4me3/chip/4bb4f490-ec93-4323-84fd-5fe2be230663/
        call-reproducibility_overlap/execution/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz
        """

    @logger.catch
    def r_seeker_anno(self):
        """
        使用ChIPSeeker进行peak注释
        """
        pass

    @logger.catch
    def find_motif_homer(self):
        """
        为在script所在位置找到的peak使用homer进行motif分析,包括denovo和known
        在"/media/Data/zhangz/chip/analysis"运行
        完成后转移回usb1
        """
        tmp_anno_dir = os.path.join(
            "/media/Data/zhangz/chip/analysis", self.species, "anno", self.tissue
        )
        subprocess.run(
            " ".join(["mkdir -p", tmp_anno_dir]), shell=True, stdout=subprocess.PIPE
        )
        os.chdir(tmp_anno_dir)
        if self.ac != 0:
            ac_peak_cp_cmd = " ".join(["cp -f", self.ac_file, tmp_anno_dir])
            ac_peak_cp = subprocess.run(
                ac_peak_cp_cmd, shell=True, stdout=subprocess.PIPE
            )
            ac_peak_tmp = os.path.join(
                tmp_anno_dir,
                "_".join(
                    [
                        self.species,
                        self.tissue,
                        "H3K27ac_overlap.optimal_peak.narrowPeak",
                    ]
                ),
            )
            ac_find_motif_out = os.path.join(
                tmp_anno_dir, "_".join([self.tissue, "H3K27ac_optimal"])
            )
            ac_find_motif_cmd = " ".join(
                [
                    "findMotifsGenome.pl",
                    ac_peak_tmp,
                    self.genome_fa,
                    ac_find_motif_out,
                    "-size given -mask -mset vertebrates -p 4",
                ]
            )
            ac_find_motif = subprocess.run(
                ac_find_motif_cmd, shell=True, stdout=subprocess.PIPE
            )
            ac_motif_cp_cmd = " ".join(["cp -r -f", ac_find_motif_out, self.anno_dir])
            ac_motif_cp = subprocess.run(
                ac_motif_cp_cmd, shell=True, stdout=subprocess.PIPE
            )
            ac_motif_tmp_rm_cmd = " ".join(["rm -r", ac_find_motif_out, ac_peak_tmp])
            ac_motif_tmp_rm = subprocess.run(
                ac_motif_tmp_rm_cmd, shell=True, stdout=subprocess.PIPE
            )
            print(
                "find motif for "
                + self.species
                + self.tissue
                + "H3K27ac_optimal successed"
            )
        if self.me != 0:
            me_peak_cp_cmd = " ".join(["cp -f", self.me_file, tmp_anno_dir])
            me_peak_cp = subprocess.run(
                me_peak_cp_cmd, shell=True, stdout=subprocess.PIPE
            )
            me_peak_tmp = os.path.join(
                tmp_anno_dir,
                "_".join(
                    [
                        self.species,
                        self.tissue,
                        "H3K4me3_overlap.optimal_peak.narrowPeak",
                    ]
                ),
            )
            me_find_motif_out = os.path.join(
                tmp_anno_dir, "_".join([self.tissue, "H3K4me3_optimal"])
            )
            me_find_motif_cmd = " ".join(
                [
                    "findMotifsGenome.pl",
                    me_peak_tmp,
                    self.genome_fa,
                    me_find_motif_out,
                    "-size given -mask -mset vertebrates -p 4",
                ]
            )
            me_find_motif = subprocess.run(
                me_find_motif_cmd, shell=True, stdout=subprocess.PIPE
            )
            me_motif_cp_cmd = " ".join(["cp -r -f", me_find_motif_out, self.anno_dir])
            me_motif_cp = subprocess.run(
                me_motif_cp_cmd, shell=True, stdout=subprocess.PIPE
            )
            me_motif_tmp_rm_cmd = " ".join(["rm -r", me_find_motif_out, me_peak_tmp])
            me_motif_tmp_rm = subprocess.run(
                me_motif_tmp_rm_cmd, shell=True, stdout=subprocess.PIPE
            )
            print(
                "find motif for "
                + self.species
                + self.tissue
                + "H3K4me3_optimal successed"
            )

    @logger.catch
    def data_keep_rm(self):
        """
        deal with files in /media/usb1/chip/ to save storage for batch 4
        """
        self.stoDir = os.path.join(
            self.ana_dir, "analysis", self.species, "chip", self.tissue
        )

        # if os.popen('hostname').read().strip() == "Genomics":
        #     self.destDir = "/media/usb1/chip/analysis"
        # elif os.popen('hostname').read().strip() == "Adaptation":
        #     self.destDir = "/media/Data/zhangz/usb1/chip/analysis"
        # self.keep_list = dict()
        # self.keepDir = os.path.join()
        def keepFile(ana_dir, species, tissue, mod, rep, chip1, inp1):
            keep_list = {}
            res_dir = os.path.join(self.stoDir, mod, "chipRes")
            if not os.path.exists(res_dir):
                sh.mkdir(res_dir, p=True)
            mvCode = 0
            resFiles = (
                "qc.html",
                "qc.json",
                "overlap.conservative_peak.narrowPeak.gz",
                "overlap.conservative_peak.narrowPeak.bb",
                "overlap.conservative_peak.narrowPeak.hammock.gz",
                "overlap.conservative_peak.narrowPeak.starch",
                "overlap.optimal_peak.narrowPeak.gz",
                "overlap.optimal_peak.narrowPeak.bb",
                "overlap.optimal_peak.narrowPeak.hammock.gz",
                "overlap.optimal_peak.narrowPeak.starch",
                "rep.pooled_x_ctl.pooled.fc.signal.bigwig",
                "rep.pooled_x_ctl.pooled.pval.signal.bigwig",
            )
            for fName in resFiles:
                if fName in [
                    "rep.pooled_x_ctl.pooled.fc.signal.bigwig",
                    "rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                ]:
                    if rep == 1:
                        fName = (
                            "_".join(
                                [
                                    chip1,
                                    "1_val_1.srt.nodup_x",
                                    inp1,
                                    "1_val_1.srt.nodup.fc.signal.bigwig",
                                ]
                            )
                            if fName == "rep.pooled_x_ctl.pooled.fc.signal.bigwig"
                            else "_".join(
                                [
                                    chip1,
                                    "1_val_1.srt.nodup_x",
                                    inp1,
                                    "1_val_1.srt.nodup.pval.signal.bigwig",
                                ]
                            )
                        )
                        maxDepth = "5"
                else:
                    maxDepth = "4"
                # print("find maxdepth is "+maxDepth)
                logger.info("find maxdepth is " + maxDepth)
                findCmd = " ".join(
                    [
                        "find",
                        os.path.join(self.stoDir, mod, "chip"),
                        "-mindepth 3 -maxdepth",
                        maxDepth,
                        "-name",
                        fName,
                        "! -type l",
                    ]
                )
                findRe = subprocess.run(
                    findCmd,
                    shell=True,
                    stderr=subprocess.STDOUT,
                    stdout=subprocess.PIPE,
                    encoding="utf-8",
                )
                # print(findRe)
                logger.info(findRe)
                if findRe.returncode == 0 and len(findRe.stdout.strip()) != 0:
                    keep_list[fName] = findRe.stdout.strip()
                    # print(" ".join(['successfully found', species, tissue, mod, fName, ':', keep_list[fName]]))
                    logger.info(
                        " ".join(
                            [
                                "successfully found",
                                species,
                                tissue,
                                mod,
                                fName,
                                ": \n",
                                keep_list[fName],
                            ]
                        )
                    )
                    mvCmd = " ".join(
                        [
                            "mv",
                            keep_list[fName],
                            os.path.join(
                                res_dir, "_".join([species, tissue, mod, fName])
                            ),
                        ]
                    )
                    mvRes = subprocess.run(
                        mvCmd,
                        shell=True,
                        stderr=subprocess.STDOUT,
                        stdout=subprocess.PIPE,
                        encoding="utf-8",
                    )
                    # print(mvRes)
                    logger.info(mvRes)
                    mvCode += mvRes.returncode
                else:
                    mvCode += 1
            failCount = 0
            for f in keep_list.keys():
                # print(" ".join(['checking if', os.path.join(res_dir, "_".join([species, tissue, mod, f])), 'exist']))
                logger.info(
                    " ".join(
                        [
                            "checking if",
                            os.path.join(res_dir, "_".join([species, tissue, mod, f])),
                            "exist",
                        ]
                    )
                )
                if not os.path.exists(
                    os.path.join(res_dir, "_".join([species, tissue, mod, f]))
                ):
                    failCount += 1
                    # print(" ".join([species, tissue, mod, f, 'cannot be found in', res_dir, '!!! please check manually.']))
                    logger.error(
                        " ".join(
                            [
                                species,
                                tissue,
                                mod,
                                f,
                                "cannot be found in",
                                res_dir,
                                "!!! please check manually.",
                            ]
                        )
                    )
                else:
                    # print(os.path.join(res_dir, f)+' exists')
                    logger.info(os.path.join(res_dir, f) + " exists")
            res = "fail"
            if failCount == 0:
                rmRes = subprocess.run(
                    " ".join(["rm -r", os.path.join(self.stoDir, mod, "chip")]),
                    shell=True,
                    stderr=subprocess.STDOUT,
                    stdout=subprocess.PIPE,
                    encoding="utf-8",
                )
                if rmRes.returncode == 0:
                    res = " ".join(
                        [
                            "File selection and rm successfully done for",
                            species,
                            tissue,
                            mod,
                        ]
                    )
            else:
                dirTree = subprocess.run(
                    " ".join(["tree -L 3", os.path.join(self.stoDir, mod)]),
                    shell=True,
                    stderr=subprocess.STDOUT,
                    stdout=subprocess.PIPE,
                    encoding="utf-8",
                )
                res = (
                    "something terribly wrong, please check the file tree of "
                    + dirTree.stdout
                )
            return res

        if self.ac != 0:
            # print(keepFile(self.ana_dir,self.species,self.tissue,'H3K27ac',self.ac, self.ac_nm1, self.inp_nm1))
            logger.info(
                keepFile(
                    self.ana_dir,
                    self.species,
                    self.tissue,
                    "H3K27ac",
                    self.ac,
                    self.ac_nm1,
                    self.inp_nm1,
                )
            )
        if self.me != 0:
            # print(keepFile(self.ana_dir,self.species,self.tissue,'H3K4me3',self.me, self.me_nm1, self.inp_nm1))
            logger.info(
                keepFile(
                    self.ana_dir,
                    self.species,
                    self.tissue,
                    "H3K4me3",
                    self.me,
                    self.me_nm1,
                    self.inp_nm1,
                )
            )
        keepFileRes = " ".join(
            [
                "result files keeping and removement done for",
                self.species,
                self.tissue,
                "but there might be errors, please check log above carefully",
            ]
        )
        return keepFileRes

    @logger.catch
    def gzip_peak(self):
        def mv_peak(peak_file):
            mvCmd = " ".join(["mv", peak_file, peak_file + ".gz"])
            mvRes = subprocess.run(
                mvCmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            return mvRes.returncode

        def gz_peak(peak_file):
            gzipCmd = " ".join(["gzip -d -f", peak_file])
            gzipRes = subprocess.run(
                gzipCmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            return gzipRes.returncode

        for peak in self.ac_peak, self.me_peak:
            if os.path.exists(peak):
                if mv_peak(peak) == 0:
                    logger.info(" ".join([peak, "moved successfully"]))
                else:
                    logger.error(
                        " ".join([peak, "cannot be moved, please check manually"])
                    )
                if gz_peak(peak) == 0:
                    logger.info(" ".join([peak, "gzipped successfully"]))
            else:
                logger.error(" ".join([peak, "cannot be found, please check manually"]))

    @logger.catch
    def merge_peak(self):
        def sort_peak(peak_file):
            sortCmd = " ".join(
                ["sort -k1,1 -k2,2n", peak_file, ">", peak_file + ".sorted"]
            )
            sortRes = subprocess.run(
                sortCmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            return sortRes.returncode

        for peak in self.ac_peak, self.me_peak:
            if os.path.exists(peak):
                if sort_peak(peak) == 0:
                    logger.info(" ".join([peak, "sorted successfully"]))
                else:
                    logger.error(
                        " ".join([peak, "cannot be sorted, please check manually"])
                    )
            else:
                logger.error(" ".join([peak, "cannot be found, please check manually"]))
        self.intersectPeak = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [
                    self.species,
                    self.tissue,
                    "_overlap_optimal_intersect_peak.narrowPeak",
                ]
            ),
        )
        intersectCmd = " ".join(
            [
                "bedtools intersect -a",
                self.ac_peak + ".sorted",
                "-b",
                self.me_peak + ".sorted",
                "-f 0.5 -r -wa -wb",
                ">",
                self.intersectPeak,
            ]
        )
        intersectRes = subprocess.run(
            intersectCmd,
            shell=True,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        if intersectRes.returncode == 0:
            logger.info(
                " ".join(["intersect peak done for", self.species, self.tissue])
            )
        else:
            logger.error(
                " ".join(["intersect peak failed for", self.species, self.tissue])
            )
        self.overlapPeak = pd.read_csv(self.intersectPeak, sep="\t", header=None)
        self.overlapPeak.columns = [
            "chr1",
            "start1",
            "end1",
            "name1",
            "score1",
            "strand1",
            "signalValue1",
            "pValue1",
            "qValue1",
            "peak1",
            "chr2",
            "start2",
            "end2",
            "name2",
            "score2",
            "strand2",
            "signalValue2",
            "pValue2",
            "qValue2",
            "peak2",
        ]
        self.overlapPeak["chr"] = self.overlapPeak["chr1"]
        self.overlapPeak["start"] = self.overlapPeak[["start1", "start2"]].min(axis=1)
        self.overlapPeak["end"] = self.overlapPeak[["end1", "end2"]].max(axis=1)
        self.overlapPeak["score"] = (
            self.overlapPeak["score1"].astype(str)
            + "."
            + self.overlapPeak["score2"].astype(str)
        )
        self.overlapPeak["score"] = self.overlapPeak["score"].astype(float)
        self.overlapPeak["name"] = (
            self.overlapPeak["name1"] + "_" + self.overlapPeak["name2"]
        )
        self.overlapPeak["strand"] = self.overlapPeak["strand1"]
        self.overlapMerge = pd.DataFrame(
            {
                "chr": self.overlapPeak["chr"],
                "start": self.overlapPeak["start"],
                "end": self.overlapPeak["end"],
                "name": self.overlapPeak["name"],
                "score": self.overlapPeak["score"],
                "strand": self.overlapPeak["strand"],
            }
        )
        self.actPromPeak = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [self.species, self.tissue, "_overlap_optimal_actProm.narrowPeak"]
            ),
        )
        self.overlapMerge.to_csv(self.actPromPeak, sep="\t", index=False, header=False)
        self.promoterPeak = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [self.species, self.tissue, "_overlap_optimal_promoter.narrowPeak"]
            ),
        )
        self.enhancerPeak = os.path.join(
            self.anno_dir,
            "peaks",
            "_".join(
                [self.species, self.tissue, "_overlap_optimal_enhancer.narrowPeak"]
            ),
        )
        enhancerCmd = " ".join(
            [
                "bedtools intersect -a",
                self.ac_peak + ".sorted",
                "-b",
                self.me_peak + ".sorted",
                "-v -f 0.5 -r",
                ">",
                self.enhancerPeak,
            ]
        )
        enhancerRes = subprocess.run(
            enhancerCmd,
            shell=True,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        if enhancerRes.returncode == 0:
            logger.info(" ".join(["get enhancer done for", self.species, self.tissue]))
        else:
            logger.error(
                " ".join(["get enhancer failed for", self.species, self.tissue])
            )
        promoterCmd = " ".join(
            [
                "bedtools intersect -a",
                self.me_peak + ".sorted",
                "-b",
                self.ac_peak + ".sorted",
                "-v -f 0.5 -r",
                ">",
                self.promoterPeak,
            ]
        )
        promoterRes = subprocess.run(
            promoterCmd,
            shell=True,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        if promoterRes.returncode == 0:
            logger.info(" ".join(["get promoter done for", self.species, self.tissue]))
        else:
            logger.error(
                " ".join(["get promoter failed for", self.species, self.tissue])
            )

    @logger.catch
    def actProm_peak(self):
        def sort_peak(peak_file):
            sortCmd = " ".join(
                ["sort -k1,1 -k2,2n", peak_file, ">", peak_file + ".sorted"]
            )
            sortRes = subprocess.run(
                sortCmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            return sortRes.returncode

        for peak in self.ac_peak, self.me_peak:
            if os.path.exists(peak):
                if sort_peak(peak) == 0:
                    logger.info(" ".join([peak, "sorted successfully"]))
                else:
                    logger.error(
                        " ".join([peak, "cannot be sorted, please check manually"])
                    )
            else:
                logger.error(" ".join([peak, "cannot be found, please check manually"]))
        self.intersectPeak = os.path.join(
            "/media/Data/zhangz/usb1/chip/analysis",
            self.species,
            "anno/peaks",
            "_".join(
                [
                    self.species,
                    self.tissue,
                    "overlap_optimal_intersect_peak.narrowPeak",
                ]
            ),
        )
        intersectCmd = " ".join(
            [
                "bedtools intersect -a",
                self.ac_peak + ".sorted",
                "-b",
                self.me_peak + ".sorted",
                "-f 0.5 -r -wa -wb",
                ">",
                self.intersectPeak,
            ]
        )
        intersectRes = subprocess.run(
            intersectCmd,
            shell=True,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        if intersectRes.returncode == 0:
            logger.info(
                " ".join(["intersect peak done for", self.species, self.tissue])
            )
        else:
            logger.error(
                " ".join(["intersect peak failed for", self.species, self.tissue])
            )
        self.overlapPeak = pd.read_csv(self.intersectPeak, sep="\t", header=None)
        self.overlapPeak.columns = [
            "chr1",
            "start1",
            "end1",
            "name1",
            "score1",
            "strand1",
            "signalValue1",
            "pValue1",
            "qValue1",
            "peak1",
            "chr2",
            "start2",
            "end2",
            "name2",
            "score2",
            "strand2",
            "signalValue2",
            "pValue2",
            "qValue2",
            "peak2",
        ]
        self.overlapPeak["chr"] = self.overlapPeak["chr1"]
        self.overlapPeak["start"] = self.overlapPeak[["start1", "start2"]].min(axis=1)
        self.overlapPeak["end"] = self.overlapPeak[["end1", "end2"]].max(axis=1)
        self.overlapPeak["score"] = (
            self.overlapPeak["score1"].astype(str)
            + "."
            + self.overlapPeak["score2"].astype(str)
        )
        self.overlapPeak["score"] = self.overlapPeak["score"].astype(float)
        self.overlapPeak["name"] = (
            self.overlapPeak["name1"] + "_" + self.overlapPeak["name2"]
        )
        self.overlapPeak["strand"] = self.overlapPeak["strand1"]
        self.overlapPeak["foldchange"] = (
            self.overlapPeak["signalValue1"].astype(float)
            + self.overlapPeak["signalValue2"].astype(float)
        ) / 2
        self.overlapMerge = pd.DataFrame(
            {
                "chr": self.overlapPeak["chr"],
                "start": self.overlapPeak["start"],
                "end": self.overlapPeak["end"],
                "name": self.overlapPeak["name"],
                "score": self.overlapPeak["score"],
                "strand": self.overlapPeak["strand"],
                "foldchange": self.overlapPeak["foldchange"],
            }
        )
        self.actPromPeak = os.path.join(
            "/media/Data/zhangz/usb1/chip/analysis/",
            self.species,
            "anno/peaks",
            "_".join([self.species, self.tissue, "overlap_optimal_actProm.narrowPeak"]),
        )
        self.overlapMerge.to_csv(self.actPromPeak, sep="\t", index=False, header=False)

    @logger.catch
    def qc_sum(self):
        if os.path.exist(self.ac_qc_file):
            pass

    @logger.catch
    def cp_all(self):
        cp_all_cmd = " ".join(["cp -r -n", self.base_dir, self.store_dir])
        cp_all_res = subprocess.run(
            cp_all_cmd,
            shell=True,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        if cp_all_res.returncode == 0:
            logger.info(" ".join(["cp all files done for", self.species, self.tissue]))
        else:
            logger.error(
                " ".join(["cp all files failed for", self.species, self.tissue])
            )

    @logger.catch
    def meme_pipe(self):
        """
        使用meme-chip进行motif分析
        """

        def bed2fasta(species, tissue, peak_type, genome_fa):
            """
            将bed文件转换为fasta文件
            """
            logger.info(
                " ".join(["starting bed converting for", species, tissue, peak_type])
            )
            bind_fasta_dir = os.path.join(self.store_dir, species, "anno/peaks")
            bind_motif_dir = os.path.join(self.store_dir, species, "anno/motifs")
            if not os.path.exists(bind_motif_dir):
                sh.mkdir(bind_motif_dir, p=True)
                logger.info(" ".join(["mkdir -p", bind_motif_dir]))
            bind_genome_dir = "/media/Data/zhangz/chip/genomes"
            meme_sif = "/media/Data/zhangz/chip/scripts/meme.sif"
            bind_list = ",".join(
                [
                    ":".join([bind_fasta_dir, "/data"]),
                    ":".join([bind_motif_dir, "/out"]),
                    ":".join([bind_genome_dir, "/genomes"]),
                ]
            )
            out_fasta = os.path.join(
                "/out", "_".join([species, tissue, peak_type]) + ".fasta"
            )
            peak_file = os.path.join(
                "/data",
                "_".join([species, tissue, "_overlap_optimal", peak_type])
                + ".narrowPeak",
            )
            genome_file = genome_fa.replace(
                "/media/Data/zhangz/chip/genomes/", "/genomes/"
            )
            b2f_cmd = " ".join(
                [
                    "singularity exec -B",
                    bind_list,
                    meme_sif,
                    "bed2fasta",
                    "-both -o ",
                    out_fasta,
                    peak_file,
                    genome_file,
                ]
            )
            logger.info(" ".join(["bed2fasta for", species, tissue, peak_type]))
            b2f_res = subprocess.run(
                b2f_cmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            logger.info(b2f_res.stdout)
            if b2f_res.returncode == 0:
                logger.info(
                    " ".join(["bed2fasta done for", species, tissue, peak_type])
                )
            else:
                logger.error(
                    " ".join(["bed2fasta failed for", species, tissue, peak_type])
                )

        # bed2fasta
        def meme_chip(species, tissue, peak_type):
            """
            使用meme-chip进行motif分析
            """
            logger.info(
                " ".join(["starting motif analysis for", species, tissue, peak_type])
            )
            bind_data_dir = os.path.join(self.store_dir, species, "anno/motifs")
            bind_genome_dir = "/media/Data/zhangz/chip/genomes"
            meme_sif = "/media/Data/zhangz/chip/scripts/meme.sif"
            bind_list = ",".join(
                [
                    ":".join([bind_data_dir, "/data"]),
                    ":".join([bind_genome_dir, "/genomes"]),
                ]
            )
            out_dir = os.path.join("/data", "_".join([species, tissue, peak_type]))
            peak_fa_file = os.path.join(
                "/data", "_".join([species, tissue, peak_type]) + ".fasta"
            )
            genome_file = os.path.join("/genomes", species, species + ".fa")
            log_file = os.path.join(
                bind_data_dir,
                "_".join(["meme_chip", species, tissue, peak_type]) + ".log",
            )
            logger.info("log file is " + log_file)
            # if not os.path.exists(out_dir):
            #     sh.mkdir(out_dir, p = True)
            #     logger.info(" ".join(['mkdir -p', out_dir]))
            meme_cmd = " ".join(
                [
                    "singularity exec -B",
                    bind_list,
                    meme_sif,
                    "meme-chip -o",
                    out_dir,
                    "-dna -meme-p 6 -db /genomes/motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt -maxw 20 -meme-nmotifs 10",
                    peak_fa_file,
                    ">",
                    log_file,
                ]
            )
            meme_res = subprocess.run(
                meme_cmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            logger.info(meme_res.stdout)
            if meme_res.returncode == 0:
                logger.info(
                    " ".join(["meme-chip done for", species, tissue, peak_type])
                )
            else:
                logger.error(
                    " ".join(["meme-chip failed for", species, tissue, peak_type])
                )

        # singularity meme-chip
        for peak_type in ["enhancer", "promoter", "actProm"]:
            bed2fasta(self.species, self.tissue, peak_type, self.genome_fa)
            meme_chip(self.species, self.tissue, peak_type)

    @logger.catch
    def meme_fimo(self):
        """
        使用fimo进行已知motif扫描
        """

        def fimo(species, tissue, peak_type):
            data_dir = os.path.join(self.store_dir, species, "anno/motifs")
            ana_dir = os.path.join(
                "/media/Data/zhangz/analysis/", species, "anno/motifs"
            )
            sh.mkdir(ana_dir, p=True)
            if not os.path.exists(
                os.path.join(ana_dir, "_".join([species, tissue, peak_type]) + ".fasta")
            ):
                sh.cp(
                    os.path.join(
                        data_dir, "_".join([species, tissue, peak_type]) + ".fasta"
                    ),
                    ana_dir,
                )
            genome_dir = "/media/Data/zhangz/chip/genomes"
            out_dir = os.path.join(
                "/data", "_".join([species, tissue, peak_type, "fimo"])
            )
            bind_list = ",".join(
                [
                    ":".join([ana_dir, "/data"]),
                    ":".join([genome_dir, "/genomes"]),
                    "/media/Data/zhangz/tmp:/media/Data/zhangz/tmp",
                ]
            )
            meme_sif = "/media/Data/zhangz/chip/scripts/meme.sif"
            input_fasta = os.path.join(
                "/data", "_".join([species, tissue, peak_type]) + ".fasta"
            )
            jaspar_file = "/genomes/motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
            log_file = os.path.join(
                ana_dir, "_".join(["fimo", species, tissue, peak_type]) + ".log"
            )
            fimo_cmd = " ".join(
                [
                    "singularity exec -B",
                    bind_list,
                    meme_sif,
                    "fimo --parse-genomic-coord -oc",
                    out_dir,
                    jaspar_file,
                    input_fasta,
                ]
            )
            logger.info(" ".join(["starting FIMO for", species, tissue, peak_type]))
            meme_res = subprocess.run(
                fimo_cmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            with open(log_file, "w") as logf:
                print(meme_res.stdout, logf)
            logger.info("log file is " + log_file)
            if meme_res.returncode == 0:
                logger.info(" ".join(["FIMO done for", species, tissue, peak_type]))
                sh.cp(
                    os.path.join(
                        ana_dir, "_".join([species, tissue, peak_type, "fimo"])
                    ),
                    os.path.join(
                        data_dir, "_".join([species, tissue, peak_type, "fimo"])
                    ),
                    r=True,
                )
                sh.cp(log_file, os.path.join(data_dir, log_file))
            else:
                logger.error(" ".join(["FIMO failed for", species, tissue, peak_type]))

        # singularity fimo
        for peak_type in ["enhancer", "promoter", "actProm"]:
            fimo(self.species, self.tissue, peak_type)

    @logger.catch
    def chipseeker(self):
        """
        使用chipseeker进行gene注释
        """
        anno_dir = os.path.join(self.store_dir, self.species, "anno/peaks")
        enhancer_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_enhancer"])
            + ".narrowPeak",
        )
        promoter_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_promoter"])
            + ".narrowPeak",
        )
        actProm_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_actProm"])
            + ".narrowPeak",
        )
        skr_cmd = " ".join(
            [
                "Rscript",
                "/media/Data/zhangz/chip/scripts/clusterPeak.R -s",
                self.species,
                "-t",
                self.tissue,
                "-g",
                self.genome_gff,
                "-w",
                anno_dir,
                "-e",
                enhancer_file,
                "-p",
                promoter_file,
                "-a",
                actProm_file,
            ]
        )
        skr_res = subprocess.run(
            skr_cmd,
            shell=True,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        logger.info(" ".join(["starting chipseeker for", self.species, self.tissue]))
        skr_log = os.path.join(
            anno_dir, "_".join(["chipseeker", self.species, self.tissue]) + ".log"
        )
        logger.info("log file: " + skr_log)
        with open(skr_log, "w") as skr_log_file:
            print(skr_res.stdout, file=skr_log_file)
        if skr_res.returncode == 0:
            logger.info(" ".join(["chipseeker done for", self.species, self.tissue]))
        else:
            logger.error(" ".join(["chipseeker failed for", self.species, self.tissue]))

    @logger.catch
    def specific_chipseeker(self):
        """
        使用chipseeker进行gene注释
        """
        anno_dir = os.path.join(
            self.store_dir, self.species, "anno/peaks", self.tissue + "specific"
        )
        enhancer_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_enhancer_specific"])
            + ".narrowPeak",
        )
        promoter_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_promoter_specific"])
            + ".narrowPeak",
        )
        actProm_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_actProm_specific"])
            + ".narrowPeak",
        )
        skr_cmd = " ".join(
            [
                "Rscript",
                "/media/Data/zhangz/chip/scripts/clusterSpecific.R -s",
                self.species,
                "-t",
                self.tissue,
                "-g",
                self.genome_gff,
                "-w",
                anno_dir,
                "-e",
                enhancer_file,
                "-p",
                promoter_file,
                "-a",
                actProm_file,
            ]
        )
        skr_log = os.path.join(
            anno_dir,
            "_".join(["tissue_specific_chipseeker", self.species, self.tissue])
            + ".log",
        )
        logger.info(" ".join(["starting chipseeker for", self.species, self.tissue]))
        with open(skr_log, "w") as skr_log_file:
            skr_res = subprocess.run(
                skr_cmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=skr_log_file,
                encoding="utf-8",
            )
        logger.info("log file: " + skr_log)
        if skr_res.returncode == 0:
            logger.info(" ".join(["chipseeker done for", self.species, self.tissue]))
        else:
            logger.error(" ".join(["chipseeker failed for", self.species, self.tissue]))

    @logger.catch
    def cord_conv(self):
        """
        使用liftOver将peak坐标转换到hg38基因组
        """

        anno_dir = os.path.join(self.store_dir, self.species, "anno/peaks")
        enhancer_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_enhancer"])
            + ".narrowPeak",
        )
        promoter_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_promoter"])
            + ".narrowPeak",
        )
        actProm_file = os.path.join(
            anno_dir,
            "_".join([self.species, self.tissue, "_overlap_optimal_actProm"])
            + ".narrowPeak",
        )
        nameList = pd.read_csv("/data1/Genome/Chain/list_all", sep="\t", header=None)
        abvName = (
            nameList[nameList[0] == self.species][1].values[0]
            if self.species != "Macaca_mulatta"
            else "macFas"
        )

        def name_awk(ele):
            awk_cmd = " ".join(
                ["awk -va=" + abvName + ".", "'{print a$0}'", ele, ">", ele + ".abv"]
            )
            logger.info(
                " ".join(
                    ["starting name modification for", self.species, self.tissue, ele]
                )
            )
            awk_res = subprocess.run(
                awk_cmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            logger.info(awk_res.stdout)
            if awk_res.returncode == 0:
                logger.info(
                    " ".join(
                        ["name modification done for", self.species, self.tissue, ele]
                    )
                )
            else:
                logger.error(
                    " ".join(
                        ["name modification failed for", self.species, self.tissue, ele]
                    )
                )

        def liftOver(ele):
            chain_file = (
                os.path.join(
                    "/data1/Genome/Chain", self.species + ".Homo_sapiens.chain"
                )
                if self.species != "Macaca_mulatta"
                else os.path.join(
                    "/data1/Genome/Chain", "Macaca_fascicularis.Homo_sapiens.chain"
                )
            )
            cord_cmd = " ".join(
                [
                    "liftOver -bedPlus=6",
                    ele + ".abv",
                    chain_file,
                    ele + ".hg38",
                    ele + ".unmapped",
                ]
            )
            logger.info(
                " ".join(["starting cord_conv for", self.species, self.tissue, ele])
            )
            cord_res = subprocess.run(
                cord_cmd,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
                encoding="utf-8",
            )
            logger.info(cord_res.stdout)
            if cord_res.returncode == 0:
                logger.info(
                    " ".join(["cord_conv done for", self.species, self.tissue, ele])
                )
            else:
                logger.error(
                    " ".join(["cord_conv failed for", self.species, self.tissue, ele])
                )

        for ele in [enhancer_file, promoter_file, actProm_file]:
            name_awk(ele)
            liftOver(ele)

    @logger.catch
    def close_logger(self):
        logger.info(
            "jobs for {0} {1} done, logger closed".format(self.species, self.tissue)
        )
        logger.remove(self.log_file)


# TODO: 增加每个步骤的完成标记，0为未完成，1为完成，表格中默认为0，主脚本调用各方法时，根据该标记决定是否执行。
