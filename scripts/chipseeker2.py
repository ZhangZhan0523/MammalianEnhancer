#!/usr/bin/env python
# coding=utf-8

# 列表文件名
# info_df = pd.read_csv("batch3fulana.csv")
# info_df = pd.read_csv(sys.argv[1])

import pandas as pd
import sys
import main


info_df = pd.read_csv(sys.argv[1])

for idx, row in info_df.iterrows():
    print(idx, row["species"], row["tissue"])
    # if idx > 1:
    TEMP = main.SpTi(
        species=row["species"],
        tissue=row["tissue"],
        dd=row["dd"],
        ac=row["ac"],
        me=row["me"],
        inp=row["inp"],
        ac_nm1=row["ac_nm1"],
        ac_nm2=row["ac_nm2"],
        me_nm1=row["me_nm1"],
        me_nm2=row["me_nm2"],
        inp_nm1=row["inp_nm1"],
        inp_nm2=row["inp_nm2"],
        ac_rep1_trimmed1=row["ac_rep1_trimmed1"],
        ac_rep1_trimmed2=row["ac_rep1_trimmed2"],
        ac_rep2_trimmed1=row["ac_rep2_trimmed1"],
        ac_rep2_trimmed2=row["ac_rep2_trimmed2"],
        me_rep2_trimmed1=row["me_rep2_trimmed1"],
        me_rep2_trimmed2=row["me_rep2_trimmed2"],
        inp_rep1_trimmed1=row["inp_rep1_trimmed1"],
        inp_rep1_trimmed2=row["inp_rep1_trimmed2"],
        inp_rep2_trimmed1=row["inp_rep2_trimmed1"],
        inp_rep2_trimmed2=row["inp_rep2_trimmed2"],
        genome_fa=row["genome_fa"],
        genome_tsv=row["genome_tsv"],
        genome_gff=row["genome_gff"],
        ana_dir=row["ana_dir"],
    )
    globals()[f'{row["species"]}_{row["tissue"]}'] = TEMP
    print(TEMP.chipseeker2())
    print(TEMP.close_logger())
