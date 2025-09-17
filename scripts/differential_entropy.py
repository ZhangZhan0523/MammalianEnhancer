#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os
from scipy.stats import differential_entropy, norm
from loguru import logger

"""
This script is used to calculate the differential entropy of a given dataset.
the dataset should be peak info of each species, and the columns should be the same.
"""


@logger.catch
def main():
    lib_factor = pd.read_csv(
        "/media/Data/zhangz/chip/scripts2/info/fc_library_factor.csv", header=0
    )
    # load the data
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
        "Rhinopithecus_roxellana",
        "Tupaia_belangeri",
        "Procavia_capensis",
        "Petaurus_breviceps",
    ]
    ent = pd.DataFrame()
    for sp in species:
        info = pd.read_csv(
            f"/media/Data/zhangz/chip/analysis/summary2/abc_1M/{sp}_overlap_anno.csv",
            header=0,
            usecols=["unique_id", "peak", "species", "tissue", "element", "lib_fc"],
        )
        ## select rep_factor from lib_factor dataframe according to species and tissue of
        ## each line, calculate raw_p as lib_fc/rep_factor
        info["raw_p"] = info.apply(
            lambda row: row["lib_fc"]
            / lib_factor.loc[
                (lib_factor["species"] == row["species"])
                & (lib_factor["tissue"] == row["tissue"]),
                "rep_factor",
            ].values[0],
            axis=1,
        )
        # calculate differential entropy for raw_p in each tissue and enhancer\promoter
        # ent_info = info.groupby(['species', 'tissue'])['raw_p'].apply(lambda x: differential_entropy(x)).reset_index()
        #     ent_info = (
        #         info.groupby(["species", "tissue"])["lib_fc"]
        #         .apply(lambda x: differential_entropy(x))
        #         .reset_index()
        #     )
        #     ent_info = ent_info.pivot(index="species", columns="tissue", values="lib_fc")
        #     ent_info = ent_info.reset_index()
        #     ent_info2 = (
        #         info.groupby(["species", "tissue", "element"])["lib_fc"]
        #         .apply(lambda x: differential_entropy(x))
        #         .reset_index()
        #     )
        #     ent_info2["te"] = ent_info2["tissue"] + "_" + ent_info2["element"]
        #     ent_info2 = ent_info2.drop(["tissue", "element"], axis=1)
        #     ent_info2 = ent_info2.pivot(index="species", columns=["te"], values="lib_fc")
        #     ent_info2 = ent_info2.reset_index()
        #     ent_info = pd.merge(ent_info, ent_info2, on="species", how="left")
        #     ent = pd.concat([ent, ent_info], axis=0)
        # ent.to_csv(
        #     "/media/Data/zhangz/chip/analysis/summary2/tissue_ent/differential_entropy.csv",
        #     index=False,
        # )
        ent_info = (
            info.groupby(["species", "tissue"])["raw_p"]
            .apply(lambda x: differential_entropy(x))
            .reset_index()
        )
        ent_info = ent_info.pivot(index="species", columns="tissue", values="raw_p")
        ent_info = ent_info.reset_index()
        ent_info2 = (
            info.groupby(["species", "tissue", "element"])["raw_p"]
            .apply(lambda x: differential_entropy(x))
            .reset_index()
        )
        ent_info2["te"] = ent_info2["tissue"] + "_" + ent_info2["element"]
        ent_info2 = ent_info2.drop(["tissue", "element"], axis=1)
        ent_info2 = ent_info2.pivot(index="species", columns=["te"], values="raw_p")
        ent_info2 = ent_info2.reset_index()
        ent_info = pd.merge(ent_info, ent_info2, on="species", how="left")
        ent = pd.concat([ent, ent_info], axis=0)
    ent.to_csv(
        "/media/Data/zhangz/chip/analysis/summary2/tissue_ent/differential_entropy_raw_p.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
