#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import multiprocessing as mp
from loguru import logger
import subprocess
import sys
from Bio import SeqIO
from functools import partial
from blast_all_mcl import make_similarity_matrix

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
def main():
    make_similarity_matrix()


if __name__ == "__main__":
    main()
