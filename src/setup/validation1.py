"""
This is to setup ppi network analysis 
"""
import src.results.util as ut
import argparse
from pathlib import Path
import os
from tfitpy.datasets.trap_cache import generate_and_cache_trap_scores

from tfitpy.datasets.trap_cache import generate_promoter_reference, get_promoter_reference,_encode_sequence_numba, _execute_trap_matrix_numba,list_human_tfs_pyjaspar


import numpy as np
import pandas as pd
from pathlib import Path
from numba import njit, prange

def generate_files(env,args):
    """"""
    data_path = Path( os.path.expandvars( env.get("DATA_PATH")))
    # generate_promoter_reference(data_path)
    generate_and_cache_trap_scores(data_path)



def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate setup files for ppi network analysis",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--env",   required=True,
                        metavar="PATH", help="Path to the .env file.")
    

    # Optional arguments — add more as needed
    parser.add_argument("--njobs", type=int, default=-1,
                        help="Number of parallel jobs. Defaults to -1 (all cores).")
    parser.add_argument("--rerun", action="store_true", default=False,
                        help="Force rerun even if output already exists.")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    env = ut.get_env(args.env)
    generate_files(env,args)


if __name__ == "__main__":
    main()