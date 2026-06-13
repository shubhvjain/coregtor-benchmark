#!/usr/bin/env python3
"""CoRegTor CLI."""

import argparse
from pathlib import Path
import sys
import os
from src.pipeline.core import CoRegTorPipeline
import subprocess
import sys
import argparse
import sys
import pandas as pd
import json
from dotenv import dotenv_values
from src.pipeline.run import reset_claimed
import src.pipeline.util as utp
import src.results.util as ut
from coregtor import RunCoRegTor
import joblib
import pickle

def run_bulk(extra_args):
    """
    access to snakemake workflow for bulk processing
    """
    snakefile = Path(__file__).parent / "Snakefile"
    
    # Base command
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--cores", "all",
        "--nolock"
    ]
    
    # Append the extra arguments (init, batch, --config, etc.)
    if extra_args:
        cmd.extend(extra_args)
        
    print(f"Executing: {' '.join(cmd)}") # Helpful for debugging
    subprocess.run(cmd, check=True)

def reset_claimed_genes(exp, config,worker):
    """
    """
    reset_claimed(exp, config,worker)

def run(env,input,args):
    """run small set of targets"""

    if input["type"] != "coregtor_run":
        raise ValueError("invalid exp file type.must be 'coregtor_run' ")

    rerun = args.rerun
    # folder setup 
    out_path, temp_path = ut.get_exp_path(input, env)
    out_path.mkdir(exist_ok=True, parents=True)

    ge_data = utp.read_dataset(input["dataset"],env)
    # print(ge_data)
    tf_list = utp.get_source_list(input["sources"].get("type","tf"),env)

    target_list = input.get("targets",[])
    # print(tf_list)
    if target_list is None :
        raise ValueError("no target genes provided")
    
    for target in target_list : 
        inter_file = out_path/ f"{target}_detail.pkl.z"
        if inter_file.exists() and not rerun:
            print(inter_file,"exists")
            continue
        print(f"starting {target}")
        res,intermediate = RunCoRegTor(ge_data,tf_list,target,        
            ge_sparsity_threshold=input.get("ge_sparsity_threshold",0.01),
            ensemble_method=input.get("ensemble_method","rf") ,
            ensemble_options= input.get("ensemble_options",None) ,
            distance_measure=input.get("distance_measure","canberra_distance"),
            cluster_method=input.get("cluster_method","community_detection"),
            cluster_options=input.get("cluster_options",None)
        )
        res.to_csv(out_path/f"{target}.csv",index=False)

        joblib.dump(intermediate, inter_file, compress=3)
        print(f"{target} done")


    # read ge data 
    # read source 
    # folder set up,
    # check if target is already run 



def main():
    parser = argparse.ArgumentParser(prog="coregtor")
    subparsers = parser.add_subparsers(dest="command")

    # The 'run' command for standard/simple use
    run_parser = subparsers.add_parser("run", help="Predict co-regulators for genes")
    run_parser.add_argument("--env",required=True)
    run_parser.add_argument("--input",required=True)
    run_parser.add_argument("--rerun", action="store_true", default=False,
                        help="Force rerun even if output already exists.")


    # The 'bulk' command explicitly for Snakemake
    bulk_parser = subparsers.add_parser("bulk", help="Run optimized Snakemake pipeline")
    
    # Use REMAINDER to catch everything after 'bulk' (init, batch, etc.)
    bulk_parser.add_argument("snakemake_args", nargs=argparse.REMAINDER, 
                             help="Subcommands and arguments to pass to Snakemake")


    # The 'reset_claimed' command explicitly for Snakemake
    reset_claimed_parser = subparsers.add_parser("reset_claimed", help="Reset exp")
    
    # Use REMAINDER to catch everything after 'bulk' (init, batch, etc.)
    reset_claimed_parser.add_argument("--env",required=True)
    reset_claimed_parser.add_argument("--input",required=True)
    reset_claimed_parser.add_argument("--worker",default=None)


    test_parser = subparsers.add_parser("hi", help="Just a test")


    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    
    if args.command == "run":
        env = ut.get_env(args.env)
        input = ut.get_input(args.input)
        run(env,input,args)    

    elif args.command == "bulk":
        run_bulk(args.snakemake_args)

    elif args.command == "reset_claimed":
        env = ut.get_env(args.env)
        input = ut.get_input(args.input)
        worker= args.worker
        reset_claimed_genes(input,env,worker)

    elif args.command == "hi":
        print("hi")
        print("bye")



if __name__ == "__main__":
    main()