import argparse
import json
import sys
from pathlib import Path
from dotenv import dotenv_values
import src.results.util as ut
from joblib import Parallel, delayed, load
import pandas as pd

from src.results.explore import gene_frequency_analysis, distance_measure_analysis
from src.results.indices import performance_indices1,performance_indices_combined,performance_indices_ppi_network,performance_indices_dcorr
from src.results.target import target_info,target_run_diagram
from src.results.parameter_selection import parameter_selection_plot
from src.results.frequent_coregs import find_frequent_coregs
from src.results.clustering import identify_coregulators

def hello(input, env, options, args):
    target = input.get("target", "world")
    print(f"Hello, {target}!")
    print(f"  env keys : {list(env.keys())}")
    print(f"  verbose  : {args.verbose}")


def get_target_run(input, env):
    ""


def _get_temp_file(temp_path, target):
    """"""
    file_path = temp_path/"results"/f"{target}.pkl"
    data = load(file_path)
    return data

def get_clustering_results(input,target,cluster_option_name,temp_file,rerun):
    """get clusters for a given config, generate if already not generated"""
    if cluster_option_name in temp_file["results"]["clusters"].keys() and not rerun:
        return temp_file["results"]["clusters"][cluster_option_name]
    else:
        coptions = input["clustering"]
        c = next((d for d in coptions if d.get("id") == cluster_option_name), None)
        if c is  None:
            raise ValueError("Invalid cluster options name")
        id = c["id"]
        dist_matrix = temp_file["results"]["distance"][c["distance"]]
        #print(dist_matrix)
        c_method = c["method"]
        c_options = c["options"]
        c_notes = c["note"]
        res = identify_coregulators(dist_matrix, target, c_method, c_options, c_notes)
        return res


def _process_single_target(t, temp_path, input, result_list, rerun):
    """Worker function to process all cluster configurations for a single target."""
    data = _get_temp_file(temp_path, t)
    target_results = {}
    for r in result_list:
        target_cluster_df = get_clustering_results(input, t, r, data, rerun)
        target_results[r] = target_cluster_df
    return target_results


def generate_result_file(input, env, options, args):
    """
    generate results files
    """
    out_path, temp_path = ut.get_exp_path(input, env)

    n_jobs = args.njobs if args.njobs is not None else -1
    rerun = args.rerun

    run_done = temp_path / "run_success.txt"

    # if  not run_done.exists():
    #     print("Experiment not run successfully")
    #     return

    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)

    # with open(out_path/"input.json") as f:
    #     input_data = json.load(f)

    success_target = ut.get_exp_target_list(out_path)
    #print(len(success_target))

    all_results = [cl["id"] for cl in input["clustering"]]
    result_list_input = options.get("cluster_list", all_results)

    result_list = []  # only that are not already done
    for r in result_list_input:
        result_file_path = cluster_folder / f"{r}.csv"
        if result_file_path.exists() and not rerun:
            continue
        result_list.append(r)

    # # simple non parallel code
    # result_objs = {r: [] for r in result_list}
    # for t in success_target:
    #     data = _get_temp_file(temp_path, t)
    #     for r in result_list:
    #         target_cluster_df = get_clustering_results(input,t,r,data,rerun)
    #         #print(target_cluster_df)
    #         result_objs[r].append(target_cluster_df)

    # --- PARALLEL EXECUTION WITH JOBLIB ---
    # Run the worker function across multiple cores
    parallel_outputs = Parallel(n_jobs=n_jobs)(
        delayed(_process_single_target)(t, temp_path, input, result_list, rerun)
        for t in success_target
    )
    # Reconstruct the result_objs dict from the parallel workers' outputs
    result_objs = {r: [] for r in result_list}
    for target_res in parallel_outputs:
        for r in result_list:
            result_objs[r].append(target_res[r])
    # ---------------------------------------

    master_file = []
    for r in result_list:
        df = pd.concat(result_objs[r], ignore_index=True)
        master_file.append(df)
        result_file_path = cluster_folder / f"{r}.csv"
        df.to_csv(result_file_path, index=False)

    # this creates a folder with the name and a file inside it titled clusters.csv.gz
    save_to_output = options.get("save_to_output",True)    
    if save_to_output:
        result_folder_out = out_path/f"{options.get('result_name','default')}"
        result_folder_out.mkdir(exist_ok=True, parents=True)
        final_df = pd.concat(master_file, ignore_index=True)
        final_df.to_csv(result_folder_out/"clusters.csv.gz", index=False)



METHOD_REGISTRY = {
    "hello": hello,
    "generate_result_file": generate_result_file,
    "gene_frequency_analysis": gene_frequency_analysis,
    "performance_indices1": performance_indices1,
    "distance_measure_analysis":distance_measure_analysis,
    "performance_indices_combined":performance_indices_combined,
    "target_info":target_info,
    "target_run_diagram":target_run_diagram,
    "parameter_selection_plot":parameter_selection_plot,
    "performance_indices_ppi_network":performance_indices_ppi_network,
    "find_frequent_coregs":find_frequent_coregs,
    "performance_indices_dcorr":performance_indices_dcorr
}

# CLI entry point

def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate experiment results",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--env",   required=True,
                        metavar="PATH", help="Path to the .env file.")
    parser.add_argument("--input", required=True,
                        metavar="PATH", help="Path to the JSON input file.")
    parser.add_argument(
        "--id",
        required=True,
        metavar="ID",
        help="ID of the result object in input['results'] to run.",
    )

    # Optional arguments — add more as needed
    parser.add_argument("--njobs", type=int, default=-1,
                        help="Number of parallel jobs. Defaults to -1 (all cores).")
    parser.add_argument("--rerun", action="store_true", default=False,
                        help="Force rerun even if output already exists.")
    parser.add_argument("--verbose", action="store_true",
                        default=False, help="Enable verbose output.")
    parser.add_argument("--dry-run", action="store_true",
                        default=False, help="Load files but skip execution.")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    env = ut.get_env(args.env)
    input_data = ut.get_input(args.input)

    options = ut.get_result_by_id(input_data, args.id)

    result_type = options.get("type")
    if result_type not in METHOD_REGISTRY:
        print(
            f"[error] unknown type '{result_type}'. Registered types: {', '.join(METHOD_REGISTRY.keys())}", file=sys.stderr)
        sys.exit(1)

    fn = METHOD_REGISTRY[result_type]

    if args.dry_run:
        print(
            f"[dry-run] would call '{result_type}' for id '{args.id}' with {len(env)} env key(s).")
        return

    fn(input_data, env, options, args)


if __name__ == "__main__":
    main()
