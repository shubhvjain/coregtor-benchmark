from tfitpy import compute_indices
import json
import os
from pathlib import Path
import pandas as pd
import  src.results.util as ut 


def compute_indices_core(data, bio_data_path):
    """Core computation - single dataframe in, scored dataframe out."""
    if data is None or data.empty:
        raise ValueError("No data provided")

    # Filter to clusters with multiple sources
    data = data[data['sources'].str.split(";").str.len() > 1].copy()

    n_job_count = 1
    
    df_score, add_data = compute_indices(
        df=data, data_path=bio_data_path, new_methods_only=True, n_jobs=n_job_count)
    return df_score, add_data


def compute_indices_batch(data, parallel=True, n_jobs=4, batch_size=1000,bio_data_path=None):
    """
    Compute indices with batching and optional parallelization.

    Args:
        data: Input dataframe
        bio_data_path: Path to bio data (defaults to $DATA_PATH env var)
        parallel: Whether to process batches in parallel
        n_jobs: Number of parallel workers (if parallel=True, -1 for all cores)
        batch_size: Rows per batch (if None, processes all at once)

    Returns:
        df_score: Scored dataframe
        add_data: Additional metadata
    """

    # bio_data_path = Path(os.path.expandvars(os.getenv("DATA_PATH")))
    #print(bio_data_path)
    # Split into batches

    if parallel:
        from joblib import Parallel, delayed
        batches = [data.iloc[i:i+batch_size]
                   for i in range(0, len(data), batch_size)]
        results = Parallel(n_jobs=n_jobs)(
            delayed(compute_indices_core)(batch, bio_data_path)
            for batch in batches
        )
        # Combine results
        df_scores = [r[0] for r in results]
        combined_df = pd.concat(df_scores, ignore_index=True)

        combined_add_data = [r[1] for r in results]
        return combined_df, combined_add_data

    else:
        return compute_indices_core(data, bio_data_path)


def performance_indices1(input, env,options ,args):
    """
    compute performance indices for cluster results 
    """
    out_path, temp_path = ut.get_exp_path(input,env)
    rerun = args.rerun
    n_jobs = args.njobs if args.njobs is not None else 4
    print(n_jobs)
    bio_path = env["DATA_PATH"]

    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)

    all_results = [cl["id"] for cl in input["clustering"]]
    result_list_input = options.get("result_list", all_results)

    result_list = []  # only that are not already done
    for r in result_list_input:
        result_file_path = cluster_folder / f"index_{r}.csv"
        if result_file_path.exists() and not rerun:
            continue
        result_list.append(r)


    for r in result_list:
        print(r)
        cluster_file = cluster_folder / f"{r}.csv"
        index_file = cluster_folder / f"index_{r}.csv"
        data = pd.read_csv(cluster_file)
        data.rename(columns={"uid":"cluster_uid"}, inplace=True)

        run_parallel = len(data) > 10
        batch_size = int(len(data)/n_jobs) + 1
        print()
        df, add = compute_indices_batch(
            data, run_parallel, n_jobs, batch_size,bio_path)
        df.to_csv(index_file,index=False)
        print("saved")



def performance_indices_combined(input, env,options ,args):
    """
    compute performance indices for cluster results 
    """
    out_path, temp_path = ut.get_exp_path(input,env)
    rerun = args.rerun
    cluster_folder = temp_path/"clusters"

    dataset_name  = input["dataset"]["name"]
    
    file_name = options.get("name","default")

    result_file =  out_path/ f"{file_name}_indices.csv.gz"
    if result_file.exists() and not rerun:
        print("already exists")
        return

    all_results = [cl["id"] for cl in input["clustering"]]
    result_list = options.get("result_list", all_results)

    res = []

    for r in result_list:
        print(r)
        index_file = cluster_folder / f"index_{r}.csv"
        if not index_file.exists():
            raise ValueError(f"file {index_file} note yet generated")
        
        data = pd.read_csv(index_file)
        data["config_name"] = r
        res.append(data)

    combined_df = pd.concat(res, ignore_index=True)
    combined_df["dataset"] = dataset_name
    combined_df.to_csv(result_file,index=False)
    print("saved")