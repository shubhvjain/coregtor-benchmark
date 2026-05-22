from tfitpy import compute_indices, load_all_network_data
import json
import os
from pathlib import Path
import pandas as pd
import  src.results.util as ut 
import numpy as np
from joblib import Parallel, delayed

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
    # print(n_jobs)
    bio_path = env["DATA_PATH"]

    # first get the list of cluster files in the result name 
    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))

    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)

    all_results = [cl["id"] for cl in input["clustering"]]
    
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list
    
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
        # data.rename(columns={"uid":"cluster_uid"}, inplace=True)

        run_parallel = len(data) > 10
        batch_size = int(len(data)/n_jobs) + 1
        df, add = compute_indices_batch(
            data, run_parallel, n_jobs, batch_size,bio_path)
        df.to_csv(index_file,index=False)
        print("saved")


def compute_module_network_metrics(M,t, ppi, null_graphs, mapping, debug=False):
    """
    M: List of gene symbols
    ppi: Original aligned igraph object
    null_graphs: List of  pre-built igraph objects
    mapping: Symbol -> Int dictionary
    debug: If True, prints mapping and subgraph samples
    """
    # 1. Map symbols to indices
    m_indices = [mapping[g] for g in M if g in mapping]
    n_m = len(m_indices)
    target_idx = mapping.get(t)
        
    if n_m < 2:
        # lcc is 1/1
        return 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, round(n_m/len(M),4)

    possible_edges = (n_m * (n_m - 1)) / 2

    # --- 2. Observed Metrics ---
    sub_obs = ppi.induced_subgraph(m_indices)
    e_obs = sub_obs.ecount()
    
    obs_density = e_obs / possible_edges if possible_edges > 0 else 0.0
    obs_lcc_count = sub_obs.connected_components().giant().vcount() if e_obs > 0 else 1
    obs_lcc_ratio = obs_lcc_count / n_m

    obs_tc_count = 0
    if target_idx is not None:
        # Get neighbors of target and find intersection with module indices
        target_neighbors = set(ppi.neighbors(target_idx))
        obs_tc_count = len(target_neighbors.intersection(m_indices))
    
    obs_tc_ratio = obs_tc_count / n_m

    # --- 3. Null Distribution ---
    density_hits = 0
    lcc_hits = 0
    tc_hits = 0
    R = len(null_graphs)
    
    for i, g_null in enumerate(null_graphs):
        sub_null = g_null.induced_subgraph(m_indices)
        e_null = sub_null.ecount()
        
        # Null Density
        null_density = e_null / possible_edges if possible_edges > 0 else 0.0
        if null_density >= obs_density:
            density_hits += 1
            
        # Null LCC Ratio
        null_lcc_count = sub_null.connected_components().giant().vcount() if e_null > 0 else 1
        null_lcc_ratio = null_lcc_count / n_m
        
        if null_lcc_ratio >= obs_lcc_ratio:
            lcc_hits += 1

        if target_idx is not None:
            null_target_neighbors = set(g_null.neighbors(target_idx))
            null_tc_count = len(null_target_neighbors.intersection(m_indices))
            if (null_tc_count / n_m) >= obs_tc_ratio:
                tc_hits += 1
            

    # --- 4. Final Scores (-ln(p)) ---
    p_density = (density_hits + 1) / (R + 1)
    density_score = max(0.0, -np.log(p_density))
    
    p_lcc = (lcc_hits + 1) / (R + 1)
    lcc_score = max(0.0, -np.log(p_lcc))

    p_tc = 0
    target_connectivity_score = 0
    if target_idx is not None:
        p_tc = (tc_hits + 1) / (R + 1)
        target_connectivity_score = max(0.0,-np.log(p_tc))        


    nodes_found_ratio = n_m/ len(M) 

    return obs_density, density_score, obs_lcc_ratio, lcc_score, obs_tc_ratio, target_connectivity_score, nodes_found_ratio


def compute_network_scores(df, ppi_key, data_path, n_models=500):
    """
    Computes PPI3, PPI4, and PPI5 for a single PPI source.
    Returns a DataFrame containing only the new metrics, indexed identically to input df.
    """
    # 1. Load the specific PPI data
    G, NULL_MODELS, MAPPING = load_all_network_data(data_path, ppi_key, n_models)
    
    # 2. Compute results row by row
    zero_result = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0)
    results = []
    print(f"starting {len(df)} rows")
    
    for _, row in df.iterrows():
        sources = row['sources'].split(';')
        target = row['target']
        
        # Check Condition 1: Module size must be at least 3
        if len(sources) < 2:
            results.append(zero_result)
            continue
            
        
        metrics = compute_module_network_metrics(
            M=row['sources'].split(';'),
            t=row['target'],
            ppi=G,
            null_graphs=NULL_MODELS,
            mapping=MAPPING
        )
        results.append(metrics)
    
    # 3. Create a result DataFrame with specific column names
    col_names = [
        f"density_{ppi_key}", f"density_score_{ppi_key}", 
        f"lcc_{ppi_key}", f"lcc_score_{ppi_key}",
        f"tc_{ppi_key}", f"tc_score_{ppi_key}",f"node_found_ratio_{ppi_key}" 
    ]
    
    result_df = pd.DataFrame(results, columns=col_names, index=df.index)
    return result_df.round(5)

def compute_scores_preloaded(df, ppi_key, G, NULL_MODELS, MAPPING):
    """Modified version of compute_network_scores that takes pre-loaded data"""
    results = []
    zero_result = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    
    for _, row in df.iterrows():
        sources = str(row['sources']).split(';')
        if len(sources) < 2:
            results.append(zero_result)
            continue
            
        metrics = compute_module_network_metrics(
            M=sources,
            t=row['target'],
            ppi=G,
            null_graphs=NULL_MODELS,
            mapping=MAPPING
        )
        results.append(metrics)
    
    col_names = [
        f"density_{ppi_key}", f"density_score_{ppi_key}", 
        f"lcc_{ppi_key}", f"lcc_score_{ppi_key}",
        f"tc_{ppi_key}", f"tc_score_{ppi_key}", f"node_found_ratio_{ppi_key}" 
    ]
    return pd.DataFrame(results, columns=col_names, index=df.index).round(5)

def performance_indices_ppi_network1(input, env,options ,args):
    """
    compute performance indices for cluster results 
    """
    out_path, temp_path = ut.get_exp_path(input,env)
    rerun = args.rerun
    n_jobs = args.njobs if args.njobs is not None else 4
    # print(n_jobs)
    bio_path = env["DATA_PATH"]

    # first get the list of cluster files in the result name 
    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))

    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)

    all_results = [cl["id"] for cl in input["clustering"]]
    
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list

    # cluster_folder = temp_path/"clusters"
    # cluster_folder.mkdir(exist_ok=True, parents=True)

    # all_results = [cl["id"] for cl in input["clustering"]]
    # result_list_input = options.get("result_list", all_results)

    result_list = []  # only that are not already done
    for r in result_list_input:
        result_file_path = cluster_folder / f"index_network_{r}.csv"
        if result_file_path.exists() and not rerun:
            continue
        result_list.append(r)


    for r in result_list:
        print(r)
        cluster_file = cluster_folder / f"{r}.csv"
        index_file = cluster_folder / f"index_network_{r}.csv"
        data = pd.read_csv(cluster_file)
        #data.rename(columns={"uid":"cluster_uid"}, inplace=True)

        ppi_list = ['hippie', 'stringdb', 'biogrid']

        # This returns a list of 3 DataFrames
        # list_of_result_dfs = Parallel(n_jobs=3)(
        #     delayed(compute_network_scores)(data, ppi, bio_path,500) for ppi in ppi_list
        # )

        list_of_result_dfs = []
        for ppi in ppi_list:
            print(f"starting {ppi}.")
            r = compute_network_scores(data, ppi, bio_path,500)
            list_of_result_dfs.append(r)
            print(f"Finished {ppi}.")

        # Combine horizontally: result will have original df cols + 18 new PPI cols
        # axis=1 tells pandas to align by the index (cluster_uid)
        final_df = pd.concat([data] + list_of_result_dfs, axis=1)
        # result = compute_network_scores(data,"stringdb",bio_path,10)
        final_df.to_csv(index_file,index=False)


        # run_parallel = len(data) > 10
        # batch_size = int(len(data)/n_jobs) + 1
        # print()
        # df, add = compute_indices_batch(
        #     data, run_parallel, n_jobs, batch_size,bio_path)
        # df.to_csv(index_file,index=False)
        # print("saved")


def performance_indices_ppi_network(input, env, options, args):
    out_path, temp_path = ut.get_exp_path(input, env)
    bio_path = env["DATA_PATH"]
    rerun = args.rerun
    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))

    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)

    all_results = [cl["id"] for cl in input["clustering"]]
    
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list
        
    ppi_list = ['hippie', 'stringdb', 'biogrid']
    
    # 1. Initialize a dictionary to hold DataFrames for each result ID
    # This prevents constant disk I/O for the cluster files
    df_store = {}
    for r in result_list_input:
        result_file_path = cluster_folder / f"index_network_{r}.csv"
        if result_file_path.exists() and not rerun:
            continue
        cluster_file = cluster_folder / f"{r}.csv"
        df_store[r] = pd.read_csv(cluster_file)
        #df_store[r].rename(columns={"uid": "cluster_uid"}, inplace=True)

    # 2. Outer Loop: PPI Networks (The "Heavy" objects)
    for ppi_key in ppi_list:
        print(f"--- Loading Network Data for: {ppi_key} ---")
        # Load the 500 models ONCE for this PPI
        G, NULL_MODELS, MAPPING = load_all_network_data(bio_path, ppi_key, 200)
        print("now running scores")
        # Inner Loop: Process all cluster files against this loaded network
        for r, data_df in df_store.items():
            print(f"Computing {ppi_key} scores for cluster set: {r}")
            
            # We call the computation logic directly here or via a modified helper
            # that accepts the pre-loaded objects
            scores_df = compute_scores_preloaded(data_df, ppi_key, G, NULL_MODELS, MAPPING)
            
            # Join the new scores to the existing DataFrame in our dictionary
            df_store[r] = pd.concat([df_store[r], scores_df], axis=1)
        
        # Optional: Clear memory of the large graph objects before next PPI
        del G, NULL_MODELS, MAPPING

    # 3. Final Step: Save consolidated results
    for r, final_df in df_store.items():
        index_file = cluster_folder / f"index_network_{r}.csv"
        final_df.to_csv(index_file, index=False)
        print(f"Saved consolidated results to {index_file}")


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
        
        network_file = cluster_folder / f"index_network_{r}.csv"
        if not network_file.exists():
            raise ValueError(f"file {network_file} note yet generated")
        
        data = pd.read_csv(index_file)

        network_score = pd.read_csv(network_file)
        network_cols = [
            "density_hippie", "density_score_hippie", "lcc_hippie", "lcc_score_hippie", 
            "tc_hippie", "tc_score_hippie", "node_found_ratio_hippie", "density_stringdb", 
            "density_score_stringdb", "lcc_stringdb", "lcc_score_stringdb", "tc_stringdb", 
            "tc_score_stringdb", "node_found_ratio_stringdb", "density_biogrid", 
            "density_score_biogrid", "lcc_biogrid", "lcc_score_biogrid", "tc_biogrid", 
            "tc_score_biogrid", "node_found_ratio_biogrid"
        ]

        # Merge the network columns into 'data'
        # If network_score has 1 row and data has many, this broadcasts the values
        # If they align row-for-row, this joins them side-by-side
        data = pd.concat([data, network_score[network_cols]], axis=1)

        data["config_name"] = r
        res.append(data)

    combined_df = pd.concat(res, ignore_index=True)
    combined_df["dataset"] = dataset_name
    combined_df.to_csv(result_file,index=False)
    print("saved")