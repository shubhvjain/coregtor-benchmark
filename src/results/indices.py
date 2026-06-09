from tfitpy import compute_indices, load_all_network_data
import json
import os
from pathlib import Path
import pandas as pd
import  src.results.util as ut 
import src.pipeline.util as put
import numpy as np
from joblib import Parallel, delayed
from src.results.dcorr import compute_module_scores,generate_combined_dcor_cache, compute_common_targets_EC_scores
import pyarrow as pa
import pyarrow.parquet as pq
from src.results.frequent_index import compute_GO_ORA,compute_PPI_shared_partners

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


def compute_dcorr_scores(
    data_df: pd.DataFrame, 
    ge_data: pd.DataFrame, 
    src: list, 
    tgt: list, 
    flat_cache, 
    gene_map, 
    total_genes,
    n_jobs=-1):
    """
    Computes distance correlation metrics (dCor_source, dCor_target, DC1, DC2) 
    parallelly across all modules defined in data_df.
    
    Parameters:
    -----------
    data_df : pd.DataFrame
        DataFrame containing cluster/module definitions. Must contain columns
        representing the target gene and the list of module member genes.
        Expected columns: 'target' and 'genes' (where 'genes' is a comma-separated string or list).
    ge_data : pd.DataFrame
        The raw gene expression DataFrame.
    src : list of str
        The source gene list.
    tgt : list of str
        The target gene list.
    ss_cache : pd.DataFrame
        Precomputed source-to-source distance correlation matrix.
    st_cache : pd.DataFrame
        Precomputed source-to-target distance correlation matrix.
    n_jobs : int
        Number of worker threads/processes to use.
        
    Returns:
    --------
    scores_df : pd.DataFrame
        DataFrame with columns ['dCor_source', 'dCor_target', 'DC1', 'DC2']
        aligned precisely with the indices of data_df.
    """
    print(f"Parallelizing module scoring over {len(data_df)} items using {n_jobs} cores...")
    
    # Standardize data_df input columns
    # Assumes 'target' contains the target gene and 'genes' contains the module members
    targets = data_df['target'].tolist()
    
    # Handle potentially stringified lists from CSV format ("geneA,geneB,geneC")
    module_lists = []
    for item in data_df['sources']:
        if isinstance(item, str):
            module_lists.append([g.strip() for g in item.split(';')])
        else:
            module_lists.append(list(item))
            
    # Execute scoring functions in parallel
    # Using default backend ('loky') to safely handle process isolation and numpy memmapping
    results = Parallel(n_jobs=n_jobs)(
        delayed(compute_module_scores)(
            module_genes=mod, 
            target_gene=t, 
            source_gene_set = src,
            flat_cache=flat_cache, 
            gene_to_idx=gene_map, 
            N=total_genes
        ) for mod, t in zip(module_lists, targets)
    )
    
    # Reconstruct back into a clean structured DataFrame matching the input index
    scores_df = pd.DataFrame(results, index=data_df.index)
    return scores_df

def get_dataset_dcorr_cache(data_path,dataset):
    """
    """
    if dataset.get("type",None) != "gct" :
        raise ValueError("not found")
    
    file1 =  data_path /'dcorr_cache' /f"{dataset['name']}_pc.parquet"
    file2 = data_path /'dcorr_cache' /f"{dataset['name']}_metadata.json"

    with open(file2, 'r') as f:
        metadata = json.load(f)
        
    num_genes = metadata["num_genes"]
    gene_to_idx = metadata["gene_to_idx"]

    table = pq.read_table(file1)
    # Convert the PyArrow column back to a standard flat 1D NumPy array
    flat_dcor_cache = table.column('dcor').to_numpy()
    
    print(f"Cache loaded successfully. Array shape: {flat_dcor_cache.shape}")
    return flat_dcor_cache, gene_to_idx, num_genes



def performance_indices_dcorr(input, env, options, args):
    out_path, temp_path = ut.get_exp_path(input, env)
    bio_path = Path(os.path.expandvars(env["DATA_PATH"]))
    rerun = args.rerun
    n_jobs = args.njobs if args.njobs is not None else 4
    
    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))
    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)
    all_results = [cl["id"] for cl in input["clustering"]]
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list

    ge_data = put.read_dataset(input["dataset"],env)
    src,tgt = ut.get_input_source_list(out_path)
    flat_cache, gene_map, total_genes = get_dataset_dcorr_cache(bio_path,input["dataset"])

    # 1. Initialize a dictionary to hold DataFrames for each result ID
    # This prevents constant disk I/O for the cluster files
    df_store = {}
    for r in result_list_input:
        result_file_path = cluster_folder / f"index_dcorr_{r}.csv"
        if result_file_path.exists() and not rerun:
            continue
        cluster_file = cluster_folder / f"{r}.csv"
        df_store[r] = pd.read_csv(cluster_file)
    
    for r, data_df in df_store.items():
        print(f"Computing dcor scores for cluster set: {r}")
        scores_df = compute_dcorr_scores(data_df,ge_data,src,tgt,flat_cache, gene_map, total_genes)
        df_store[r] = pd.concat([df_store[r], scores_df], axis=1)

    # Final Step: Save consolidated results
    for r, final_df in df_store.items():
        index_file = cluster_folder / f"index_dcorr_{r}.csv"
        final_df.to_csv(index_file, index=False)
        print(f"Saved consolidated results to {index_file}")


def compute_tfbs_score(target_gene, gene_list, cache_df):
    """
    Computes the pure TF found percentage and average TFBS affinity 
    for a single target gene and a discrete list of input genes.
    """
    target_clean = str(target_gene).strip().upper()
    input_set = {str(g).strip().upper() for g in gene_list if pd.notna(g)}
    total_input_count = len(input_set)

    if total_input_count == 0:
        return np.nan

    cache_columns = set(cache_df.columns)
    matched_tfs = input_set.intersection(cache_columns)
    n_matched = len(matched_tfs)

    if target_clean not in cache_df.index or n_matched == 0:
        return np.nan 

    observed_values = cache_df.loc[target_clean, list(matched_tfs)]
    mean_affinity = np.nanmean(observed_values)

    return mean_affinity


def _process_tfbs_chunk(
    chunk_df,
    cache_df,
    all_gene_sources,
    n_permutations,
    seed
):
    """
    Processes a sub-batch of data_df rows. Evaluates the observed TFBS affinity,
    TF found percentage, and runs a gene-set level permutation test.
    """
    rng = np.random.default_rng(seed)
    global_pool = np.array(list(all_gene_sources))
    cache_columns = set(cache_df.columns)
    results = []

    for _, row in chunk_df.iterrows():
        target = str(row["target"]).strip().upper()
        sources_str = row["sources"]

        if pd.isna(sources_str) or not str(sources_str).strip():
            results.append({
                "TFBS_affinity": np.nan,
                "TF_found_per": 0.0,
                "TFBS_affinity_score": np.nan
            })
            continue

        # Parse unique items from this row's source configuration
        row_sources = list({g.strip().upper() for g in str(sources_str).split(";") if g.strip()})
        total_source_count = len(row_sources)

        # Pre-evaluate intersection count for accurate TF_found percentage
        matched_tfs = [g for g in row_sources if g in cache_columns]
        n_matched = len(matched_tfs)
        
        tf_found_pct = n_matched / total_source_count if total_source_count > 0 else 0.0
        
        # Call your original function safely
        observed_mean = compute_tfbs_score(target, row_sources, cache_df)

        if pd.isna(observed_mean) or n_matched == 0:
            results.append({
                "TFBS_affinity": np.nan,
                "TF_found_per": tf_found_pct,
                "TFBS_affinity_score": np.nan
            })
            continue

        # Permutation Test over All Unique Gene Sources
        extreme_count = 0
        for _ in range(n_permutations):
            simulated_genes = rng.choice(global_pool, size=n_matched, replace=False)
            simulated_mean = compute_tfbs_score(target, simulated_genes, cache_df)
            
            if pd.notna(simulated_mean) and simulated_mean >= observed_mean:
                extreme_count += 1

        p_value = (extreme_count + 1) / (n_permutations + 1)
        score = max(0.0, -np.log(p_value))
        results.append({
            "TFBS_affinity": round(observed_mean,5),
            "TF_found_per": round(tf_found_pct,2),
            "TFBS_affinity_score": round(score,5)
        })

    return pd.DataFrame(results)


def compute_tfbs_scores_df(data_df, cache, all_gene_sources, n_permutations=200, n_jobs=-1, seed=42):
    """
    Parallel computation of TFBS affinity metrics where TF_found is calculated 
    out of the parsed 'sources' row elements.
    """
    # 1. Standardize cache labels ONCE globally before passing to workers
    cache_clean = cache.copy()
    cache_clean.index = cache_clean.index.astype(str).str.strip().str.upper()
    cache_clean.columns = cache_clean.columns.astype(str).str.strip().str.upper()

    # Pre-clean the global pool string casing
    cleaned_global_sources = {str(g).strip().upper() for g in all_gene_sources if pd.notna(g)}

    # 2. Chunk input dataframe based on available workers
    n_workers = n_jobs if n_jobs > 0 else 4
    chunks = np.array_split(data_df, n_workers)

    # 3. Leverage Joblib to process batches independently
    processed_dfs = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_process_tfbs_chunk)(
            chunk_df=chunk,
            cache_df=cache_clean,  # Pass the correctly pre-cleaned DataFrame here
            all_gene_sources=cleaned_global_sources,
            n_permutations=n_permutations,
            seed=seed + i  # Shift the seed per block to avoid identical random states across cores
        )
        for i, chunk in enumerate(chunks)
    )

    # 4. Merge results and preserve original indices
    scores_df = pd.concat(processed_dfs, ignore_index=True)
    scores_df.index = data_df.index
    return scores_df


def performance_indices_tfbs(input, env, options, args):
    out_path, temp_path = ut.get_exp_path(input, env)
    bio_path = Path(os.path.expandvars(env["DATA_PATH"]))
    rerun = args.rerun
    n_jobs = args.njobs if args.njobs is not None else 4
    src,tgt = ut.get_input_source_list(out_path)


    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))
    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)
    all_results = [cl["id"] for cl in input["clustering"]]
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list

    cache = pd.read_parquet(bio_path / "tfbs"/ "trap_scores.parquet")
    # print(cache)
    # return 
    # 1. Initialize a dictionary to hold DataFrames for each result ID
    # This prevents constant disk I/O for the cluster files
    df_store = {}
    for r in result_list_input:
        result_file_path = cluster_folder / f"index_tfbs_{r}.csv"
        if result_file_path.exists() and not rerun:
            continue
        cluster_file = cluster_folder / f"{r}.csv"
        df_store[r] = pd.read_csv(cluster_file)
    
    for r, data_df in df_store.items():
        print(f"Computing tfbs scores for cluster set: {r}")
        scores_df = compute_tfbs_scores_df(data_df,cache,src,200,n_jobs)
        df_store[r] = pd.concat([df_store[r], scores_df], axis=1)

    # Final Step: Save consolidated results
    for r, final_df in df_store.items():
        index_file = cluster_folder / f"index_tfbs_{r}.csv"
        final_df.to_csv(index_file, index=False)
        print(f"Saved consolidated results to {index_file}")


# =======

def compute_freq_scores(
    data_df: pd.DataFrame, 
    ge_data: pd.DataFrame, 
    src: list, 
    tgt: list, 
    flat_cache, 
    gene_map, 
    total_genes,
    data_path,
    n_jobs=4
    ):
    """
    Computes distance correlation metrics (dCor_source, dCor_target, DC1, DC2) 
    parallelly across all modules defined in data_df.
    
    Parameters:
    -----------
    data_df : pd.DataFrame
        DataFrame containing cluster/module definitions. Must contain columns
        representing the target gene and the list of module member genes.
        Expected columns: 'target' and 'genes' (where 'genes' is a comma-separated string or list).
    ge_data : pd.DataFrame
        The raw gene expression DataFrame.
    src : list of str
        The source gene list.
    tgt : list of str
        The target gene list.
    ss_cache : pd.DataFrame
        Precomputed source-to-source distance correlation matrix.
    st_cache : pd.DataFrame
        Precomputed source-to-target distance correlation matrix.
    n_jobs : int
        Number of worker threads/processes to use.
        
    Returns:
    --------
    scores_df : pd.DataFrame
        DataFrame with columns ['dCor_source', 'dCor_target', 'DC1', 'DC2']
        aligned precisely with the indices of data_df.
    """
    print(f"Parallelizing module scoring over {len(data_df)} items using {n_jobs} cores...")
    
    # Standardize data_df input columns
    # Assumes 'target' contains the target gene and 'genes' contains the module members
    # targets = data_df['target'].tolist()
    
    # Handle potentially stringified lists from CSV format ("geneA,geneB,geneC")
    module_lists = []
    for item in data_df['common_targets']:
        if isinstance(item, str):
            module_lists.append([g.strip() for g in item.split(';')])
        else:
            module_lists.append(list(item))
            
    results1_raw = Parallel(n_jobs=n_jobs)(
        delayed(compute_common_targets_EC_scores)(
            common_target=mod, 
            flat_cache=flat_cache, 
            gene_to_idx=gene_map, 
            N=total_genes
        ) for mod in module_lists )
    results1 = pd.DataFrame(results1_raw, index=data_df.index)
  
    results2 = compute_GO_ORA(
        data=data_df,
        background_gene_list=tgt, 
        data_path=data_path,
        batch=n_jobs
    )

    results3 = compute_PPI_shared_partners(data=data_df,data_path=data_path)
    # print(results3)

    # Reconstruct back into a clean structured DataFrame matching the input index
    scores_df = pd.concat([results1, results2,results3], axis=1)
    return scores_df


def performance_indices_freq(input, env, options, args):
    out_path, temp_path = ut.get_exp_path(input, env)
    bio_path = Path(os.path.expandvars(env["DATA_PATH"]))
    rerun = args.rerun
    n_jobs = args.njobs if args.njobs is not None else 4
    
    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))
    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)
    all_results = [cl["id"] for cl in input["clustering"]]
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list

    ge_data = put.read_dataset(input["dataset"],env)
    src,tgt = ut.get_input_source_list(out_path)
    flat_cache, gene_map, total_genes = get_dataset_dcorr_cache(bio_path,input["dataset"])

    # 1. Initialize a dictionary to hold DataFrames for each result ID
    # This prevents constant disk I/O for the cluster files
    df_store = {}
    for r in result_list_input:
        result_file_path = cluster_folder / f"freq_index_{r}.csv"
        if result_file_path.exists() and not rerun:
            continue
        freq_file = cluster_folder / f"freq_coreg_{r}.csv"
        df_store[r] = pd.read_csv(freq_file)
    
    for r, data_df in df_store.items():
        print(f"Computing dcor scores for cluster set: {r}")
        scores_df = compute_freq_scores(data_df,ge_data,src,tgt,flat_cache, gene_map, total_genes,bio_path)
        df_store[r] = pd.concat([df_store[r], scores_df], axis=1)

    # Final Step: Save consolidated results
    for r, final_df in df_store.items():
        index_file = cluster_folder / f"freq_index_{r}.csv"
        final_df.to_csv(index_file, index=False)
        print(f"Saved consolidated results to {index_file}")
# =======

def performance_indices_combined(input, env,options ,args):
    """
    compute performance indices for cluster results 
    """
    out_path, temp_path = ut.get_exp_path(input,env)
    rerun = args.rerun
    cluster_folder = temp_path/"clusters"

    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))
    all_results = [cl["id"] for cl in input["clustering"]]
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list


    dataset_name  = input["dataset"]["name"]
    
    file_name = options.get("name","default")

    result_file =  out_path/ f"cluster_indices_{file_name}.parquet"
    if result_file.exists() and not rerun:
        print("already exists")
        return


    res = []

    for r in result_list_input:
        print(r)
        index_file = cluster_folder / f"index_{r}.csv"
        if not index_file.exists():
            raise ValueError(f"file {index_file} note yet generated")
        
        data = pd.read_csv(index_file)

        network_file = cluster_folder / f"index_network_{r}.csv"
        if not network_file.exists():
            raise ValueError(f"file {network_file} note yet generated")
        
        network_score = pd.read_csv(network_file)
        network_cols = [
            "density_hippie", "density_score_hippie", "lcc_hippie", "lcc_score_hippie", 
            "tc_hippie", "tc_score_hippie", "node_found_ratio_hippie", "density_stringdb", 
            "density_score_stringdb", "lcc_stringdb", "lcc_score_stringdb", "tc_stringdb", 
            "tc_score_stringdb", "node_found_ratio_stringdb", "density_biogrid", 
            "density_score_biogrid", "lcc_biogrid", "lcc_score_biogrid", "tc_biogrid", 
            "tc_score_biogrid", "node_found_ratio_biogrid"
        ]

        dcor_file = cluster_folder / f"index_dcorr_{r}.csv"
        if not dcor_file.exists():
            raise ValueError(f"file {dcor_file} note yet generated")
        
        dcorr_score = pd.read_csv(dcor_file)
        dcorr_cols = ["dCor_sources","dCor_target","DC1","DC2"]


        tfbs_file = cluster_folder / f"index_tfbs_{r}.csv"
        if not tfbs_file.exists():
            raise ValueError(f"file {tfbs_file} note yet generated")
        
        tfbs_score = pd.read_csv(tfbs_file)
        tfbs_cols = ["TFBS_affinity","TF_found_per","TFBS_affinity_score"]
        # cluster_uid,target,sources,n_source,silhouette_score,cluster_density,cluster_diameter,note,cluster_note,TFBS_affinity,TF_found_per,TFBS_affinity_score

        # Merge the network columns into 'data'
        # If network_score has 1 row and data has many, this broadcasts the values
        # If they align row-for-row, this joins them side-by-side
        data = pd.concat([data, network_score[network_cols],dcorr_score[dcorr_cols],tfbs_score[tfbs_cols]], axis=1)

        data["config_name"] = r
        res.append(data)

    combined_df = pd.concat(res, ignore_index=True)
    combined_df["dataset"] = dataset_name
    combined_df.to_csv(result_file,index=False)
    combined_df.to_parquet(result_file, index=False, compression="gzip")
    print("saved")


def performance_indices_freq_combined(input, env,options ,args):
    """
    compute performance indices for frequently found clusters  
    """
    out_path, temp_path = ut.get_exp_path(input,env)
    rerun = args.rerun
    cluster_folder = temp_path/"clusters"

    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))
    all_results = [cl["id"] for cl in input["clustering"]]
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list


    dataset_name  = input["dataset"]["name"]
    
    file_name = options.get("name","default")

    result_file =  out_path/ f"freq_coreg_indices_{file_name}.csv.gz"
    if result_file.exists() and not rerun:
        print("already exists")
        return

    res = []

    target_columns = [
        "source", "support_size", "common_targets", "note",
        "dCor_targets", "dCor_targets_score", 
        "go_score", "go_pvalue", "go_classes",
        "ppi_shared_partners_hippie", "ppi_shared_partners_stringdb", "ppi_shared_partners_biogrid"
    ]

    for r in result_list_input:
        print(r)
        index_file = cluster_folder / f"freq_index_{r}.csv"
        if not index_file.exists():
            raise ValueError(f"file {index_file} note yet generated")
        
        data = pd.read_csv(index_file)
        data = data.loc[:, ~data.columns.duplicated(keep="first")]

        # Filter down strictly to the clean target schema columns that exist
        available_cols = [c for c in target_columns if c in data.columns]
        data = data[available_cols]

        data["config_name"] = r
        res.append(data)

    combined_df = pd.concat(res, ignore_index=True)
    combined_df["dataset"] = dataset_name
    combined_df.to_csv(result_file,index=False)
    print("saved")