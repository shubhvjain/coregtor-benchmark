import numpy as np
import pandas as pd
import dcor
from joblib import Parallel, delayed
import os
import json
import pyarrow as pa
import pyarrow.parquet as pq

def compute_distance_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """
    Computes the sample distance correlation (dCor) between two vectors.
    """
    return float(dcor.distance_correlation(x, y))


def _compute_combined_triangular_row(i: int, combined_matrix: np.ndarray) -> np.ndarray:
    """
    Computes a single row of the strict upper triangle for the combined master matrix.
    Excludes the diagonal to match analytical 1D index mapping requirements.
    """
    num_genes = combined_matrix.shape[1]
    
    # If it's the last gene, there are no remaining upper triangle elements
    if i >= num_genes - 1:
        return np.array([], dtype=np.float64)
        
    # Allocate only for strict upper triangle items (excluding diagonal)
    row_len = num_genes - 1 - i
    row_results = np.zeros(row_len)
    reference_vector = combined_matrix[:, i]
    
    for idx, j in enumerate(range(i + 1, num_genes)):
        row_results[idx] = float(dcor.distance_correlation(reference_vector, combined_matrix[:, j], method='mergesort'))
    return row_results


def generate_combined_dcor_cache(df_expression: pd.DataFrame, source_gene_list: list, target_gene_list: list, n_jobs: int = -1):
    """
    Computes a flat 1D upper-triangular array of pairwise dCor values across all unique 
    genes. Returns the 1D cache and a dictionary mapping gene names to matrix indices.
    """
    # Handle duplicate gene columns by keeping the one with the highest variance
    if df_expression.columns.duplicated().any():
        print("Duplicate genes detected. Selecting representatives with highest variance...")
        
        variances = df_expression.var(axis=0)
        unique_highest_var_genes = variances.sort_values(ascending=False).index.drop_duplicates(keep='first')
        
        col_positions = [df_expression.columns.get_loc(gene) for gene in unique_highest_var_genes]
        final_positions = [pos[0] if isinstance(pos, np.ndarray) or isinstance(pos, slice) else pos for pos in col_positions]
        
        df_expression = df_expression.iloc[:, final_positions]

    # 2. Consolidate a unique master pool of all genes involved
    master_genes = list(set(source_gene_list + target_gene_list))
    valid_genes = [g for g in master_genes if g in df_expression.columns]
    num_genes = len(valid_genes)
    
    # Map gene strings to integer array positions for fast downstream lookup
    gene_to_idx = {gene: idx for idx, gene in enumerate(valid_genes)}
    
    print(f"Total Master Genes for Permutation Pool: {num_genes}")
    print(f"Parallelizing upper triangle generation using {n_jobs} core(s)...")
    
    # Extract contiguous NumPy expression matrix (Samples x Genes)
    expr_matrix = np.ascontiguousarray(df_expression[valid_genes].to_numpy())
    
    # 3. Parallel Pass over Rows using the exact method name provided
    raw_rows = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_compute_combined_triangular_row)(i, expr_matrix) for i in range(num_genes)
    )
    
    # 4. Concatenate directly into a single flat 1D array
    flat_dcor_cache = np.concatenate([row for row in raw_rows if len(row) > 0])
    
    print(f"Cache generated successfully. Size in RAM: {flat_dcor_cache.nbytes / (1024**3):.2f} GB")
    
    return flat_dcor_cache, gene_to_idx, num_genes


def compute_module_scores(
    module_genes: list, 
    target_gene: str, 
    source_gene_set: list,
    flat_cache: np.ndarray, 
    gene_to_idx: dict, 
    N: int, 
    R: int = 500):
    """
    Computes observed dCor scores and executes permutation testing to calculate
    DC1(M) and DC2(M), drawing random modules strictly from the provided source_gene_set.
    
    Args:
        module_genes (list): List of gene names in the module M.
        target_gene (str): Gene name of the target t.
        source_gene_set (list): The master list of valid source genes (e.g., TFs/co-regulators).
        flat_cache (np.ndarray): Loaded 1D flat upper-triangular array.
        gene_to_idx (dict): Gene-to-index mapping dictionary from JSON metadata.
        N (int): Total number of unique genes in the cache pool (num_genes).
        R (int): Number of permutations.
    """
    # Helper to convert gene pair (g1, g2) to the strict flat 1D array index
    def get_flat_idx(g1, g2):
        i, j = gene_to_idx[g1], gene_to_idx[g2]
        if i > j: 
            i, j = j, i  # Symmetrical adjustment for strict upper triangle
        return int(i * N - (i * (i + 1)) // 2 + (j - i - 1))

    # 1. Clean and validate inputs against our cache universe
    all_genes_universe = list(gene_to_idx.keys())
    module_genes = [g for g in module_genes if g in gene_to_idx]
    n = len(module_genes)
    
    # Filter the background permutation pool to include ONLY valid source genes present in the cache
    valid_source_pool = np.array([g for g in source_gene_set if g in gene_to_idx])
    universe_arr = np.array(all_genes_universe)
    
    if n < 2 or target_gene not in gene_to_idx or len(valid_source_pool) < n:
        return {"dCor_sources": 0.0, "dCor_target": 0.0, "DC1": 0.0, "DC2": 0.0, "error": "Invalid pool size or missing indices"}

    # ==========================================
    # 2. COMPUTE OBSERVED METRICS FOR MODULE M
    # ==========================================
    # Ms: Source-to-Source pairs (i < j)
    ms_indices = [get_flat_idx(module_genes[k], module_genes[m]) 
                  for k in range(n) for m in range(k + 1, n)]
    obs_dcor_source = np.mean(flat_cache[ms_indices])
    
    # Mt: Source-to-Target pairs
    mt_indices = [get_flat_idx(g, target_gene) for g in module_genes]
    obs_dcor_target = np.mean(flat_cache[mt_indices])
    
    # ==========================================
    # 3. PERMUTATION TESTING (R ITERATIONS)
    # ==========================================
    better_source_count = 0
    better_target_count = 0
    
    for _ in range(R):
        # CORRECTED STEP: Sample random background module M' STRICTLY from the regulatory source pool
        perm_module = np.random.choice(valid_source_pool, size=n, replace=False)
        
        # Sample random background target t' from the entire gene universe
        perm_target = np.random.choice(universe_arr)
        
        # Calculate Permuted Source-to-Source Mean
        perm_ms_idx = [get_flat_idx(perm_module[k], perm_module[m]) 
                       for k in range(n) for m in range(k + 1, n)]
        perm_dcor_source = np.mean(flat_cache[perm_ms_idx])
        
        # Calculate Permuted Source-to-Target Mean
        perm_mt_idx = [get_flat_idx(g, perm_target) for g in perm_module]
        perm_dcor_target = np.mean(flat_cache[perm_mt_idx])
        
        # Evaluate Indicator Condition: I[m'_i >= m]
        if perm_dcor_source >= obs_dcor_source:
            better_source_count += 1
        if perm_dcor_target >= obs_dcor_target:
            better_target_count += 1
            
    # ==========================================
    # 4. CALCULATE P-VALUES AND DC SCORES
    # ==========================================
    p_source = (better_source_count + 1) / (R + 1)
    p_target = (better_target_count + 1) / (R + 1)
    
    dc1_score = -np.log(p_target)
    dc2_score = -np.log(p_source)
    
    return {
        "dCor_sources": round(obs_dcor_source,5),
        "dCor_target": round(obs_dcor_target,5),
        "DC1": round(dc1_score,5),
        "DC2": round(dc2_score,5)
    }



def compute_common_targets_EC_scores(
    common_target: list, 
    flat_cache: np.ndarray, 
    gene_to_idx: dict, 
    N: int, 
    R: int = 500
) -> dict:
    """
    Computes the internal Expression Coherence (EC) score for a set of common targets
    by comparing their mean pairwise dCor against a randomized background distribution.
    
    Args:
        common_target (list): List of target gene names belonging to a frequent itemset module.
        flat_cache (np.ndarray): Flattened 1D upper-triangular array of pairwise dCor values.
        gene_to_idx (dict): Dictionary mapping gene strings to matrix indices.
        N (int): Total number of unique genes in the cache pool.
        R (int): Number of permutations.
        
    Returns:
        dict: Containing observed mean target dCor, the empirical P-value, and the ECCT score.
    """
    # Helper to convert gene pair (g1, g2) to the strict flat 1D array index
    def get_flat_idx(g1, g2):
        i, j = gene_to_idx[g1], gene_to_idx[g2]
        if i > j: 
            i, j = j, i  # Symmetrical adjustment for strict upper triangle
        return int(i * N - (i * (i + 1)) // 2 + (j - i - 1))

    # 1. Clean and validate target list against the cache pool universe
    all_genes_universe = list(gene_to_idx.keys())
    universe_arr = np.array(all_genes_universe)
    
    valid_targets = [g for g in common_target if g in gene_to_idx]
    n_targets = len(valid_targets)
    
    # Statistical control: Pairwise correlation requires at least 2 elements
    if n_targets < 2:
        return {
            "obs_mean_dCor_targets": 0.0,
            "p_value": 1.0,
            "ECCT": 0.0,
            "error": "Insufficient valid targets in cache universe (minimum 2 required)"
        }

    # ==========================================
    # 2. COMPUTE OBSERVED MEAN PAIRWISE TARGET dCor
    # ==========================================
    # Get all distinct pairs i < j within the target set
    target_pair_indices = [
        get_flat_idx(valid_targets[k], valid_targets[m]) 
        for k in range(n_targets) for m in range(k + 1, n_targets)
    ]
    obs_mean_target_dcor = np.mean(flat_cache[target_pair_indices])
    
    # ==========================================
    # 3. PERMUTATION TESTING (R ITERATIONS)
    # ==========================================
    better_target_set_count = 0
    
    for _ in range(R):
        # Sample random background targets of size n_targets from the whole universe
        perm_targets = np.random.choice(universe_arr, size=n_targets, replace=False)
        
        # Calculate Permuted Target-to-Target Mean
        perm_pair_indices = [
            get_flat_idx(perm_targets[k], perm_targets[m]) 
            for k in range(n_targets) for m in range(k + 1, n_targets)
        ]
        perm_mean_target_dcor = np.mean(flat_cache[perm_pair_indices])
        
        # Evaluate indicator criteria: Did a random shuffle score equal or better?
        if perm_mean_target_dcor >= obs_mean_target_dcor:
            better_target_set_count += 1
            
    # ==========================================
    # 4. CALCULATE EMPIRICAL P-VALUE AND ECCT
    # ==========================================
    p_value = (better_target_set_count + 1) / (R + 1)
    ecct_score = -np.log(p_value)
    
    return {
        "dCor_targets": round(obs_mean_target_dcor, 5),
        "dCor_targets_score": round(ecct_score, 5)
    }