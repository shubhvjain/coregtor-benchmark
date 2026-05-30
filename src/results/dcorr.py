import numpy as np
import pandas as pd
import dcor
from joblib import Parallel, delayed

def compute_distance_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """
    Computes the sample distance correlation (dCor) between two vectors.
    """
    return float(dcor.distance_correlation(x, y))


def generate_dcor_cache(df_expression: pd.DataFrame, source_gene_list: list, target_gene_list: list):
    """
    Generates a 2D source-target cache matrix and a 2D source-source cache matrix 
    directly from the expression DataFrame.
    
    Parameters:
    -----------
    df_expression : pd.DataFrame
        DataFrame where index is sample_name and columns are gene_name strings.
    source_gene_list : list of str
        List of gene names representing the ~3k pool of source genes.
    target_gene_list : list of str
        List of gene names representing the target genes.
        
    Returns:
    --------
    source_target_cache : np.ndarray
        2D array of shape (num_sources, num_targets), where entry [i, j]
        is dCor(source_gene_list[i], target_gene_list[j])
    source_source_cache : np.ndarray
        2D symmetric matrix of shape (num_sources, num_sources) for pairwise 
        source relationships.
    """
    # 1. Filter lists to ensure genes actually exist in the loaded dataset columns
    valid_sources = [g for g in source_gene_list if g in df_expression.columns]
    valid_targets = [g for g in target_gene_list if g in df_expression.columns]
    
    num_sources = len(valid_sources)
    num_targets = len(valid_targets)
    
    print(f"Generating caches for {num_sources} valid sources and {num_targets} valid targets.")
    
    # Initialize cache arrays
    source_target_cache = np.zeros((num_sources, num_targets))
    source_source_cache = np.zeros((num_sources, num_sources))
    
    # Extract underlying NumPy matrices for fast indexing inside loops
    # .to_numpy() bypasses slower pandas overhead inside the loops
    source_matrix = df_expression[valid_sources].to_numpy()
    target_matrix = df_expression[valid_targets].to_numpy()
    
    # 2. Populate Source-to-Target Cache Matrix
    for j in range(num_targets):
        target_vector = target_matrix[:, j]
        
        for i in range(num_sources):
            source_vector = source_matrix[:, i]
            source_target_cache[i, j] = compute_distance_correlation(source_vector, target_vector)
            
    # 3. Populate Source-to-Source Cache Matrix
    for i in range(num_sources):
        gene_i_vector = source_matrix[:, i]
        source_source_cache[i, i] = 1.0  # Self-correlation identity
        
        for j in range(i + 1, num_sources):
            gene_j_vector = source_matrix[:, j]
            score = compute_distance_correlation(gene_i_vector, gene_j_vector)
            
            # Symmetrical assignment
            source_source_cache[i, j] = score
            source_source_cache[j, i] = score
            
    return source_source_cache,source_target_cache


def _compute_target_row(i: int, source_matrix: np.ndarray, target_matrix: np.ndarray) -> np.ndarray:
    """
    Worker function to compute all target distance correlations for a single source gene.
    Computing an entire row at once reduces parallel scheduling overhead.
    """
    num_targets = target_matrix.shape[1]
    row_results = np.zeros(num_targets)
    source_vector = source_matrix[:, i]
    
    for j in range(num_targets):
        row_results[j] = float(dcor.distance_correlation(source_vector, target_matrix[:, j]))
    return row_results


def _compute_source_row(i: int, source_matrix: np.ndarray) -> np.ndarray:
    """
    Worker function to compute distance correlations for a single source gene i
    against all other source genes j where j >= i (upper triangle).
    """
    num_sources = source_matrix.shape[1]
    row_results = np.zeros(num_sources)
    source_i_vector = source_matrix[:, i]
    
    # Self-correlation identity
    row_results[i] = 1.0
    
    # Only calculate the upper triangle
    for j in range(i + 1, num_sources):
        row_results[j] = float(dcor.distance_correlation(source_i_vector, source_matrix[:, j]))
    return row_results


def generate_dcor_caches_parallel(df_expression: pd.DataFrame, source_gene_list: list, target_gene_list: list, n_jobs: int = -1):
    """
    Generates a 2D source-target cache matrix and a 2D source-source cache matrix 
    by distributing computations across cores using joblib.
    
    Parameters:
    -----------
    df_expression : pd.DataFrame
        DataFrame where rows are samples and columns are gene_name strings.
    source_gene_list : list of str
        List of gene names representing the pool of source genes (~3k).
    target_gene_list : list of str
        List of gene names representing the target genes.
    n_jobs : int, default=-1
        The number of CPUs to use. -1 means use all available processors.
        
    Returns:
    --------
    source_target_cache : np.ndarray
        2D array of shape (num_sources, num_targets)
    source_source_cache : np.ndarray
        2D symmetric matrix of shape (num_sources, num_sources)
    valid_sources : list of str
    valid_targets : list of str
    """
    # Filter lists to ensure genes exist in dataset columns
    valid_sources = [g for g in source_gene_list if g in df_expression.columns]
    valid_targets = [g for g in target_gene_list if g in df_expression.columns]
    
    num_sources = len(valid_sources)
    num_targets = len(valid_targets)
    
    print(f"Parallelizing cache generation via joblib using {n_jobs} core(s)...")
    print(f"Sources: {num_sources}, Targets: {num_targets}")
    
    # Convert to contiguous NumPy arrays. 
    # Joblib automatically memmaps these when spawning processes.
    source_matrix = np.ascontiguousarray(df_expression[valid_sources].to_numpy())
    target_matrix = np.ascontiguousarray(df_expression[valid_targets].to_numpy())
    
    # Initialize cache structures
    source_target_cache = np.zeros((num_sources, num_targets))
    source_source_cache = np.zeros((num_sources, num_sources))
    
    # 1. Compute Source-to-Target Cache
    print("Computing source-to-target cache...")
    target_rows = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_compute_target_row)(i, source_matrix, target_matrix) for i in range(num_sources)
    )
    # Stack row results into the final 2D array
    source_target_cache = np.vstack(target_rows)
    
    # 2. Compute Source-to-Source Cache (Upper Triangle)
    print("Computing source-to-source cache...")
    source_rows = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_compute_source_row)(i, source_matrix) for i in range(num_sources)
    )
    
    #Reconstruct symmetric matrix from the upper triangular row outputs
    for i, row in enumerate(source_rows):
        # Assign upper triangle values calculated by the worker
        source_source_cache[i, i:] = row[i:]
        # Mirror values to the lower triangle
        source_source_cache[i:, i] = row[i:]
            
    # return source_source_cache,source_target_cache
    # Wrap raw matrices into beautifully indexed DataFrames
    df_source_target = pd.DataFrame(
        source_target_cache, 
        index=valid_sources, 
        columns=valid_targets
    )
    
    df_source_source = pd.DataFrame(
        source_source_cache, 
        index=valid_sources, 
        columns=valid_sources
    )
    
    return df_source_source,df_source_target

def score_single_module(module_genes: list, target_gene: str, df_src_src: pd.DataFrame, df_src_tgt: pd.DataFrame, R: int = 1000) -> dict:
    """
    Computes the observed scores, runs permutation testing, and returns 
    the four metrics (dCor_source, dCor_target, DC1, DC2) for a single module.
    
    Parameters:
    -----------
    module_genes : list of str
        List of gene names belonging to the co-regulatory module M.
    target_gene : str
        The target gene name t.
    df_src_src : pd.DataFrame
        The precomputed source-to-source cache DataFrame.
    df_src_tgt : pd.DataFrame
        The precomputed source-to-target cache DataFrame.
    R : int
        Number of permutations.
    """
    # 1. Filter module genes against cache coverage
    valid_m = [g for g in module_genes if g in df_src_src.index]
    m_size = len(valid_m)
    
    if m_size < 2 or target_gene not in df_src_tgt.columns:
        return {"dCor_source": np.nan, "dCor_target": np.nan, "DC1": np.nan, "DC2": np.nan}
    
    # 2. Extract underlying raw matrices and coordinate index integers
    src_src_matrix = df_src_src.to_numpy()
    src_tgt_matrix = df_src_tgt.to_numpy()
    
    # Map names to integer positions
    src_indices = [df_src_src.index.get_loc(g) for g in valid_m]
    tgt_idx = df_src_tgt.columns.get_loc(target_gene)
    all_possible_sources = np.arange(src_src_matrix.shape[0])
    
    # 3. Compute Observed Metrics
    # Source-to-Target Observed Mean
    obs_target = np.mean(src_tgt_matrix[src_indices, tgt_idx])
    
    # Source-to-Source Observed Mean (extracting unique upper triangular pairs)
    sub_matrix = src_src_matrix[np.ix_(src_indices, src_indices)]
    obs_source = np.sum(np.triu(sub_matrix, k=1)) / (m_size * (m_size - 1) / 2)
    
    # 4. Permutation Testing Loop
    better_target_count = 0
    better_source_count = 0
    
    # Pre-generate random index choices to maximize internal CPU caching
    # Drawing m_size random source genes uniformly at random without replacement
    for _ in range(R):
        perm_indices = np.random.choice(all_possible_sources, size=m_size, replace=False)
        
        # Calculate Permuted Target Mean
        perm_target = np.mean(src_tgt_matrix[perm_indices, tgt_idx])
        if perm_target >= obs_target:
            better_target_count += 1
            
        # Calculate Permuted Source Mean
        perm_sub = src_src_matrix[np.ix_(perm_indices, perm_indices)]
        perm_source = np.sum(np.triu(perm_sub, k=1)) / (m_size * (m_size - 1) / 2)
        if perm_source >= obs_source:
            better_source_count += 1
            
    # 5. Evaluate final P-values and Log Scores
    p_target = (better_target_count + 1) / (R + 1)
    p_source = (better_source_count + 1) / (R + 1)
    
    dc1 = -np.log(p_target)
    dc2 = -np.log(p_source)
    
    return {
        "dCor_source": obs_source,
        "dCor_target": obs_target,
        "DC1": dc1,
        "DC2": dc2
    }

