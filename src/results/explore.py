
import  src.results.util as ut 
import pandas as pd
import numpy as np
from scipy.stats import skew

import seaborn as sns
import matplotlib.pyplot as plt

### Gene frequency histogram stats

import numpy as np
import pandas as pd
from scipy.stats import skew

def _get_gene_stats(df):
    """
    Analyzes a gene frequency histogram DataFrame.
    Returns a flat dictionary of results including context density,
    entropy, evenness, peer density, and peer ranks.
    """
    # 1. Dimensions & Global Counts
    total_rows = df.shape[0]  # r (Number of root genes)
    total_cols = df.shape[1]  # s (Number of source genes)
    total_cells = df.size
    
    # 2. Activity & Sparsity
    active_cols_count = (df != 0).any(axis=0).sum()
    active_source_pcr = (active_cols_count / total_cols * 100) if total_cols > 0 else 0
    
    zero_cells_count = (df == 0).sum().sum()
    sparsity_pcr = (zero_cells_count / total_cells * 100) if total_cells > 0 else 0
    
    # 3. Frequency Value Distribution (Existing Logic)
    flat_values = df.values.flatten()
    active_vals = flat_values[flat_values > 0]
    
    c_min = c_max = c_avg = c_med = c_std = c_skew = 0
    if active_vals.size > 0:
        c_max = active_vals.max()
        c_avg = active_vals.mean()
        c_med = np.median(active_vals)
        c_std = active_vals.std()
        try:
            c_skew = skew(active_vals) if active_vals.size > 2 else 0
        except:
            c_skew = 0

    # =========================================================================
    #  Row-wise Analysis (Context, Entropy, Peer behaviors)
    # =========================================================================
    
    # Identify which columns (source genes) are also root genes (rows)
    # This assumes df.index contains root gene names and df.columns contains source gene names
    root_genes_set = set(df.index)
    is_peer_col = np.array([col in root_genes_set for col in df.columns])
    
    # Pre-allocate arrays for row-wise metrics
    context_densities = []
    entropies = []
    evennesses = []
    peer_densities = []
    mean_peer_ranks = []
    
    # Calculate H_max = log2(s)
    h_max = np.log2(total_cols) if total_cols > 0 else 0
    
    # Convert DataFrame to numpy for faster row-by-row iteration
    matrix = df.values 
    
    for i in range(total_rows):
        row_frequencies = matrix[i, :]
        row_sum = row_frequencies.sum()
        
        # --- Context Density ---
        non_zero_mask = row_frequencies > 0
        non_zero_count = non_zero_mask.sum()
        c_density_i = (non_zero_count / total_cols) if total_cols > 0 else 0
        context_densities.append(c_density_i)
        
        # --- Shannon Entropy & Evenness ---
        if row_sum > 0 and total_cols > 0:
            # Normalize to discrete probability distribution p_ij
            p_i = row_frequencies / row_sum
            # Filter out zeros to avoid log2(0)
            p_i_active = p_i[p_i > 0]
            entropy_i = -np.sum(p_i_active * np.log2(p_i_active))
            evenness_i = (entropy_i / h_max) if h_max > 0 else 0
        else:
            entropy_i = 0
            evenness_i = 0
            
        entropies.append(entropy_i)
        evennesses.append(evenness_i)
        
        # --- Peer Density & Peer Ranks ---
        # Peer columns with non-zero values in this specific row (R_i)
        active_peers_mask = non_zero_mask & is_peer_col
        r_i_count = active_peers_mask.sum()
        
        # peer_density(i) = |R_i| / r
        p_density_i = (r_i_count / total_rows) if total_rows > 0 else 0
        peer_densities.append(p_density_i)
        
        # Mean Peer Rank
        if r_i_count > 0:
            # Rank all genes in the row (method='dense', descending order: highest freq gets rank 1)
            # Scipy/Pandas equivalents can be slow in a loop, so we use a clean pandas rank approach per row
            row_series = pd.Series(row_frequencies)
            # We rank descending so the highest frequency has rank 1.
            ranks = row_series.rank(method='dense', ascending=False).values
            
            # Extract ranks of active peer genes
            peer_ranks_i = ranks[active_peers_mask]
            mean_peer_rank_i = peer_ranks_i.mean()
        else:
            mean_peer_rank_i = 0
            
        mean_peer_ranks.append(mean_peer_rank_i)
        
    # --- Matrix-level aggregations (Averages over all r root nodes) ---
    mean_context_density = np.mean(context_densities) if total_rows > 0 else 0
    mean_entropy = np.mean(entropies) if total_rows > 0 else 0
    mean_evenness = np.mean(evennesses) if total_rows > 0 else 0
    mean_peer_density = np.mean(peer_densities) if total_rows > 0 else 0
    mean_peer_rank_m = np.mean(mean_peer_ranks) if total_rows > 0 else 0

    # =========================================================================
    # 5. Final Flat Dictionary
    # =========================================================================
    stats = {
        "active_source_pcr": round(float(active_source_pcr), 2),
        "sparsity_pcr": round(float(sparsity_pcr), 2),
        "total_root_nodes": int(total_rows),
        "total_sources": int(total_cols),
        "count_min": round(float(c_min), 5),
        "count_max": round(float(c_max), 5),
        "count_avg": round(float(c_avg), 5),
        "count_median": round(float(c_med), 5),
        "count_std": round(float(c_std), 5),
        "count_skew": round(float(c_skew), 5),
        
        # New Additions
        "mean_context_density": round(float(mean_context_density), 5),
        "mean_entropy": round(float(mean_entropy), 5),
        "mean_evenness_index": round(float(mean_evenness), 5),
        "mean_peer_density": round(float(mean_peer_density), 5),
        "mean_peer_rank": round(float(mean_peer_rank_m), 5)
    }
    
    return stats

# def _get_gene_stats(df):
#     """
#     Analyzes a gene frequency histogram DataFrame.
#     Returns a flat dictionary of results.
#     """
#     # 1. Dimensions & Global Counts
#     total_rows = df.shape[0]
#     total_cols = df.shape[1]
#     total_cells = df.size
    
#     # 2. Activity & Sparsity
#     # active_source_pcr: % of columns that have at least one non-zero
#     active_cols_count = (df != 0).any(axis=0).sum()
#     active_source_pcr = (active_cols_count / total_cols * 100) if total_cols > 0 else 0
    
#     # sparsity_pcr: % of total cells in the matrix that are zero
#     zero_cells_count = (df == 0).sum().sum()
#     sparsity_pcr = (zero_cells_count / total_cells * 100) if total_cells > 0 else 0
    
#     # 3. Frequency Value Distribution
#     flat_values = df.values.flatten()
#     active_vals = flat_values[flat_values > 0]
    
#     # Default values for empty datasets
#     c_min = c_max = c_avg = c_med = c_std = c_skew = 0
    
#     if active_vals.size > 0:
#         c_max = active_vals.max()
#         c_avg = active_vals.mean()
#         c_med = np.median(active_vals)
#         c_std = active_vals.std()
        
#         # Skewness requires at least 3 values to be meaningful
#         try:
#             c_skew = skew(active_vals) if active_vals.size > 2 else 0
#         except:
#             c_skew = 0

#     # 4. Final Flat Dictionary
#     stats = {
#         "active_source_pcr": round(float(active_source_pcr),2),
#         "sparsity_pcr": round(float(sparsity_pcr),2),
#         "total_root_nodes": int(total_rows),
#         "total_sources": int(total_cols),
#         "count_min": round(float(c_min),5),
#         "count_max": round(float(c_max),5),
#         "count_avg": round(float(c_avg),5),
#         "count_median": round(float(c_med),5),
#         "count_std": round(float(c_std),5),
#         "count_skew": round(float(c_skew),5)
#     }
    
#     return stats

## old 
def gene_frequency_analysis(input, env,options ,args):
  """
  generate per target gene frequency stats for the whole dataset
  """
  out_path, temp_path = ut.get_exp_path(input,env)
  rerun = args.rerun

  gf_file_path = out_path/"gene_frequency.csv"
  if gf_file_path.exists() and not rerun:
      print("File already exists")
      return

  stats = []
  success_target = ut.get_exp_target_list(out_path)
  for t in success_target:
      data = ut.get_temp_file(temp_path,t)
      gf = data["results"]["gene_frequency"]
      stat = _get_gene_stats(gf)
      stat["target"] = t
      stats.append(stat)

  res = pd.DataFrame(stats)
  res["exp"] = input["id"]
  res["dataset"] = input["dataset"]["name"] 
  # print(res)
  cols_to_front = ["exp", "dataset", "target"]
  other_cols = [c for c in res.columns if c not in cols_to_front]
  res = res[cols_to_front + other_cols]
  res.to_csv(gf_file_path,index=False)
  save_plot_gene_frequency_dataset(res,out_path)



def plot_skewness_dist(df, ax=None, col_name="count_skew",save_path=None):
    """
    Returns a KDE plot of skewness. 
    Integrates into subplots if 'ax' is provided.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    
    sns.kdeplot(data=df, x=col_name, hue='dataset', fill=True, ax=ax)
    ax.set_title("Gene Frequency Skewness")
    
    if save_path:
        plt.savefig(save_path)
    return ax

def plot_metric_boxplot(df, y_col, title=None, ylabel=None,ax=None, save_path=None):
    """
    Returns a boxplot for any metric (active_source_pcr, root_count, etc.)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    
    sns.boxplot(data=df, x='dataset', y=y_col, ax=ax)

    if title is None:
      title = f"Target-wise {y_col.replace('_', ' ').title()}"
    ax.set_title(title)
    
    if ylabel is None:
        ylabel = y_col.replace('_', ' ')
    ax.set_ylabel(ylabel)

    if save_path:
        plt.savefig(save_path)
    
    return ax



def save_plot_gene_frequency_dataset(res, out_file):
    # Setup 4 subplots with width ratios [2, 2, 2, 4] to match [20, 20, 20, 40]
    fig, axes = plt.subplots(
        1, 4, 
        figsize=(12,4), 
        gridspec_kw={'width_ratios': [2, 2, 2, 4]}
    )

    # 1. Total Root Nodes (Boxplot)
    plot_metric_boxplot(res, y_col='total_root_nodes',title="Root node count", ax=axes[0])
    
    # 2. Active Source Percentage (Boxplot)
    plot_metric_boxplot(res, y_col='active_source_pcr',title="Active sources in context(%)", ax=axes[1])
    
    # 3. Sparsity Percentage (Boxplot)
    plot_metric_boxplot(res, y_col='sparsity_pcr',title="Zero count Percentage", ax=axes[2])
    
    # 4. Skewness Distribution (KDE)
    plot_skewness_dist(res, ax=axes[3])

    plt.tight_layout()
    # Save as SVG 
    plt.savefig(out_file / "gene_frequency.svg")


### Distance stats

def _get_dist_stats(df):
    # Ensure we are working with a NumPy array
    matrix = df.values if isinstance(df, pd.DataFrame) else np.array(df)
    
    # Get indices for the upper triangle, excluding the diagonal (k=1)
    # Number of elements = [n * (n - 1)] / 2
    tri_indices = np.triu_indices_from(matrix, k=1)
    upper_tri_values = matrix[tri_indices]
    
    # Check if the array is empty (happens if n_roots < 2)
    if upper_tri_values.size == 0:
        return {
            "n_roots": int(matrix.shape[0]),
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "std": None
        }
    
    # Calculate stats and force conversion to standard Python types for JSON/CSV
    stats = {
        "n_roots": int(matrix.shape[0]),
        "min": round(float(np.min(upper_tri_values)),5),
        "max": round(float(np.max(upper_tri_values)),5),
        "mean": round(float(np.mean(upper_tri_values)),5),
        "median": round(float(np.median(upper_tri_values)),5),
        "std": round(float(np.std(upper_tri_values)),5)
    }   
    return stats


def distance_measure_analysis(input, env,options ,args):
    """
    """
    out_path, temp_path = ut.get_exp_path(input,env)
    rerun = args.rerun

    dm_file_path = out_path/"distance.csv"
    if dm_file_path.exists() and not rerun:
        print("File already exists")
        return

    stats = []
    success_target = ut.get_exp_target_list(out_path)
    dmeasures = input["context_comparison"]["methods"]
    for t in success_target:
        data = ut.get_temp_file(temp_path,t)
        for dms in dmeasures: 
            dmatrix = data["results"]["distance"][dms]
            stat = _get_dist_stats(dmatrix)
            stat["target"] = t
            stat["distance"] = dms
            stats.append(stat)
            
        
    res = pd.DataFrame(stats)
    res["exp"] = input["id"]
    res["dataset"] = input["dataset"]["name"] 
    # print(res)
    cols_to_front = ["exp", "dataset", "target"]
    other_cols = [c for c in res.columns if c not in cols_to_front]
    res = res[cols_to_front + other_cols]
    res.to_csv(dm_file_path,index=False)
    save_plot_distance_dataset(res,dmeasures,out_path)


def save_plot_distance_dataset(res, dmeasures, out_file):
    """
    Generates a grid of boxplots: 
    Rows = Distance Measures (e.g., Euclidean, Cosine)
    Cols = Min, Mean, and Max distributions across all targets
    """
    n_measures = len(dmeasures)
    
    # Grid: Rows = n_measures, Cols = 3 (Min, Mean, Max)
    fig, axes = plt.subplots(
        nrows=n_measures, 
        ncols=3, 
        figsize=(12, 4 * n_measures), 
        squeeze=False 
    )

    # Define columns to plot and their visual settings
    cols_to_plot = [
        ('min', 'Min Distance', 'skyblue'),
        ('mean', 'Average Distance', 'salmon'),
        ('max', 'Max Distance', 'lightgreen')
    ]

    for i, dms in enumerate(dmeasures):
        subset = res[res['distance'] == dms]
        
        for j, (col_name, label, color) in enumerate(cols_to_plot):
            ax = axes[i, j]
            sns.boxplot(data=subset, x='dataset', y=col_name, ax=ax, color=color)
            
            # Formatting
            ax.set_title(f"{dms}: {label}")
            ax.set_ylabel(label)
            ax.set_xlabel("") # Keeps the x-axis clean
            
            # Optional: Add a grid for better readability of distance scales
            ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    
    # Save as SVG
    save_path = out_file / "distance_measures.svg"
    plt.savefig(save_path)
    print(f"Comprehensive plot saved to {save_path}")
    plt.close()