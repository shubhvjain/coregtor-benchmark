
import  src.results.util as ut 
import pandas as pd
import numpy as np
from scipy.stats import skew

import seaborn as sns
import matplotlib.pyplot as plt

def _get_gene_stats(df):
    """
    Analyzes a gene frequency histogram DataFrame.
    Returns a flat dictionary of results.
    """
    # 1. Dimensions & Global Counts
    total_rows = df.shape[0]
    total_cols = df.shape[1]
    total_cells = df.size
    
    # 2. Activity & Sparsity
    # active_source_pcr: % of columns that have at least one non-zero
    active_cols_count = (df != 0).any(axis=0).sum()
    active_source_pcr = (active_cols_count / total_cols * 100) if total_cols > 0 else 0
    
    # sparsity_pcr: % of total cells in the matrix that are zero
    zero_cells_count = (df == 0).sum().sum()
    sparsity_pcr = (zero_cells_count / total_cells * 100) if total_cells > 0 else 0
    
    # 3. Frequency Value Distribution
    flat_values = df.values.flatten()
    active_vals = flat_values[flat_values > 0]
    
    # Default values for empty datasets
    c_min = c_max = c_avg = c_med = c_std = c_skew = 0
    
    if active_vals.size > 0:
        c_max = active_vals.max()
        c_avg = active_vals.mean()
        c_med = np.median(active_vals)
        c_std = active_vals.std()
        
        # Skewness requires at least 3 values to be meaningful
        try:
            c_skew = skew(active_vals) if active_vals.size > 2 else 0
        except:
            c_skew = 0

    # 4. Final Flat Dictionary
    stats = {
        "active_source_pcr": round(float(active_source_pcr),2),
        "sparsity_pcr": round(float(sparsity_pcr),2),
        "total_root_nodes": int(total_rows),
        "total_sources": int(total_cols),
        "count_min": float(c_min),
        "count_max": float(c_max),
        "count_avg": float(c_avg),
        "count_median": float(c_med),
        "count_std": round(float(c_std),5),
        "count_skew": round(float(c_skew),5)
    }
    
    return stats

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



def plot_skewness_dist(df, ax=None, save_path=None):
    """
    Returns a KDE plot of skewness. 
    Integrates into subplots if 'ax' is provided.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    
    sns.kdeplot(data=df, x='count_skew', hue='dataset', fill=True, ax=ax)
    ax.set_title("Gene Frequency Skewness Distribution")
    
    if save_path:
        plt.savefig(save_path)
    return ax

def plot_metric_boxplot(df, y_col, title=None ,ax=None, save_path=None):
    """
    Returns a boxplot for any metric (active_source_pcr, root_count, etc.)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    
    sns.boxplot(data=df, x='dataset', y=y_col, ax=ax)

    if title is None:
      title = f"Target-wise {y_col.replace('_', ' ').title()}"
    ax.set_title(title)
    
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
        "min": float(np.min(upper_tri_values)),
        "max": float(np.max(upper_tri_values)),
        "mean": float(np.mean(upper_tri_values)),
        "median": float(np.median(upper_tri_values)),
        "std": float(np.std(upper_tri_values))
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