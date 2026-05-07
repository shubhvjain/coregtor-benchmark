
import  src.results.util as ut 
import pandas as pd
import numpy as np
from scipy.stats import skew

import seaborn as sns
import matplotlib.pyplot as plt

def get_gene_stats(df):
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
      stat = get_gene_stats(gf)
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