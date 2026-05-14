import  src.results.util as ut 
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt


import re

def natural_keys(text):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text)]

def parameter_selection_plot(input, env, options, args):
  """
  """
  out_path, temp_path = ut.get_exp_path(input, env)  
  index_file = out_path/ f"{options['index_file_name']}_indices.csv.gz"
  diagram_folder = out_path / "ps"
  diagram_folder.mkdir(exist_ok=True,parents=True)
  diagram_pdf = diagram_folder/ f"{options['name']}.pdf"
  rerun = args.rerun

  if not index_file.exists():
    print("index file not found")
    return 
  
  if diagram_pdf.exists() and not rerun:
    print("file already exists")
    return 
  
  optimal_results_only = options.get("optimal_results_only",False)
  default_metrics = ['silhouette_score',  'n_source','shared_PPI_partners_score_stringdb', 'shortest_PPI_path_score_stringdb', 'goa_similarity_jc']
  metrics_list = options.get("metrics_list",default_metrics)

  data = pd.read_csv(index_file)

  plot_landscape_metrics(data,metrics_list,diagram_pdf,optimal_results_only)

   


def plot_landscape_metrics(df, metrics_list, output_name=None, choose_optimal=False):
    """
    Plots a multi-row boxplot comparison for configs.
    Optional: choose_optimal filters df for records where sil_score_optimal is True.
    """
    # --- 1. FILTER FOR OPTIMAL RECORDS IF REQUESTED ---
    if choose_optimal:
        if 'sil_score_optimal' in df.columns:
            # Filtering the dataframe to only include 'True' optimal records
            df = df[df['sil_score_optimal'] == True].copy()
        else:
            print("Warning: 'sil_score_optimal' column not found. Plotting full dataframe.")

    n_rows = len(metrics_list)
    
    INDEX_META = {
        "silhouette_score": {"key": "SIL", "label": "Silhouette score"},
        "n_source": {"key": "NC", "label": "Cluster size"},
        "goa_similarity_lin": {"key": "GO1", "label": "GO Similarity (Lin)"},
        "goa_similarity_resnik": {"key": "GO2", "label": "GO Similarity (Resnik)"},
        "goa_similarity_jc": {"key": "GO3", "label": "GO Similarity (JC)"},
        "shortest_PPI_path_score_hippie": {"key": "PPI2 (H)", "label": "PPI Shortest Path (Hippie)"},
        "shortest_PPI_path_score_stringdb": {"key": "PPI2 (S)", "label": "PPI Shortest Path (StringDB)"},
        "shortest_PPI_path_score_biogrid": {"key": "PPI2 (B)", "label": "PPI Shortest Path (Biogrid)"},
        "shared_PPI_partners_score_hippie": {"key": "PPI1 (H)", "label": "PPI Shared Partners (Hippie)"},
        "shared_PPI_partners_score_stringdb": {"key": "PPI1 (S)", "label": "PPI Shared Partners (StringDB)"},
        "shared_PPI_partners_score_biogrid": {"key": "PPI1 (B)", "label": "PPI Shared Partners (Biogrid)"},
        "grn_collectri_precision": {"key": "PRC", "label": "GRN Precision (Collectri)"},
        "grn_collectri_recall": {"key": "RCL", "label": "GRN Recall (Collectri)"},
        "grn_collectri_jaccard": {"key": "JAC", "label": "GRN Jaccard (Collectri)"},
    }
    
    fig, axes = plt.subplots(nrows=n_rows, ncols=1, figsize=(11.69, 8.27))
    if n_rows == 1: axes = [axes]

    # Ensure configs are sorted/consistent even after filtering
    unique_configs = sorted(df['config_name'].unique(), key=natural_keys)

    for i, metric in enumerate(metrics_list):
        ax = axes[i]
        
        sns.boxplot(
            data=df, 
            x='config_name', 
            y=metric, 
            ax=ax, 
            color='skyblue', 
            showfliers=False,
            linewidth=0.4,
            order=unique_configs
        )
        
        # --- LOGIC TO GET THE "KEY" ---
        meta = INDEX_META.get(metric, metric)
        if isinstance(meta, dict):
            y_label = meta.get('key', meta.get('label', metric))
        else:
            y_label = meta

        # Set Y-axis with the short Key
        ax.set_ylabel(y_label, rotation=90, labelpad=10, 
                      fontsize=7, verticalalignment='center')
        ax.set_xlabel('')
        ax.tick_params(axis='y', labelsize=5)
        
        # --- X-AXIS LABELS ---
        ax.set_xticks(range(len(unique_configs)))
        ax.set_xticklabels(unique_configs, rotation=90, fontsize=3) # Slightly smaller for 60 cols
        ax.tick_params(axis='x', labelbottom=True)

    # Tightened hspace slightly for a cleaner landscape look
    plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.15, hspace=0.6)
    
    if output_name:
        plt.savefig(output_name, bbox_inches='tight', dpi=300)
    
    # plt.show()
    # return 

# --- Example Call ---
# plot_landscape_metrics(df, target_metrics, choose_optimal=True)