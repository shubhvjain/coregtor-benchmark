import  src.results.util as ut 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

def natural_keys(text):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text)]

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_all_configs_raw(df, freq_df=None, output_name=None):
    """
    Plots a multi-row boxplot comparison for configs.
    """
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
        "density_score_hippie": {"key": "PPI3 (H)", "label": "Density Score (Hippie)"},
        "lcc_score_hippie": {"key": "PPI4 (H)", "label": "LCC Score (Hippie)"},
        "tc_score_hippie": {"key": "PPI5 (H)", "label": "TC Score (Hippie)"},
        "density_score_stringdb": {"key": "PPI3 (S)", "label": "Density Score (StringDB)"},
        "lcc_score_stringdb": {"key": "PPI4 (S)", "label": "LCC Score (StringDB)"},
        "tc_score_stringdb": {"key": "PPI5 (S)", "label": "TC Score (StringDB)"},
        "density_score_biogrid": {"key": "PPI3 (B)", "label": "Density Score (Biogrid)"},
        "lcc_score_biogrid": {"key": "PPI4 (B)", "label": "LCC Score (Biogrid)"},
        "tc_score_biogrid": {"key": "PPI5 (B)", "label": "TC Score (Biogrid)"},
        "dCor_targets_score":{"key":"CT DC Score", "label":""},
        "ppi_shared_partners_stringdb":{"key":"CT PPI Score (S)", "label":""},
        "ppi_shared_partners_biogrid":{"key":"CT PPI Score (B)", "label":""},
        "ppi_shared_partners_hippie":{"key":"CT PPI Score (H)", "label":""},
        "go_score":{"key":"CT GO Score", "label":""}

    }
    
    # Internal hardcoded lists provided by user
    metrics_list = [
        'silhouette_score', 'n_source', 'cluster_density', 'cluster_diameter', 
        'shortest_PPI_path_score_hippie', 'shortest_PPI_path_score_stringdb',
        'shortest_PPI_path_score_biogrid', 'shared_PPI_partners_score_hippie',
        'shared_PPI_partners_score_stringdb', 'shared_PPI_partners_score_biogrid', 
        'goa_similarity_lin', 'goa_similarity_resnik', 'goa_similarity_jc', 
        'density_hippie', 'density_score_hippie', 'lcc_hippie', 'lcc_score_hippie', 
        'tc_hippie', 'tc_score_hippie', 'node_found_ratio_hippie', 'density_stringdb',
        'density_score_stringdb', 'lcc_stringdb', 'lcc_score_stringdb',
        'tc_stringdb', 'tc_score_stringdb', 'node_found_ratio_stringdb',
        'density_biogrid', 'density_score_biogrid', 'lcc_biogrid',
        'lcc_score_biogrid', 'tc_biogrid', 'tc_score_biogrid',
        'node_found_ratio_biogrid', 'dCor_sources', 'dCor_target', 'DC1', 'DC2',
        'TFBS_affinity', 'TF_found_per', 'TFBS_affinity_score'
    ]
    
    freq_metrics = [
        'dCor_targets_score', 'ppi_shared_partners_hippie', 
        'ppi_shared_partners_biogrid','ppi_shared_partners_stringdb', 'go_score'
    ]
    
    # Define tasks to process sequentially
    # Format: (associated_dataframe, metric_name)
    tasks = []
    for m in metrics_list:
        if m in df.columns:  # Safety check to avoid KeyErrors if df missing a metric
            tasks.append((df, m))
            
    for m in freq_metrics:
        if freq_df is not None and m in freq_df.columns:
            tasks.append((freq_df, m))
            
    n_rows = len(tasks)
    if n_rows == 0:
        print("No valid metrics found in dataframes.")
        return

    unique_configs = sorted(df["config_name"].unique(), key=natural_keys)
    n_cols = len(unique_configs)
    
    # --- PURE UNCONSTRAINED DYNAMIC SIZING ---
    inch_per_col = 3.0 * 0.393701
    inch_per_row = 3.5  # Room for Max, Med, and Min line info labels

    dynamic_width = n_cols * inch_per_col
    dynamic_height = n_rows * inch_per_row

    fig, axes = plt.subplots(
        nrows=n_rows, ncols=1, figsize=(dynamic_width, dynamic_height)
    )

    if n_rows == 1: 
        axes = [axes]

    # --- SINGLE UNIFIED LOOP ---
    for idx, (current_df, metric) in enumerate(tasks):
        ax = axes[idx]
        box_color = ut.color_palette["ice_blue"]
        
        sns.boxplot(
            data=current_df, 
            x='config_name', 
            y=metric, 
            ax=ax, 
            color=box_color, 
            showfliers=False,
            linewidth=0.4,
            order=unique_configs,
            width=0.75
        )

        ax.grid(True, axis='y', linewidth=0.4, color='#e0e0e0', zorder=0)

        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

        # --- MIN, MEDIAN, AND MAX CALCULATIONS ---
        medians = current_df.groupby('config_name')[metric].median().dropna()
        absolute_max_val = current_df[metric].max()
        absolute_min_val = current_df[metric].min()
        
        if not medians.empty and not pd.isna(absolute_max_val) and not pd.isna(absolute_min_val):
            best_config_med = medians.idxmax()
            max_median_val = medians.max()
            
            best_config_max = current_df.loc[current_df[metric] == absolute_max_val, 'config_name'].iloc[0]
            best_config_min = current_df.loc[current_df[metric] == absolute_min_val, 'config_name'].iloc[0]
            
            line_color = ut.color_palette["dark_navy"]
            
            # 1. Max Median Line (dashed)
            ax.axhline(y=max_median_val, color=line_color, linestyle='--', linewidth=0.5, alpha=0.6)
            
            # 2. Absolute Max Line (dotted)
            ax.axhline(y=absolute_max_val, color=line_color, linestyle=':', linewidth=0.5, alpha=0.6)

            # 3. Absolute Min Line (dash-dot)
            ax.axhline(y=absolute_min_val, color=line_color, linestyle='-.', linewidth=0.5, alpha=0.6)
            
            # 4. Text Label Block
            label_text = (
                f"Max: {best_config_max} ({absolute_max_val:.3f})\n"
                f"Med: {best_config_med} ({max_median_val:.3f})\n"
                f"Min: {best_config_min} ({absolute_min_val:.3f})"
            )
            
            text_y_position = (absolute_max_val + absolute_min_val) / 2
            
            ax.text(
                x=1.01, 
                y=text_y_position, 
                s=label_text, 
                color=line_color, 
                fontsize=5, 
                verticalalignment='center',
                transform=ax.get_yaxis_transform()
            )
        
        # --- LOGIC TO GET THE "KEY" ---
        meta = INDEX_META.get(metric, metric)
        if isinstance(meta, dict):
            y_label = meta.get('key', meta.get('label', metric))
        else:
            # Fallback formatting if metric isn't in INDEX_META dictionary
            y_label = str(meta).replace('_', ' ').title()

        ax.set_ylabel(y_label, rotation=90, labelpad=10, 
                      fontsize=7, verticalalignment='center')
        ax.set_xlabel('')
        ax.tick_params(axis='y', labelsize=5)
        
        # --- X-AXIS LABELS ---
        ax.set_xticks(range(len(unique_configs)))
        ax.set_xticklabels(unique_configs, rotation=90, fontsize=4) 
        ax.tick_params(axis='x', labelbottom=True)

    plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.15, hspace=0.6)
    
    if output_name:
        plt.savefig(output_name, bbox_inches='tight', dpi=300)


def generate_heatmap(ax, pivot_df, title, letter, show_y_label=True):
    """
    Generates a single heatmap panel on a provided matplotlib axis using a pre-pivoted dataframe.
    """
    # Skip rendering if pivot data is missing or empty
    if pivot_df is None or pivot_df.empty:
        ax.axis('off')
        return

    # 1. Render Heatmap directly from the passed matrix
    cbar_options = {"shrink": 0.75, "pad": 0.03}
    annot_options = {"size": 4}
    
    sns.heatmap(pivot_df, annot=True, fmt=".2f", cmap="viridis", ax=ax, 
                cbar_kws=cbar_options, annot_kws=annot_options)

    # 2. Apply Local Typographic Styles
    ax.set_title(title)
    ax.set_xlabel("Clustering Method")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # 3. Handle Left-Column Specific Formatting Independently
    if show_y_label:
        ax.set_ylabel("Distance Measure")
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])  # Hide tick labels safely
    
    panel_label_x = -0.25

    # 4. Add Panel Identifier Label
    panel_label_y = 1.15 
    ax.text(panel_label_x, panel_label_y, letter, transform=ax.transAxes, 
            fontsize=10, va="top")



def generate_ps_plot1(df, freq_df, config_details, options, output_path):
    """
    Generates an asymmetric grid layout (Row 1: 3 equal cols | Row 2: 2 cols at 25%/75% width)
    by calculating matrices centrally and passing them to independent component calls.
    """
    # Set standard clean aesthetic suitable for journals
    sns.set_theme(style="whitegrid", context="paper", font_scale=0.9)
    plt.rcParams.update(ut.default_sns_configs)
    
    # Data Preprocessing & Aggregation
    working_df, cluster_labels = ut.process_and_rename_dataset(df, config_details)

    working_freq_df , freq_cluster_labels = ut.process_and_rename_freq(freq_df, config_details)

    # Establish Figure Container
    fig = plt.figure(figsize=(8.2, 6.5))
    
    # Define a 2-row, 12-column Grid Spec allocation matrix
    gs = fig.add_gridspec(nrows=3, ncols=12)
    
    # -------------------------------------------------------------
    # ROW 1 PLOTS: 3 Equal Columns (Each spans 4 grid units)
    # -------------------------------------------------------------
    # Panel A: PPI(S)
    pivot_a = working_df.pivot_table(
        columns ="distance_measure", index="clustering_method2", values="PPI(S)", aggfunc="mean"
    ) if "PPI(S)" in working_df.columns else None
    
    ax_a = fig.add_subplot(gs[0, 0:4])   
    generate_heatmap(
        ax=ax_a, pivot_df=pivot_a, title="Overall PPI Score (StringDB)", letter="(A)", show_y_label=True
    )
    
    # Panel B: GO
    pivot_b = working_df.pivot_table(
        columns ="distance_measure", index="clustering_method2", values="GO", aggfunc="mean"
    ) if "GO" in working_df.columns else None
    
    ax_b = fig.add_subplot(gs[0, 4:8])   
    generate_heatmap(
        ax=ax_b, pivot_df=pivot_b, title="Overall GO Score", letter="(B)", show_y_label=False
    )
    
    # Panel C: DC
    pivot_c = working_df.pivot_table(
        columns ="distance_measure", index="clustering_method2", values="DC", aggfunc="mean"
    ) if "DC" in working_df.columns else None
    
    ax_c = fig.add_subplot(gs[0, 8:12])  
    generate_heatmap(
        ax=ax_c, pivot_df=pivot_c, title="Overall Dist Corr Score", letter="(C)", show_y_label=False
    )
    
    # -------------------------------------------------------------
    # ROW 2 PLOTS: Asymmetric Columns (25% vs 75% width split)
    # -------------------------------------------------------------
    # Panel D: PPI(H) - Left Panel: 25% Width (Spans 3 units: index 0 to 2)
    pivot_d = working_df.pivot_table(
        columns ="distance_measure", index ="clustering_method2", values="TFBS", aggfunc="mean"
    ) if "TFBS" in working_df.columns else None
    
    ax_d = fig.add_subplot(gs[1, 0:4])   
    generate_heatmap(
        ax=ax_d, pivot_df=pivot_d, title="TFBS Affinity Score", letter="(D)", show_y_label=True
    )
    
    # Panel E: sil score 

    pivot_e = working_df.pivot_table(
       columns ="distance_measure", index="clustering_method2", values="silhouette_score", aggfunc="mean"
    ).fillna(0) if working_df is not None and "silhouette_score" in working_df.columns else None
    # pivot_e = None

    ax_e = fig.add_subplot(gs[1, 4:8])  
    generate_heatmap(
        ax=ax_e, pivot_df=pivot_e, title="Silhouette Score", letter="(E)", show_y_label=False
    )

    # Panel F: go_score (from freq_df) - Right Panel: 75% Width (Spans 9 units: index 3 to 11)

    pivot_f = working_freq_df.pivot_table(
       columns ="distance_measure", index="clustering_method2", values="CTS", aggfunc="mean"
    ).fillna(0) if working_freq_df is not None and "CTS" in working_freq_df.columns else None
    # pivot_e = None

    ax_e = fig.add_subplot(gs[1, 8:12])  
    generate_heatmap(
        ax=ax_e, pivot_df=pivot_f, title="Overall Common Targets Score", letter="(F)", show_y_label=False
    )

    
    # # Panel F: database vs overall score

    # pivot_f = working_df.pivot_table(
    #    columns ="dataset", index="clustering_method2", values="OVERALL_SCORE", aggfunc="mean"
    # ).fillna(0) if working_df is not None and "OVERALL_SCORE" in working_df.columns else None
    # # pivot_e = None

    # ax_f = fig.add_subplot(gs[1, 7:12])  
    # generate_heatmap(
    #     ax=ax_f, pivot_df=pivot_f, title="Overall Score", letter="(F)", show_y_label=False
    # )

    # Final spatial compilation configurations
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()