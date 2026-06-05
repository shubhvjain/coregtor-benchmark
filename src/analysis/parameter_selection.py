import  src.results.util as ut 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

def natural_keys(text):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text)]


def plot_all_configs_raw(df,metrics_list, freq_df=None, freq_metrics=[] ,output_name=None):
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
        "density_score_hippie":{"key": "PPI3 (H)", "label": "  "},
        "lcc_score_hippie":{"key": "PPI4 (H)", "label": "  "},
        "tc_score_hippie":{"key": "PPI5 (H)", "label": "  "},
        "density_score_stringdb":{"key": "PPI3 (S)", "label": "  "},
        "lcc_score_stringdb":{"key": "PPI4 (S)", "label": "  "},
        "tc_score_stringdb":{"key": "PPI5 (S)", "label": "  "},
        "density_score_biogrid":{"key": "PPI3 (B)", "label": "  "},
        "lcc_score_biogrid":{"key": "PPI4 (B)", "label": "  "},
        "tc_score_biogrid":{"key": "PPI5 (B)", "label": "  "}
    }
    n_rows = len(metrics_list) + len(freq_metrics)
    unique_configs = sorted(df["config_name"].unique(), key=natural_keys)
    n_cols = len(unique_configs)
    # --- PURE UNCONSTRAINED DYNAMIC SIZING ---
    # 6 cm converted to inches = 2.3622 inches per column
    inch_per_col = 3.0 * 0.393701
    inch_per_row = 3.0  # Extra breathing room for the vertical axis

    dynamic_width = n_cols * inch_per_col
    dynamic_height = n_rows * inch_per_row

    fig, axes = plt.subplots(
        nrows=n_rows, ncols=1, figsize=(dynamic_width, dynamic_height)
    )

    # fig, axes = plt.subplots(nrows=n_rows, ncols=1, figsize=(16.54, 11.69))

    if n_rows == 1: axes = [axes]

    # Ensure configs are sorted/consistent even after filtering
    unique_configs = sorted(df['config_name'].unique(), key=natural_keys)

    for i, metric in enumerate(metrics_list):
        ax = axes[i]
        box_color = ut.color_palette["ice_blue"]
        sns.boxplot(
            data=df, 
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

        # # Group by config_name, get the median of the current metric, and drop NaNs
        # medians = df.groupby('config_name')[metric].median().dropna()
        
        # if not medians.empty:
        #     best_config = medians.idxmax()
        #     max_median_val = medians.max()
        #     line_color = ut.color_palette["dark_navy"]
        #     # Draw the horizontal line
        #     ax.axhline(
        #         y=max_median_val, 
        #         color=line_color, 
        #         linestyle='--', 
        #         linewidth=0.5, 
        #         alpha=0.8
        #     )
            
        #     # Add text label right above the line on the far right of the plot
        #     # (Using transform=ax.get_yaxis_transform() keeps x-position relative to axes frame)
        #     ax.text(
        #         x=1.01, 
        #         y=max_median_val, 
        #         s=f"{best_config} ({max_median_val:.3f})", 
        #         color=line_color, 
        #         fontsize=6, 
        #         verticalalignment='center',
        #         transform=ax.get_yaxis_transform()
        #     )

        # --- MEDIAN AND MAX CALCULATIONS ---
        # Group by config_name to find the best median
        medians = df.groupby('config_name')[metric].median().dropna()
        # Find the absolute max across all configurations for this metric
        absolute_max_val = df[metric].max()
        
        if not medians.empty and not pd.isna(absolute_max_val):
            best_config_med = medians.idxmax()
            max_median_val = medians.max()
            
            # Find which config holds the absolute maximum value
            best_config_max = df.loc[df[metric] == absolute_max_val, 'config_name'].iloc[0]
            
            line_color = ut.color_palette["dark_navy"]
            
            # 1. Draw Max Median Line
            ax.axhline(
                y=max_median_val, 
                color=line_color, 
                linestyle='--', 
                linewidth=0.5, 
                alpha=0.6
            )
            
            # 2. Draw Absolute Max Line
            ax.axhline(
                y=absolute_max_val, 
                color=line_color, 
                linestyle=':', 
                linewidth=0.5, 
                alpha=0.6
            )
            
            # 3. Combined Text Label
            # Format: Config_Med (med_val) | Config_Max (max_val)
            label_text = (
                f"Med: {best_config_med} ({max_median_val:.3f}) \n"
                f"Max: {best_config_max} ({absolute_max_val:.3f})"
            )
            
            # Position the text midway between the two lines to prevent overlap issues
            text_y_position = (max_median_val + absolute_max_val) / 2
            
            ax.text(
                x=1.01, 
                y=text_y_position, 
                s=label_text, 
                color=line_color, 
                fontsize=5, 
                verticalalignment='center',
                transform=ax.get_yaxis_transform()
            )
        #==========
        
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
        ax.set_xticklabels(unique_configs, rotation=90, fontsize=4) # Slightly smaller for 60 cols
        ax.tick_params(axis='x', labelbottom=True)

    for i, metric in enumerate(freq_metrics):
        ax = axes[i]
        box_color = ut.color_palette["ice_blue"]
        sns.boxplot(
            data=freq_df, 
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

        # # Group by config_name, get the median of the current metric, and drop NaNs
        # medians = df.groupby('config_name')[metric].median().dropna()
        
        # if not medians.empty:
        #     best_config = medians.idxmax()
        #     max_median_val = medians.max()
        #     line_color = ut.color_palette["dark_navy"]
        #     # Draw the horizontal line
        #     ax.axhline(
        #         y=max_median_val, 
        #         color=line_color, 
        #         linestyle='--', 
        #         linewidth=0.5, 
        #         alpha=0.8
        #     )
            
        #     # Add text label right above the line on the far right of the plot
        #     # (Using transform=ax.get_yaxis_transform() keeps x-position relative to axes frame)
        #     ax.text(
        #         x=1.01, 
        #         y=max_median_val, 
        #         s=f"{best_config} ({max_median_val:.3f})", 
        #         color=line_color, 
        #         fontsize=6, 
        #         verticalalignment='center',
        #         transform=ax.get_yaxis_transform()
        #     )

        # --- MEDIAN AND MAX CALCULATIONS ---
        # Group by config_name to find the best median
        medians = freq_df.groupby('config_name')[metric].median().dropna()
        # Find the absolute max across all configurations for this metric
        absolute_max_val = freq_df[metric].max()
        
        if not medians.empty and not pd.isna(absolute_max_val):
            best_config_med = medians.idxmax()
            max_median_val = medians.max()
            
            # Find which config holds the absolute maximum value
            best_config_max = freq_df.loc[freq_df[metric] == absolute_max_val, 'config_name'].iloc[0]
            
            line_color = ut.color_palette["dark_navy"]
            
            # 1. Draw Max Median Line
            ax.axhline(
                y=max_median_val, 
                color=line_color, 
                linestyle='--', 
                linewidth=0.5, 
                alpha=0.6
            )
            
            # 2. Draw Absolute Max Line
            ax.axhline(
                y=absolute_max_val, 
                color=line_color, 
                linestyle=':', 
                linewidth=0.5, 
                alpha=0.6
            )
            
            # 3. Combined Text Label
            # Format: Config_Med (med_val) | Config_Max (max_val)
            label_text = (
                f"Med: {best_config_med} ({max_median_val:.3f}) \n"
                f"Max: {best_config_max} ({absolute_max_val:.3f})"
            )
            
            # Position the text midway between the two lines to prevent overlap issues
            text_y_position = (max_median_val + absolute_max_val) / 2
            
            ax.text(
                x=1.01, 
                y=text_y_position, 
                s=label_text, 
                color=line_color, 
                fontsize=5, 
                verticalalignment='center',
                transform=ax.get_yaxis_transform()
            )
        #==========
        
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
        ax.set_xticklabels(unique_configs, rotation=90, fontsize=4) # Slightly smaller for 60 cols
        ax.tick_params(axis='x', labelbottom=True)


    # Tightened hspace slightly for a cleaner landscape look
    plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.15, hspace=0.6)
    
    if output_name:
        plt.savefig(output_name, bbox_inches='tight', dpi=300)
    
    # plt.show()
    # return 


def generate_ps_plot1(df, freq_df, config_details, options, output_path):
    """
    Generates the integrated parameter selection analysis plot (Panels A, B, & C)
    and saves it as a publication-ready file formatted for a half A4 page.
    """
    # Set standard clean aesthetic suitable for journals
    sns.set_theme(style="whitegrid", context="paper", font_scale=0.9)
    plt.rcParams.update(ut.default_sns_configs)
    
    # Data Preprocessing & Aggregation
    working_df, cluster_labels = ut.process_and_rename_dataset(df, config_details)

    # Academic dimensions: ~8.2 inches wide (A4 width) and 3.2 inches high (half-A4 allocation)
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.2), sharey=True)
    
    # Common heatmap appearance configurations for consistent publication style
    cbar_options = {"shrink": 0.75, "pad": 0.03}
    
    # 1. First Heatmap: PPI(S) - Panel A
    pivot_b = working_df.pivot_table(
        index="distance_measure", 
        columns="clustering_method2", 
        values="PPI(S)", 
        aggfunc="mean"
    )
    sns.heatmap(pivot_b, annot=True, fmt=".2f", cmap="viridis", ax=axes[0], 
                cbar_kws=cbar_options, annot_kws={"size": 6})
    axes[0].set_title("Overall PPI Score (StringDB)", fontsize=9, pad=8)
    axes[0].set_xlabel("Clustering Method", fontsize=8)
    axes[0].set_ylabel("Distance Measure", fontsize=8)
    axes[0].tick_params(labelsize=7)
    # Relative positioning slightly above the top-left boundary of the panel
    axes[0].text(-0.25, 1.10, "(A)", transform=axes[0].transAxes, fontsize=10, va="top")
    
    # 2. Second Heatmap: GO - Panel B
    pivot_h = working_df.pivot_table(
        index="distance_measure", 
        columns="clustering_method2", 
        values="GO", 
        aggfunc="mean"
    )
    sns.heatmap(pivot_h, annot=True, fmt=".2f", cmap="viridis", ax=axes[1], 
                cbar_kws=cbar_options, annot_kws={"size": 6})
    axes[1].set_title("Overall GO Score", fontsize=9, pad=8)
    axes[1].set_xlabel("Clustering Method", fontsize=8)
    axes[1].set_ylabel("") 
    axes[1].tick_params(labelsize=7)
    axes[1].text(-0.08, 1.10, "(B)", transform=axes[1].transAxes, fontsize=10, va="top")
    
    # 3. Third Heatmap: DC - Panel C
    pivot_s = working_df.pivot_table(
        index="distance_measure", 
        columns="clustering_method2", 
        values="DC", 
        aggfunc="mean"
    )
    sns.heatmap(pivot_s, annot=True, fmt=".2f", cmap="viridis", ax=axes[2], 
                cbar_kws=cbar_options, annot_kws={"size": 6})
    axes[2].set_title("Overall Dist Corr Score", fontsize=9, pad=8)
    axes[2].set_xlabel("Clustering Method", fontsize=8)
    axes[2].set_ylabel("")
    axes[2].tick_params(labelsize=7)
    axes[2].text(-0.08, 1.10, "(C)", transform=axes[2].transAxes, fontsize=10, va="top")

    # Tight adjustments ensuring colorbars and structural sub-labels do not truncate
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


# def generate_ps_plot1(df,freq_df,config_details,options,output_path):
#     """
#     Generates the integrated parameter selection analysis plot (Panels A & B)
#     and saves it as a publication-ready SVG file.
#     """
    
#     # Set standard clean aesthetic suitable for journals
#     sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)
#     plt.rcParams.update(ut.default_sns_configs)
#     # print(df.columns)
#     # Data Preprocessing & Aggregation
#     working_df, cluster_labels = ut.process_and_rename_dataset(df,config_details)
#     #print(working_df.columns)
#     #print(working_df)
#     #print(cluster_labels)



#     fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    
#     # 1. First Heatmap: PPI(B)
#     pivot_b = working_df.pivot_table(
#         index="distance_measure", 
#         columns="clustering_method2", 
#         values="PPI(S)", 
#         aggfunc="mean"
#     )
#     sns.heatmap(pivot_b, annot=True, fmt=".2f", cmap="viridis", ax=axes[0], cbar_kws={"label": "Score"})
#     axes[0].set_title("Overall PPI Score (StringDB)")
#     axes[0].set_xlabel("Clustering Method")
#     axes[0].set_ylabel("Distance Measure")
    
#     # 2. Second Heatmap: PPI(H) (Inferred from your columns)
#     pivot_h = working_df.pivot_table(
#         index="distance_measure", 
#         columns="clustering_method2", 
#         values="GO", 
#         aggfunc="mean"
#     )
#     sns.heatmap(pivot_h, annot=True, fmt=".2f", cmap="viridis", ax=axes[1], cbar_kws={"label": "Score"})
#     axes[1].set_title("Overall GO")
#     axes[1].set_xlabel("Clustering Method")
#     axes[1].set_ylabel("") # Suppressed because sharey=True is active
    
#     # 3. Third Heatmap: PPI(S) (Inferred from your columns)
#     pivot_s = working_df.pivot_table(
#         index="distance_measure", 
#         columns="clustering_method2", 
#         values="DC", 
#         aggfunc="mean"
#     )
#     sns.heatmap(pivot_s, annot=True, fmt=".2f", cmap="viridis", ax=axes[2], cbar_kws={"label": "Score"})
#     axes[2].set_title("Overall Distance Correlation Score")
#     axes[2].set_xlabel("Clustering Method")
#     axes[2].set_ylabel("")

#     # Adjust layout to prevent any label or title truncation
#     plt.tight_layout()
#     plt.savefig(output_path, dpi=300)
#     plt.close()

#     return
    # Define baseline indices for Panel B tracking
    # master_indices = {
    #     'silhouette_score': 'Silhouette Score',
    #     'cluster_density': 'Cluster Density',
    #     'PPI_Composite': 'PPI Connectedness',
    #     'GO_Composite': 'GO Semantic Sim',
    #     'grn_collectri_jaccard': 'GRN Jaccard',
    #     'TFBS_affinity_score': 'TFBS Affinity'
    # }
    # # Keep only indices present in the dataframe columns
    # active_indices = {k: v for k, v in master_indices.items() if k in working_df.columns}

    # # Generate Mean Profiles grouped by key parameters
    # # Group by both hyperparameters for Panel A matrix mapping
    # interaction_profile = working_df.groupby(['distance_measure', 'clustering_method']).mean(numeric_only=True).reset_index()
    # # Group independently for Panel B baseline vectors
    # distance_profile = working_df.groupby('distance_measure').mean(numeric_only=True)
    # method_profile = working_df.groupby('clustering_method').mean(numeric_only=True)

    # # 3. Structural Canvas Construction (GridSpec Assembly)
    # fig = plt.figure(figsize=(15, 11))
    # # 2 Main rows, 6 horizontal grid divisions to easily share spans (2+2+2 for Top, 3+3 for Bottom)
    # gs = fig.add_gridspec(2, 6, height_ratios=[1, 1], hspace=0.35, wspace=0.45)

    # # Define unique ordering to keep the vertical layout uniform across axes
    # dist_order = sorted(working_df['distance_measure'].unique())
    # method_order = sorted(working_df['clustering_method'].unique())

    # # ==========================================
    # # PANEL A: Targeted Hyperparameter Matrices
    # # ==========================================
    # ax_a1 = fig.add_subplot(gs[0, 0:2])
    # ax_a2 = fig.add_subplot(gs[0, 2:4], sharey=ax_a1)
    # ax_a3 = fig.add_subplot(gs[0, 4:6], sharey=ax_a1)

    # # A1: PPI Layer
    # pivot_ppi = interaction_profile.pivot(index='distance_measure', columns='clustering_method', values='PPI_Composite').reindex(index=dist_order, columns=method_order)
    # sns.heatmap(pivot_ppi, annot=True, fmt=".3f", cmap="YlGnBu", ax=ax_a1, cbar_kws={'label': 'Mean Composite Score'})
    # ax_a1.set_title('A1: PPI Validation Layer\n(Hippie + StringDB + BioGrid)', fontsize=10, fontweight='bold')
    # ax_a1.set_ylabel('Distance Measure', fontsize=10)

    # # A2: Gene Ontology Layer
    # pivot_go = interaction_profile.pivot(index='distance_measure', columns='clustering_method', values='GO_Composite').reindex(index=dist_order, columns=method_order)
    # sns.heatmap(pivot_go, annot=True, fmt=".3f", cmap="Blues", ax=ax_a2, cbar_kws={'label': 'Mean Functional Similarity'})
    # ax_a2.set_title('A2: Gene Ontology Layer\n(Lin + Resnik)', fontsize=10, fontweight='bold')
    # ax_a2.set_xlabel('Clustering Method', fontsize=10, labelpad=8)
    # ax_a2.tick_params(left=False, labelleft=False) # Hide redundant shared labels

    # # A3: TFBS Layer
    # if 'TFBS_affinity_score' in working_df.columns:
    #     pivot_tfbs = interaction_profile.pivot(index='distance_measure', columns='clustering_method', values='TFBS_affinity_score').reindex(index=dist_order, columns=method_order)
    #     sns.heatmap(pivot_tfbs, annot=True, fmt=".3f", cmap="Oranges", ax=ax_a3, cbar_kws={'label': 'TFBS Affinity Score'})
    #     ax_a3.set_title('A3: Transcription Factor Layer\n(TFBS Affinity Only)', fontsize=10, fontweight='bold')
    #     ax_a3.tick_params(left=False, labelleft=False)

    # # Inject Panel A Heading Context
    # fig.text(0.02, 0.95, "PANEL A: Targeted Hyperparameter Interactions Across Functional Layers", fontsize=13, fontweight='bold', color='#1a1a1a')

    # # ==========================================
    # # PANEL B: Isolated Profile Vector Baselines
    # # ==========================================
    # ax_b1 = fig.add_subplot(gs[1, 0:3])
    # ax_b2 = fig.add_subplot(gs[1, 3:6])

    # # Extract ordered matrix columns based on active configuration variables
    # metric_keys = list(active_indices.keys())
    # metric_labels = list(active_indices.values())

    # # B1: Distance Profile vs Metric Matrix
    # mat_dist = distance_profile[metric_keys].reindex(dist_order)
    # sns.heatmap(mat_dist, annot=True, fmt=".3f", cmap="BuPu", ax=ax_b1, xticklabels=metric_labels, cbar_kws={'label': 'Metric Baseline Value'})
    # ax_b1.set_title('B1: Distance Measures vs All Evaluation Indices', fontsize=10, fontweight='bold')
    # ax_b1.set_ylabel('Distance Measure', fontsize=10)
    # ax_b1.tick_params(axis='x', rotation=25, labelsize=9)

    # # B2: Clustering Method Profile vs Metric Matrix
    # mat_method = method_profile[metric_keys].reindex(method_order)
    # sns.heatmap(mat_method, annot=True, fmt=".3f", cmap="BuPu", ax=ax_b2, xticklabels=metric_labels, cbar_kws={'label': 'Metric Baseline Value'})
    # ax_b2.set_title('B2: Clustering Methods vs All Evaluation Indices', fontsize=10, fontweight='bold')
    # ax_b2.set_ylabel('Clustering Method', fontsize=10)
    # ax_b2.tick_params(axis='x', rotation=25, labelsize=9)

    # # Inject Panel B Heading Context
    # fig.text(0.02, 0.47, "PANEL B: Isolated Performance Profiling Against Master Index Suite", fontsize=13, fontweight='bold', color='#1a1a1a')

    # # 4. Save Vector Blueprint Block
    # plt.savefig(output_path, dpi=300, bbox_inches='tight', format='svg')
    # plt.close(fig)
    # # plt.show()
    # print(f"Publication figure successfully exported to -> {output_path}")