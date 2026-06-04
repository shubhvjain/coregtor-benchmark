"""
"""
import src.results.util as ut
from src.results.parameter_selection import plot_landscape_metrics
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)


def _get_exp_stats(exp, output_path, temp_path):
    """"""
    dataset_info = ut.get_exp_gtex_info(exp["dataset"]["name"])

    parquet_path = output_path / "target_stats.parquet"
    json_path = output_path / "target_stats.json"

    # Check for Parquet first, then fall back to JSON
    if parquet_path.exists():
        data = pd.read_parquet(parquet_path)
    elif json_path.exists():
        with open(json_path) as f:
            target_file = json.load(f)
        data = pd.json_normalize(data=target_file)
    else:
        raise ValueError(
            f"Neither {parquet_path.name} nor {json_path.name} exists"
        )

    # Add the metadata columns
    data["dataset_name"] = dataset_info["name"]
    data["dataset_abbr"] = dataset_info["abbr"]
    data["dataset_sample_size"] = dataset_info["n_sample"]
    data["exp_id"] = exp["id"]

    return data

def _generate_save_stats_diagram(df, folder_path,options):
    # 'paper' context ensures labels and fonts are appropriately scaled for publication
    # print(db_keys)

    # Fits safely within half of an A4 page (width: ~6.5", height: 4.8")
    HALF_A4_GRID = (6.5, 4.8)

    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    # Removed sharex=True to allow manual control of x-ticks per subplot
    fig, axes = plt.subplots(
        nrows=3,
        ncols=2,
        figsize=HALF_A4_GRID,
        # Slightly increased hspace for rotated labels
        gridspec_kw={"wspace": 0.35, "hspace": 0.55}
    )

    # Flatten the 2D array of axes to a 1D list
    ax_list = axes.flatten()

    default_bar_color = ut.color_palette["bright_yellow"]
    default_box_color = ut.color_palette["pale_blue"]
    default_box_width = 0.45
    default_line_width = 0.5

    # --- 1. Samples Barplot (Top Left) ---
    sns.barplot(data=df, x="dataset_abbr", y="dataset_sample_size",
                ax=ax_list[0], color=default_bar_color)
    ax_list[0].set_ylabel("# samples")
    if ax_list[0].containers:
        ax_list[0].bar_label(ax_list[0].containers[0], padding=2, fontsize=7)

    # --- 2. Root Nodes (Top Right) ---
    sns.boxplot(data=df, x="dataset_abbr", y="gf.total_root_nodes", ax=ax_list[1],
                color=default_box_color, width=default_box_width, showfliers=False, linewidth=default_line_width)
    ax_list[1].set_ylabel("# root nodes")

    # --- 3. Context Density (Middle Left) ---
    sns.boxplot(data=df, x="dataset_abbr", y="gf.mean_context_density", ax=ax_list[2],
                color=default_box_color, width=default_box_width, showfliers=False, linewidth=default_line_width)
    ax_list[2].set_ylabel("mean context density")

    # --- 4. Evenness Index (Middle Right) ---
    sns.boxplot(data=df, x="dataset_abbr", y="gf.mean_evenness_index", ax=ax_list[3],
                color=default_box_color, width=default_box_width, showfliers=False, linewidth=default_line_width)
    ax_list[3].set_ylabel("mean evenness index")

    # --- 5. Peer Density (Bottom Left) ---
    sns.boxplot(data=df, x="dataset_abbr", y="gf.mean_peer_density", ax=ax_list[4],
                color=default_box_color, width=default_box_width, showfliers=False, linewidth=default_line_width)
    ax_list[4].set_ylabel("mean peer density")

    # --- 6. Peer Rank (Bottom Right) ---
    sns.boxplot(data=df, x="dataset_abbr", y="gf.mean_peer_rank", ax=ax_list[5],
                color=default_box_color, width=default_box_width, showfliers=False, linewidth=default_line_width)
    ax_list[5].set_ylabel("mean peer rank")

    # Toggle to show or hide x-tick labels for upper rows
    show_all_x_ticks = True

    # Layout tuning loop
    for i, ax in enumerate(ax_list):
        # Determine whether to show x-tick labels for this subplot
        should_show_label = show_all_x_ticks or (i >= 4)

        # Format X-axis ticks
        ax.tick_params(axis="x", labelbottom=should_show_label,
                       rotation=45, labelsize=7)

        # Format Y-axis ticks
        ax.tick_params(axis="y", labelsize=7)

        ax.yaxis.label.set_size(7)

        # chr(65) is 'A', chr(66) is 'B', etc.
        subplot_letter = f"({chr(65 + i)})"
        
        # Positioned slightly above and to the left of the actual axes boundary
        ax.text(-0.22, 1.05, subplot_letter, 
                transform=ax.transAxes, 
                fontsize=7, 
                va="bottom", 
                ha="right")

        # Only show the "Dataset" x-axis title on the bottom row to reduce clutter
        if i >= 4:
            ax.set_xlabel("Dataset", fontsize=7)
        else:
            ax.set_xlabel("")

    # --- GENERATE FOOTER CAPTION FOR DATASETS ---
    db_keys = ut.get_gtex_key()
    
    show_label = options.get("show_dataset_legend",True)
    if show_label:
        # Get sorted unique abbreviations present in this specific dataframe
        unique_abbrs = list(df["dataset_abbr"].dropna().unique())
        
        # Map them to "ABBR - Long Name", defaulting to the abbreviation if missing from the dict
        legend_items = [f"{abbr} - {db_keys.get(abbr, abbr)}" for abbr in unique_abbrs]
        
        # Split the items into lines (e.g., chunks of up to 4 items per line depending on width)
        chunk_size = 6
        lines = []
        for i in range(0, len(legend_items), chunk_size):
            chunk = legend_items[i:i + chunk_size]
            lines.append(" | ".join(chunk))
        
        legend_text = "\n".join(lines)

        # Push the subplots up slightly to make room for the caption text at the bottom
        plt.subplots_adjust(bottom=0.20)

        # Place the multi-line footer text across the full width of the figure
        fig.text(
            x=0.05,            # Starts near the left margin
            y=0.02,            # Placed near the absolute bottom boundary
            s=legend_text, 
            fontsize=6,        # Small font size to comfortably fit ~5 lines
            color="#333333",   # Charcoal gray to reduce visual clutter
            va="bottom", 
            ha="left",
            wrap=True          # Safe fallback if lines get long
        )

    # Save and prevent memory leaks
    fig.savefig(folder_path / "gene_frequency_histogram.pdf",
                bbox_inches="tight")
    fig.savefig(folder_path / "gene_frequency_histogram.svg",
                bbox_inches="tight")
    plt.close(fig)


def _generate_save_runtime_diagram(df, folder_path):
    # normalize runtime for 8 core cpu
    df["runtime8core"] = round( (df["stats.total_time"] * df["stats.n_cores_used"])/ 8  , 6 )
    
    
    # 'paper' context ensures labels and fonts are appropriately scaled for publication
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)

    # Fits safely within half of an A4 page (width: ~6.5", height: 4.8")
    HALF_A4_GRID = (6.5, 4.8)

    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    # Removed sharex=True to allow manual control of x-ticks per subplot
    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=HALF_A4_GRID,
        # Slightly increased hspace for rotated labels
        gridspec_kw={"wspace": 0.35, "hspace": 0.5}
    )

    # Flatten the 2D array of axes to a 1D list
    ax_list = axes.flatten()

    default_bar_color = ut.color_palette["bright_yellow"]
    default_box_color = ut.color_palette["ice_blue"]
    default_box_width = 0.45
    default_line_width = 0.5

    # --- 1. Samples Barplot (Top Left) ---
    sns.barplot(data=df, x="dataset_abbr", y="dataset_sample_size",
                ax=ax_list[0], color=default_bar_color)
    ax_list[0].set_ylabel("# samples")
    if ax_list[0].containers:
        ax_list[0].bar_label(ax_list[0].containers[0], padding=2, fontsize=7)

    # --- 2. Root Nodes (Top Right) ---
    sns.boxplot(data=df, x="dataset_abbr", y="runtime8core", ax=ax_list[1],
                color=default_box_color, width=default_box_width, showfliers=False, linewidth=default_line_width)
    ax_list[1].set_ylabel("Runtime (seconds)")

    # Toggle to show or hide x-tick labels for upper rows
    show_all_x_ticks = True

    # Layout tuning loop
    for i, ax in enumerate(ax_list):
        # Determine whether to show x-tick labels for this subplot
        should_show_label = show_all_x_ticks or (i >= 4)

        # Format X-axis ticks
        ax.tick_params(axis="x", labelbottom=should_show_label,
                       rotation=45, labelsize=7)

        # Format Y-axis ticks
        ax.tick_params(axis="y", labelsize=9)

        # Only show the "Dataset" x-axis title on the bottom row to reduce clutter
        if i >= 4:
            ax.set_xlabel("Dataset", fontsize=9)
        else:
            ax.set_xlabel("")

    # Save and prevent memory leaks
    fig.savefig(folder_path / "runtime.pdf", bbox_inches="tight")
    plt.close(fig)

def combined_target_diagram(input, env, options, args):
    """
    """
    apath = ut.get_analysis_path(env)
    folder = apath/input["id"]
    folder.mkdir(exist_ok=True, parents=True)
    rerun = args.rerun

    # print("hello")

    # out_file = folder / f"gene_freq_hist.svg"

    # if out_file.exists() and not rerun:
    #     print("File already exists")
    #     return

    data = []
    for e in input["experiment_results_paths"]:
        o, t = ut.get_exp_path({"path": e}, env)
        with open(o/"input.json") as f:
            input_data = json.load(f)
        details = _get_exp_stats(input_data, o, t)
        data.append(details)

    df = pd.concat(data, ignore_index=True)
    save_file =  options.get("saved_combined_stats_file",False)
    if save_file: 
        df.to_csv(folder/"target_stats.csv", index=False)
    sort_col = options.get("sort_by", "dataset_name")
    df.sort_values(by=sort_col, inplace=True)
    _generate_save_stats_diagram(df, folder,options)
    _generate_save_runtime_diagram(df, folder)


def parameter_config_plot(input, env, options, args):
    """"""
    apath = ut.get_analysis_path(env)
    folder = apath/input["id"]
    folder.mkdir(exist_ok=True, parents=True)
    rerun = args.rerun

    out_file = folder / f"param_selection_{options['result_name']}.svg"

    if out_file.exists() and not rerun:
        print("File already exists")
        return

    data = []
    all_exp = options.get("experiments", [])
    exp_list = options.get("exp_list", all_exp)
    for e in exp_list:
        o, t = ut.get_exp_path({"id": e}, env)

        with open(o/"input.json") as f:
            input_data = json.load(f)
        details = pd.read_csv(o/f"{options['result_name']}_indices.csv.gz")
        details["exp_id"] = e
        dataset_info = ut.get_exp_gtex_info(input_data["dataset"]["name"])
        details["dataset_name"] = dataset_info["name"]
        details["dataset_abbr"] = dataset_info["abbr"]
        data.append(details)

    df = pd.concat(data, ignore_index=True)
    sf = options.get("save_file", False)
    if sf:
        df.to_csv(
            folder/f"clusters_{options['result_name']}.csv.gz", index=False)

    default_metrics = ['silhouette_score',
                       'n_source',
                       'shared_PPI_partners_score_stringdb',
                       'shortest_PPI_path_score_stringdb',
                       'density_score_stringdb',
                       'lcc_score_stringdb',
                       'tc_score_stringdb',
                       'goa_similarity_lin',
                       'goa_similarity_resnik',
                       'goa_similarity_jc'
                       ]

    metrics_list = options.get("metrics_list", default_metrics)

    pconfig =  ut.get_ps_configs()
    plot_landscape_metrics(df, metrics_list, out_file, False)
    result_plot1 = folder / f"param_selection_summary_{options['result_name']}.pdf"
    plot_pe_results(df,pconfig,result_plot1)
    result_plot2 = folder / f"param_selection_ppi_{options['result_name']}.pdf"
    plot_ps_ppi_grid(df,pconfig,result_plot2)



def clean_algo_names(row):
    method = row['clustering_method']
    opt1 = str(row['clustering_option1'])
    
    if 'hierarchical' in method:
        return 'Hierarchical (Inc)' if 'inconsistency' in opt1 else 'Hierarchical (OC)'
    elif 'hdbscan' in method:
        return 'HDBSCAN'
    elif 'community_detection' in method:
        if 'percentile_50' in opt1: return 'CD (P-50)'
        if 'percentile_75' in opt1: return 'CD (P-75)'
        if 'knn' in opt1: return 'CD (KNN)'
    return method


def plot_pe_results1(df, config, output_name=None):
    # Merge experimental results with the config dataframe
    merged_df = df.merge(config, left_on='config_name', right_on='config_Id', how='left')
    
    # Process categorical labels for plotting
    merged_df['Clustering Method'] = merged_df.apply(clean_algo_names, axis=1)
    
    # Standardize names to match your manual list exactly
    merged_df['Distance Measure'] = merged_df['distance_measure'].str.replace('_distance', '').str.title()
    name_map = {
        'Euclidian': 'Euclidean',
        'Euclidean': 'Euclidean',
        'Chebysehv': 'Chebyshev',
        'Chebyshev': 'Chebyshev',
        'Weighted_Jaccard': 'Weight Jaccard',
        'Jensenshannon': 'Jensen-Shannon',
        'Sorensen': 'Sørensen'
    }
    merged_df['Distance Measure'] = merged_df['Distance Measure'].replace(name_map)

    # Calculate average PPI score from the five requested stringdb metrics
    ppi_cols = [
        'shared_PPI_partners_score_stringdb', 
        'shortest_PPI_path_score_stringdb',
        'density_score_stringdb', 
        'lcc_score_stringdb', 
        'tc_score_stringdb'
    ]
    merged_df['Avg PPI Score'] = merged_df[ppi_cols].mean(axis=1)

    # Initialize the 2x2 subplot matrix
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))
    
    # Set plot style configurations using the utility palette
    plt.rcParams['font.family'] = 'sans-serif'
    grid_color = ut.color_palette["light_gray"]
    border_color = ut.color_palette["medium_gray"]
    
    # Your requested exact manual sort order
    requested_dist_order = [
        'Manhattan', 'Canberra', 'Euclidean', 'Chebyshev', 
        'Cosine', 'Sørensen', 'Weight Jaccard', 'Jensen-Shannon'
    ]
    existing_dists = merged_df['Distance Measure'].unique()
    dist_order = [d for d in requested_dist_order if d in existing_dists]
    
    algo_order = ['CD (P-50)', 'CD (P-75)', 'CD (KNN)', 'Hierarchical (Inc)', 'Hierarchical (OC)', 'HDBSCAN']

    # Row 1 Column 1: Distance Measures vs Silhouette Score
    ax = axes[0, 0]
    sns.boxplot(data=merged_df, x='Distance Measure', y='silhouette_score', ax=ax,
                order=dist_order, color=ut.color_palette["classic_blue"], width=0.5, showfliers=False, linewidth=0.8)
    
    medians = merged_df.groupby('Distance Measure')['silhouette_score'].median()
    if not medians.empty:
        max_val, best_cat = medians.max(), medians.idxmax()
        ax.axhline(y=max_val, color=ut.color_palette["vibrant_orange"], linestyle='--', linewidth=1, zorder=2)
        #ax.text(1.02, max_val, f"Max Med:\n{best_cat}\n({max_val:.3f})", color=ut.color_palette["vibrant_orange"], fontsize=8, va='center', transform=ax.get_yaxis_transform(), fontweight='bold')

    ax.set_ylabel('Silhouette Score', fontsize=10,  color=ut.color_palette["dark_navy"])
    ax.set_title('1A. Silhouette Performance by Distance Metric', fontsize=11, loc='left', color=ut.color_palette["midnight_blue"], pad=10)

    # Row 1 Column 2: Clustering Algorithm vs Silhouette Score
    ax = axes[0, 1]
    sns.boxplot(data=merged_df, x='Clustering Method', y='silhouette_score', ax=ax,
                order=algo_order, color=ut.color_palette["ocean_blue"], width=0.5, showfliers=False, linewidth=0.8)
    
    medians = merged_df.groupby('Clustering Method')['silhouette_score'].median()
    if not medians.empty:
        max_val, best_cat = medians.max(), medians.idxmax()
        ax.axhline(y=max_val, color=ut.color_palette["vibrant_orange"], linestyle='--', linewidth=1, zorder=2)
        #ax.text(1.02, max_val, f"Max Med:\n{best_cat}\n({max_val:.3f})", color=ut.color_palette["vibrant_orange"], fontsize=8, va='center', transform=ax.get_yaxis_transform(), fontweight='bold')

    ax.set_ylabel('Silhouette Score', fontsize=10,  color=ut.color_palette["dark_navy"])
    ax.set_title('1B. Silhouette Performance by Clustering Variant', fontsize=11, loc='left', color=ut.color_palette["midnight_blue"], pad=10)

    # Row 2 Column 1: Distance Measures vs Average PPI Score
    ax = axes[1, 0]
    sns.boxplot(data=merged_df, x='Distance Measure', y='Avg PPI Score', ax=ax,
                order=dist_order, color=ut.color_palette["vibrant_orange"], width=0.5, showfliers=False, linewidth=0.8)
    
    medians = merged_df.groupby('Distance Measure')['Avg PPI Score'].median()
    if not medians.empty:
        max_val, best_cat = medians.max(), medians.idxmax()
        ax.axhline(y=max_val, color=ut.color_palette["bright_yellow"], linestyle='--', linewidth=1, zorder=2)
        #ax.text(1.02, max_val, f"Max Med:\n{best_cat}\n({max_val:.3f})", color=ut.color_palette["vibrant_orange"], fontsize=8, va='center', transform=ax.get_yaxis_transform(), fontweight='bold')

    ax.set_ylabel('Avg PPI Score (StringDB)', fontsize=10,  color=ut.color_palette["dark_navy"])
    ax.set_title('2A. Biological PPI Coherence by Distance Metric', fontsize=11, loc='left', color=ut.color_palette["midnight_blue"], pad=10)

    # Row 2 Column 2: Clustering Algorithm vs Average PPI Score
    ax = axes[1, 1]
    sns.boxplot(data=merged_df, x='Clustering Method', y='Avg PPI Score', ax=ax,
                order=algo_order, color=ut.color_palette["bright_yellow"], width=0.5, showfliers=False, linewidth=0.8)
    
    medians = merged_df.groupby('Clustering Method')['Avg PPI Score'].median()
    if not medians.empty:
        max_val, best_cat = medians.max(), medians.idxmax()
        ax.axhline(y=max_val, color=ut.color_palette["vibrant_orange"], linestyle='--', linewidth=1, zorder=2)
        # ax.text(1.02, max_val, f"Max Med:\n{best_cat}\n({max_val:.3f})", color=ut.color_palette["vibrant_orange"],  fontsize=8, va='center', transform=ax.get_yaxis_transform(), fontweight='bold')

    ax.set_ylabel('Avg PPI Score (StringDB)', fontsize=10, color=ut.color_palette["dark_navy"])
    ax.set_title('2B. Biological PPI Coherence by Clustering Variant', fontsize=11, loc='left', color=ut.color_palette["midnight_blue"], pad=10)

    # Apply structural cleanup across all subplots
    for row in axes:
        for subplot_ax in row:
            subplot_ax.set_xlabel('')
            subplot_ax.tick_params(axis='both', labelsize=9, colors=ut.color_palette["dark_navy"])
            subplot_ax.set_xticklabels(subplot_ax.get_xticklabels(), rotation=15, ha='right')
            subplot_ax.grid(True, axis='y', linewidth=0.5, color=grid_color, zorder=0)
            
            for spine in subplot_ax.spines.values():
                spine.set_linewidth(0.6)
                spine.set_color(border_color)
                
    plt.subplots_adjust(left=0.08, right=0.88, top=0.92, bottom=0.12, hspace=0.4, wspace=0.3)
    
    if output_name:
        plt.savefig(output_name, bbox_inches='tight', dpi=300)
    # plt.show()


def plot_pe_results(df, config, output_name=None):
    # Merge experimental results with the config dataframe
    merged_df = df.merge(config, left_on='config_name', right_on='config_Id', how='left')
    
    # Process categorical labels for plotting
    merged_df['Clustering Method'] = merged_df.apply(clean_algo_names, axis=1)
    
    # Standardize names to match your manual list exactly
    merged_df['Distance Measure'] = merged_df['distance_measure'].str.replace('_distance', '').str.title()
    name_map = {
        'Euclidian': 'Euclidean',
        'Euclidean': 'Euclidean',
        'Chebysehv': 'Chebyshev',
        'Chebyshev': 'Chebyshev',
        'Weighted_Jaccard': 'Weight Jaccard',
        'Jensenshannon': 'Jensen-Shannon',
        'Sorensen': 'Sørensen'
    }
    merged_df['Distance Measure'] = merged_df['Distance Measure'].replace(name_map)

    # Calculate average PPI score from StringDB metrics
    ppi_cols = [
        'shared_PPI_partners_score_stringdb', 
        'shortest_PPI_path_score_stringdb',
        'density_score_stringdb', 
        'lcc_score_stringdb', 
        'tc_score_stringdb'
    ]
    merged_df['Avg PPI Score'] = merged_df[ppi_cols].mean(axis=1)

    # Calculate average Functional Similarity score from GO metrics
    go_cols = ['goa_similarity_lin', 'goa_similarity_resnik', 'goa_similarity_jc']
    merged_df['Avg Func Sim Score'] = merged_df[go_cols].mean(axis=1)

    # Pivot datasets into Mean Heatmap matrices
    pivot_sil = merged_df.groupby(['Distance Measure', 'Clustering Method'])['silhouette_score'].mean().unstack()
    pivot_ppi = merged_df.groupby(['Distance Measure', 'Clustering Method'])['Avg PPI Score'].mean().unstack()
    pivot_fun = merged_df.groupby(['Distance Measure', 'Clustering Method'])['Avg Func Sim Score'].mean().unstack()

    # Define exact manual coordinate ordering
    dist_order = [d for d in ['Manhattan', 'Canberra', 'Euclidean', 'Chebyshev', 'Cosine', 'Sørensen', 'Weight Jaccard', 'Jensen-Shannon'] if d in pivot_sil.index]
    algo_order = [a for a in ['CD (P-50)', 'CD (P-75)', 'CD (KNN)', 'Hierarchical (Inc)', 'Hierarchical (OC)', 'HDBSCAN'] if a in pivot_sil.columns]

    # Reindex matrices to guarantee matching visual ordering
    pivot_sil = pivot_sil.reindex(index=dist_order, columns=algo_order)
    pivot_ppi = pivot_ppi.reindex(index=dist_order, columns=algo_order)
    pivot_fun = pivot_fun.reindex(index=dist_order, columns=algo_order)

    # Initialize 3x2 grid sized for a standard half A4 page (approx 8.3 x 5.8 inches)
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(8.3, 5.8))
    plt.rcParams['font.family'] = 'sans-serif'
    border_color = ut.color_palette["medium_gray"]
    
    # Panel 1A (Row 0, Col 0)
    ax_1a = axes[0, 0]
    sns.heatmap(pivot_sil, ax=ax_1a, cmap='viridis', annot=False, cbar_kws={'label': 'Mean Silhouette'})
    ax_1a.set_title('1A. Cluster Quality (Silhouette Score)', fontsize=9, loc='left', color=ut.color_palette["midnight_blue"], pad=8)

    # Panel 1B (Row 0, Col 1)
    ax_1b = axes[0, 1]
    sns.heatmap(pivot_ppi, ax=ax_1b, cmap='magma', annot=False, cbar_kws={'label': 'Mean PPI'})
    ax_1b.set_title('1B. Avg PPI Score', fontsize=9, loc='left', color=ut.color_palette["midnight_blue"], pad=8)

    # Panel 2A (Row 1, Col 0)
    ax_2a = axes[1, 0]
    sns.heatmap(pivot_fun, ax=ax_2a, cmap='inferno', annot=False, cbar_kws={'label': 'Mean GO Sim'})
    ax_2a.set_title('2A. Avg GO Similarity Score', fontsize=9, loc='left', color=ut.color_palette["midnight_blue"], pad=8)

    # Hide unused panels in the 3x2 matrix
    axes[1, 1].axis('off')
    axes[2, 0].axis('off')
    axes[2, 1].axis('off')

    # Layout cleanup for active heatmaps
    heatmap_axes = [ax_1a, ax_1b, ax_2a]
    for subplot_ax in heatmap_axes:
        subplot_ax.set_xlabel('Clustering Method', fontsize=8, color=ut.color_palette["dark_navy"])
        subplot_ax.set_ylabel('Distance Measure', fontsize=8, color=ut.color_palette["dark_navy"])
        subplot_ax.tick_params(axis='both', labelsize=7, colors=ut.color_palette["dark_navy"])
        subplot_ax.set_xticklabels(subplot_ax.get_xticklabels(), rotation=25, ha='right')
        subplot_ax.set_yticklabels(subplot_ax.get_yticklabels(), rotation=0)
        
        for spine in subplot_ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_color(border_color)
                
    plt.tight_layout()
    
    if output_name:
        plt.savefig(output_name, bbox_inches='tight', dpi=300)


def plot_ps_ppi_grid(df, config, output_name=None):
# Merge experimental results with the config dataframe
    merged_df = df.merge(config, left_on='config_name', right_on='config_Id', how='left')
    
    # Process categorical labels for plotting using clean helper
    merged_df['Clustering Method'] = merged_df.apply(clean_algo_names, axis=1)
    
    # Standardize names to match your manual list exactly
    merged_df['Distance Measure'] = merged_df['distance_measure'].str.replace('_distance', '').str.title()
    name_map = {
        'Euclidian': 'Euclidean',
        'Euclidean': 'Euclidean',
        'Chebysehv': 'Chebyshev',
        'Chebyshev': 'Chebyshev',
        'Weighted_Jaccard': 'Weight Jaccard',
        'Jensenshannon': 'Jensen-Shannon',
        'Sorensen': 'Sørensen'
    }
    merged_df['Distance Measure'] = merged_df['Distance Measure'].replace(name_map)

    # List of the five individual target PPI columns
    ppi_cols = [
        'shared_PPI_partners_score_stringdb', 
        'shortest_PPI_path_score_stringdb',
        'density_score_stringdb', 
        'lcc_score_stringdb', 
        'tc_score_stringdb'
    ]
    
    # Transform data from wide to long format for multi-boxplot mapping
    melted_ppi = merged_df.melt(
        id_vars=['Distance Measure', 'Clustering Method'],
        value_vars=ppi_cols,
        var_name='PPI_Metric',
        value_name='Score'
    )
    
    # Rename variables to clear short abbreviations for micro-axes
    metric_map = {
        'shared_PPI_partners_score_stringdb': 'PPI1',
        'shortest_PPI_path_score_stringdb': 'PPI2',
        'density_score_stringdb': 'PPI3',
        'lcc_score_stringdb': 'PPI4',
        'tc_score_stringdb': 'PPI5'
    }
    melted_ppi['PPI_Metric'] = melted_ppi['PPI_Metric'].replace(metric_map)

    # Establish exact row and column visual sort indexing
    dist_order = [d for d in ['Manhattan', 'Canberra', 'Euclidean', 'Chebyshev', 'Cosine', 'Sørensen', 'Weight Jaccard', 'Jensen-Shannon'] if d in melted_ppi['Distance Measure'].unique()]
    algo_order = [a for a in ['CD (P-50)', 'CD (P-75)', 'CD (KNN)', 'Hierarchical (Inc)', 'Hierarchical (OC)', 'HDBSCAN'] if a in melted_ppi['Clustering Method'].unique()]

    # Initialize FacetGrid scaled to fit standard A4 paper dimensions
    plt.rcParams['font.family'] = 'sans-serif'
    g = sns.FacetGrid(
        melted_ppi, 
        row="Distance Measure", 
        col="Clustering Method",
        row_order=dist_order,
        col_order=algo_order,
        margin_titles=True,
        height=1.2,
        aspect=1.4
    )
    
    # Map the five colorful boxplots across every single coordinate element
    metric_colors = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3"]
    g.map_dataframe(
        sns.boxplot, 
        x="PPI_Metric", 
        y="Score", 
        hue="PPI_Metric",
        palette=metric_colors,
        width=0.7,
        showfliers=False,
        linewidth=0.7
    )
    
    # Axis configuration and titles
    g.set_titles(row_template="{row_name}", col_template="{col_name}")
    g.set_axis_labels("", "")
    
    # Grid cleanup loop across cells
    tick_positions = range(5)
    for (row_val, col_val), subplot_ax in g.axes_dict.items():
        subplot_ax.tick_params(axis='both', which='major', labelsize=7, colors=ut.color_palette["dark_navy"])
        subplot_ax.grid(True, axis='y', linestyle=':', alpha=0.5, linewidth=0.5, color=ut.color_palette["light_gray"])
        
        # Safe tick setting strategy preventing deprecation warnings
        subplot_ax.set_xticks(tick_positions)
        subplot_ax.set_xticklabels(list(metric_map.values()), rotation=45, ha='right', fontsize=6)
        
        for spine in subplot_ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_color(ut.color_palette["medium_gray"])
            
    # Shrink margin text font sizes to completely prevent labels from overlapping
    for ax in g.figure.axes:
        for text in ax.texts:
            text.set_fontsize(7)

    # Adjust layout boundaries leaving a solid margin block at the base for the legend text
    plt.subplots_adjust(wspace=0.2, hspace=0.45, bottom=0.14)
    
    # Build clean line proxies to guarantee zero extraction errors from axes
    legend_patches = [plt.Line2D([0], [0], color=c, lw=4) for c in metric_colors]
    full_labels =list(metric_map.values())
    
    # Standard figure level call avoiding structural deprecation paths
    g.figure.legend(
        handles=legend_patches,
        labels=full_labels,
        loc='lower center',
        bbox_to_anchor=(0.5, 0.02),
        ncol=5,
        fontsize=9,
        title="Protein-Protein Interaction (PPI) Metrics",
        title_fontsize=10,
        frameon=True
    )
    
    if output_name:
        plt.savefig(output_name, bbox_inches='tight', dpi=300)