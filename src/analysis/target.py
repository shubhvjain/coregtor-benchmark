"""
"""
import src.results.util as ut
from src.analysis.parameter_selection import plot_all_configs_raw,generate_ps_plot1
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


def _get_exp_index_file(result_name, output_path):
    """"""

    parquet_path = output_path / f"cluster_indices_{result_name}.parquet"
    csv_path = output_path / f"cluster_indices_{result_name}.csv.gz"

    # Check for Parquet first, then fall back to JSON
    if parquet_path.exists():
        data = pd.read_parquet(parquet_path)
    elif csv_path.exists():
        data = pd.pd.read_csv(csv_path)
    else:
        raise ValueError(
            f"Neither {parquet_path.name} nor {csv_path.name} exists"
        )
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

    plt.rcParams.update(ut.default_sns_configs)

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

def combined_gene_freq_diagram(input, env, options, args):
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


def parameter_config_detailed_plot(input, env, options, args):
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
    freq_data = []
    
    for e in input["experiment_results_paths"]:
        o, t = ut.get_exp_path({"path": e}, env)

        with open(o/"input.json") as f:
            input_data = json.load(f)
        dataset_info = ut.get_exp_gtex_info(input_data["dataset"]["name"])
        
        details = _get_exp_index_file(options['result_name'],o) 
        details["exp_id"] = e
        details["dataset_name"] = dataset_info["name"]
        details["dataset_abbr"] = dataset_info["abbr"]
        data.append(details)

        # details1 = pd.read_csv(o/f"freq_coreg_indices_{options['result_name']}.csv.gz")
        # details1["exp_id"] = e
        # details1["dataset_name"] = dataset_info["name"]
        # details1["dataset_abbr"] = dataset_info["abbr"]
        # freq_data.append(details1)


    df = pd.concat(data, ignore_index=True)
    sf = options.get("save_file", False)
    if sf:
        df.to_csv(folder/f"clusters_{options['result_name']}.csv.gz", index=False)
    
    freq_df = None # pd.concat(freq_data, ignore_index=True)
    # print(freq_df)

    default_metrics = ['silhouette_score','n_source','cluster_density', 'cluster_diameter', 
       'shortest_PPI_path_score_hippie', 'shortest_PPI_path_score_stringdb',
       'shortest_PPI_path_score_biogrid', 'shared_PPI_partners_score_hippie',
       'shared_PPI_partners_score_stringdb',
       'shared_PPI_partners_score_biogrid', 'goa_similarity_lin',
       'goa_similarity_resnik', 'goa_similarity_jc', 'density_hippie',
       'density_score_hippie', 'lcc_hippie', 'lcc_score_hippie', 'tc_hippie',
       'tc_score_hippie', 'node_found_ratio_hippie', 'density_stringdb',
       'density_score_stringdb', 'lcc_stringdb', 'lcc_score_stringdb',
       'tc_stringdb', 'tc_score_stringdb', 'node_found_ratio_stringdb',
       'density_biogrid', 'density_score_biogrid', 'lcc_biogrid',
       'lcc_score_biogrid', 'tc_biogrid', 'tc_score_biogrid',
       'node_found_ratio_biogrid', 'dCor_sources', 'dCor_target', 'DC1', 'DC2',
       'TFBS_affinity', 'TF_found_per', 'TFBS_affinity_score']
    
    freq_metrics = ['dCor_targets_score','ppi_shared_partners_hippie','ppi_shared_partners_biogrid','go_score']

    metrics_list = options.get("metrics_list", default_metrics)
    
    # raw config comparison 
    plot_all_configs_raw(df,metrics_list,freq_df, [] ,out_file)

    cfilters = options.get("custom_plots",[])

    for c in cfilters:
        cd_configs = c.get("configs",[])
        out_file_cd = folder / f"param_selection_cd_{c['name']}.svg"
        dfcd  = df[df["config_name"].isin(cd_configs)]
        plot_all_configs_raw(dfcd,metrics_list,freq_df, [] ,out_file_cd)


    # result_plot1 = folder / f"param_selection_{options['result_name']}.pdf"
    #generate_ps_plot1(df,freq_df,result_plot1)
    #plot_pe_results(df,pconfig,result_plot1)
    #result_plot2 = folder / f"param_selection_ppi_{options['result_name']}.pdf"
    #plot_ps_ppi_grid(df,pconfig,result_plot2)


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
    freq_data = []
    
    for e in input["experiment_results_paths"]:
        o, t = ut.get_exp_path({"path": e}, env)

        with open(o/"input.json") as f:
            input_data = json.load(f)
        dataset_info = ut.get_exp_gtex_info(input_data["dataset"]["name"])
        
        details = _get_exp_index_file(options['result_name'],o)  
        details["exp_id"] = e
        details["dataset_name"] = dataset_info["name"]
        details["dataset_abbr"] = dataset_info["abbr"]
        data.append(details)

        # details1 = pd.read_csv(o/f"freq_coreg_indices_{options['result_name']}.csv.gz")
        # details1["exp_id"] = e
        # details1["dataset_name"] = dataset_info["name"]
        # details1["dataset_abbr"] = dataset_info["abbr"]
        # freq_data.append(details1)


    df = pd.concat(data, ignore_index=True)
    sf = options.get("save_file", False)
    if sf:
        df.to_csv(folder/f"clusters_{options['result_name']}.csv.gz", index=False)
    
    freq_df = None # pd.concat(freq_data, ignore_index=True)
    # print(freq_df)
    c_name = input.get("parameter_config_name","ps_configs")
    print(c_name)
    pconfig =  ut.get_ps_configs(c_name)

    result_plot1 = folder / f"param_{options['result_name']}.pdf"
    generate_ps_plot1(df,freq_df,pconfig,options,result_plot1)
    

    # result_plot1 = folder / f"param_selection_{options['result_name']}.pdf"
    #generate_ps_plot1(df,freq_df,result_plot1)
    #plot_pe_results(df,pconfig,result_plot1)
    #result_plot2 = folder / f"param_selection_ppi_{options['result_name']}.pdf"
    #plot_ps_ppi_grid(df,pconfig,result_plot2)
