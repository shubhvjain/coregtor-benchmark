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

    data_path = output_path / "target_stats.json"
    if not data_path.exists():
        raise ValueError("target_stats.json does not exists")
    with open(data_path) as f:
        target_file = json.load(f)
    # print(target_file)
    data = pd.json_normalize(data=target_file)
    data["dataset_name"] = dataset_info["name"]
    data["dataset_abbr"] = dataset_info["abbr"]
    data["dataset_sample_size"] = dataset_info["n_sample"]
    data["exp_id"] = exp["id"]
    return data

def _generate_save_stats_diagram(df, folder_path):
    # 'paper' context ensures labels and fonts are appropriately scaled for publication

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
        ax.tick_params(axis="y", labelsize=9)

        # Only show the "Dataset" x-axis title on the bottom row to reduce clutter
        if i >= 4:
            ax.set_xlabel("Dataset", fontsize=9)
        else:
            ax.set_xlabel("")

    # Save and prevent memory leaks
    fig.savefig(folder_path / "gene_frequency_histogram.pdf",
                bbox_inches="tight")
    fig.savefig(folder_path / "gene_frequency_histogram.svg",
                bbox_inches="tight")
    plt.close(fig)


def _generate_save_runtime_diagram(df, folder_path):
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
    sns.boxplot(data=df, x="dataset_abbr", y="stats.total_time", ax=ax_list[1],
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

    out_file = folder / f"gene_freq_hist.svg"

    if out_file.exists() and not rerun:
        print("File already exists")
        return

    data = []
    for e in input["experiments"]:
        o, t = ut.get_exp_path({"id": e}, env)
        with open(o/"input.json") as f:
            input_data = json.load(f)
        details = _get_exp_stats(input_data, o, t)
        data.append(details)

    df = pd.concat(data, ignore_index=True)
    df.to_csv(folder/"target_stats.csv", index=False)
    sort_col = options.get("sort_by", "dataset_name")
    df.sort_values(by=sort_col, inplace=True)
    _generate_save_stats_diagram(df, folder)
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

    plot_landscape_metrics(df, metrics_list, out_file, False)
