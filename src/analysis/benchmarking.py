import src.results.util as ut
from pathlib import Path
import os
from src.pipeline.util import read_dataset, get_tflist, get_protein_coding_genes, get_coreglist
from joblib import dump, load
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import pyarrow as pa
import pyarrow.parquet as pq
import src.results.dcorr as dc

import math
import matplotlib.pyplot as plt
import seaborn as sns
import src.results.util as ut


def generate_combined_indices_file(input, env, options, args):
    """"""
    apath = ut.get_analysis_path(env)
    folder = apath/input["id"]
    folder.mkdir(exist_ok=True, parents=True)
    rerun = args.rerun

    main_file = folder/"data.paraquet"

    if main_file.exists() and not rerun:
        print("main stats file already exists")
        return

    results = []
    result_paths = input.get("results_paths")
  
    exp_folder = Path(env.get("EXP_OUTPUT_PATH"))

    coregtor_name = input["coregtor_options"]["result_name"]
    include_optimal_only = True
    for ct in result_paths["coregtor"]:
        res_file = exp_folder / ct / f"{coregtor_name}_indices.csv"
        res = pd.read_csv(res_file)
        if include_optimal_only:
            score_cols = ["goa_similarity_jc","goa_similarity_lin", "goa_similarity_resnik", "TFBS_affinity_score", "TFBS_affinity_sum_score", "dCor_sources_score","dCor_target_score","shared_PPI_partners_score_stringdb","shortest_PPI_path_score_stringdb","density_score_stringdb", "lcc_score_stringdb","tc_score_stringdb"]

            res["combined_evidence_score"] = res[score_cols].sum(axis=1)
            res = res.loc[res.groupby("target")["combined_evidence_score"].idxmax()].reset_index(drop=True)

        res["tool"] = "coregtor"
        res["experiment"] = ct
        results.append(res)

    coregnet_name = input["coregnet_options"]["result_name"]
    for cn in result_paths["coregnet"]:
        res_file = exp_folder / cn / f"{coregnet_name}_indices.csv"
        res = pd.read_csv(res_file)
        res["tool"] = "coregnet"
        res["experiment"] = cn
        results.append(res)
    
    combined = pd.concat(results)
    #  save all results , unedited
    raw_file = folder / "data_raw.paraquet"
    combined.to_parquet(raw_file)

    combined["n_source"] = combined["n_source"].fillna(combined["n_sources"])
    combined["cluster_uid"]  = combined["cluster_uid"].fillna(combined["uid"])
    combined = combined.drop(columns=["n_sources","uid","Unnamed: 0"])
    combined = combined[combined["n_source"]<100]
    combined.to_parquet(main_file)





def plot_tool_index_boxplots(
    df,
    index_cols,
    tool_col="tool",
    n_cols=4,
    save_path=None
):
    sns.set_theme(style="whitegrid", rc=ut.default_sns_configs)

    tool_order = list(df[tool_col].dropna().unique())

    palette = {
        tool: ut.tool_map[tool]["color"]
        for tool in tool_order
        if tool in ut.tool_map
    }

    tick_labels = [
        ut.tool_map.get(tool, {}).get("label", tool)
        for tool in tool_order
    ]

    n_rows = math.ceil(len(index_cols) / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * 3.0, n_rows * 2.6),
        squeeze=False
    )
    axes = axes.flatten()

    for i, col in enumerate(index_cols):
        ax = axes[i]

        sns.boxplot(
            data=df,
            x=tool_col,
            y=col,
            order=tool_order,
            palette=palette,
            linewidth=0.8,
            showfliers=False,
            ax=ax
        )

        ax.set_title(ut.col_index_map.get(col, col))
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(tick_labels, rotation=0)

    for j in range(len(index_cols), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")

    # plt.show()



def tool_benchmark_combined_plot(input, env, options, args):
    """"""
    apath = ut.get_analysis_path(env)
    folder = apath/input["id"]
    folder.mkdir(exist_ok=True, parents=True)
    rerun = args.rerun

    figure_file = folder/"plot.svg"

    if figure_file.exists() and not rerun:
        print("figure already exists")
        return

    data_file = folder / "data.paraquet"
    data = pd.read_parquet(data_file)
    
    data = data[ data['n_source'] < 75 ]

    fig_list = [
        {
          "name":"plot_stringdb",
          "col":  ut.index_col2
        },
        {
          "name":"plot_biogrid",
          "col":  ut.index_col1
        },
        {
          "name":"plot_hippie",
          "col":  ut.index_col3
        },
    ]
    for fig in fig_list:
      fpath = folder/f"{fig['name']}.svg"
      if fpath.exists() and not rerun:
          continue
      plot_tool_index_boxplots(data, fig["col"],save_path=fpath)




