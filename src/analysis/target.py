"""
"""
import src.results.util as ut


def combined_target_diagram(input, env, options, args):
  """
  """
  apath = ut.get_analysis_path(env)
  folder = apath/input["id"]
  folder.mkdir(exist_ok=True, parents=True)
  rerun = args.rerun

  out_file = folder / f"gene_freq_hist_stats.svg"

  if out_file.exists() and not rerun:
    print("File already exists")
    return


def combined_dataset_target_diagram(input, env, options, args):
    """"""
    apath = ut.get_analysis_path(env)
    folder = apath/"parameter_selection"
    folder.mkdir(exist_ok=True, parents=True)
    rerun = args.rerun

    out_file = folder / f"{input['id']}_dataset_stats.svg"

    if out_file.exists() and not rerun:
        print("File already exists")
        return

    
    dfs = []
    for i in input["exp_list"]:
        print(i)
        data_path = Path(env["EXP_OUTPUT_PATH"])/i/"target_stats.json"
        if not data_path.exists():
            print("target_stats.json does not exists")
            continue

        with open(data_path) as f:
            datajs = json.load(f)
    
        data = pd.json_normalize(data=datajs)
        #print(data)
        data["dataset"]=i
        dfs.append(data)


    df = pd.concat(dfs,ignore_index=True)
    print(df)
    df.to_csv(folder/f"{input['id']}_dataset_stats.csv")

    # Setup 4 subplots with width ratios [2, 2, 2, 4] to match [20, 20, 20, 40]
    fig, axes = plt.subplots(
        1, 4, 
        figsize=(10,3), 
        gridspec_kw={'width_ratios': [2, 2, 2, 4]}
    )

    # 1. Total Root Nodes (Boxplot)
    plot_metric_boxplot(data, y_col='stats.n_unique_roots',title="Root node count", ax=axes[0],ylabel="Count")
    
    # 2. Active Source Percentage (Boxplot)
    plot_metric_boxplot(data, y_col='stats.total_time',title="Run time", ax=axes[1],ylabel="Time(Seconds)")

    data["gf_coverage_per"] = 100 - data["gf.sparsity_pcr"]
    
    # 3. Sparsity Percentage (Boxplot)
    plot_metric_boxplot(data, y_col='gf_coverage_per',title="Gene Frequency Coverage", ax=axes[2],ylabel="Percent")
    
    # 4. Skewness Distribution (KDE)
    plot_skewness_dist(data, col_name="gf.count_skew" ,ax=axes[3])

    plt.tight_layout()
    # Save as SVG 
    plt.savefig(out_file)  