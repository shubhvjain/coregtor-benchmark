import argparse
import json
import sys
from pathlib import Path
from dotenv import dotenv_values
import src.results.util as ut
from joblib import Parallel, delayed, load
import pandas as pd
import matplotlib.pyplot as plt

from src.results.explore import  plot_skewness_dist, plot_metric_boxplot

from src.analysis.gtex import generate_gtex_stats,gtex_stats_diagram,generate_gtex_dcor_cache

from src.analysis.target import combined_target_diagram, parameter_config_plot


def hello(input, env, options, args):
    target = input.get("target", "world")
    print(f"Hello, {target}!")
    print(f"  env keys : {list(env.keys())}")
    print(f"  verbose  : {args.verbose}")


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



METHOD_REGISTRY = {
    "hello": hello,
    "generate_gtex_stats":generate_gtex_stats,
    "combined_target_diagram":combined_target_diagram,
    "combined_dataset_target_diagram": combined_dataset_target_diagram,
    "parameter_config_plot":parameter_config_plot,
    "gtex_stats_diagram":gtex_stats_diagram,
    "generate_gtex_dcor_cache":generate_gtex_dcor_cache
}

# CLI entry point


def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate analysis",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--env",   required=True,
                        metavar="PATH", help="Path to the .env file.")
    parser.add_argument("--input", required=True,
                        metavar="PATH", help="Path to the JSON input file.")
    parser.add_argument(
        "--id",
        required=True,
        metavar="ID",
        help="ID of the result object in input['results'] to run.",
    )

    # Optional arguments — add more as needed
    parser.add_argument("--njobs", type=int, default=-1,
                        help="Number of parallel jobs. Defaults to -1 (all cores).")
    parser.add_argument("--rerun", action="store_true", default=False,
                        help="Force rerun even if output already exists.")
    parser.add_argument("--verbose", action="store_true",
                        default=False, help="Enable verbose output.")
    parser.add_argument("--dry-run", action="store_true",
                        default=False, help="Load files but skip execution.")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    env = ut.get_env(args.env)
    input_data = ut.get_input(args.input)

    options = ut.get_result_by_id(input_data, args.id)

    result_type = options.get("type")
    if result_type not in METHOD_REGISTRY:
        print(
            f"[error] unknown type '{result_type}'. Registered types: {', '.join(METHOD_REGISTRY.keys())}", file=sys.stderr)
        sys.exit(1)

    fn = METHOD_REGISTRY[result_type]

    if args.dry_run:
        print(
            f"[dry-run] would call '{result_type}' for id '{args.id}' with {len(env)} env key(s).")
        return

    fn(input_data, env, options, args)


if __name__ == "__main__":
    main()
