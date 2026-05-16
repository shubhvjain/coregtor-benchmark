"""
Stats for GTEx tissue data
"""
import src.results.util as ut
from pathlib import Path
import os
from coregtor.workflow.util import read_dataset, get_tflist, get_protein_coding_genes
from joblib import dump, load
import pandas as pd
import numpy as np

def generate_gtex_stats(input, env, options, args):
    """generate basic stats about GTEx tissue data"""
    apath = ut.get_analysis_path(env)
    folder = apath/"gtex_info"
    folder.mkdir(exist_ok=True, parents=True)
    rerun = args.rerun

    main_stats_file = folder/"stats.pkl"

    if main_stats_file.exists() and not rerun:
        print("main stats file already exists")
        return

    all_datasets = sorted(list(input["datasets"].keys()))
    print(all_datasets)

    CONFIG = {
        "data_path": Path(env["DATA_PATH"]),
        "temp_path": Path(env["EXP_TEMP_PATH"]),
        "out_path": Path(env["EXP_OUTPUT_PATH"])
    }
    tf = get_tflist(CONFIG)

    detailed_stats = {}
    summary_stats = []

    for d in all_datasets:
        dets = input["datasets"][d]

        dets["path"] = os.path.expandvars(f"$DATA_PATH/{dets['path']}")
        # print(os.path.expandvars(dets["path"]))
        # print(dets)
        data = read_dataset(dets, CONFIG)
        # print(data)
        stats = compute_all_genes_stats(data, tf,CONFIG=CONFIG)
        detailed_stats[d] = stats
        summary = {
            "name": dets.get("tissue_name"),
            "abbr": d,
            "title": dets.get("title"),
            "n_sample": stats.get("n_samples"),
            "n_genes": stats.get("n_genes"),
            "n_tf": stats.get("n_tf"),
            "n_pc": stats.get("n_pc"),
            'sparsity': stats.get("sparsity"),
            'sparsity_tf': stats.get("sparsity_tf"),
            'sparsity_pc': stats.get("sparsity_pc"),
        }
        summary_stats.append(summary)

    dump(detailed_stats, main_stats_file,compress=3)
    data = pd.DataFrame(summary_stats)
    data.to_csv(folder/"stats.csv", index=False)

    # generate the latex data table
    cols_to_keep = [
        'name', "abbr", 'n_sample', 'n_genes',
        'sparsity', 'n_tf', 'sparsity_tf', 'n_pc', 'sparsity_pc'
    ]
    data_latex = data[cols_to_keep].copy()

    # Rename columns
    rename_dict = {
        'name': 'Tissue Name',
        "abbr": "Key",
        'n_sample': 'Total Samples',
        'n_genes': 'Total Genes',
        'sparsity': 'Overall Sparsity',
        'n_tf': "TF Count",
        'sparsity_tf': 'TF Sparsity',
        'n_pc': "PCG Count",
        'sparsity_pc': 'PCG Sparsity'
    }
    data_latex = data_latex.rename(columns=rename_dict)

    data_latex['Overall Sparsity'] = data_latex['Overall Sparsity'].apply(
        lambda x: f"{x:.3f}")
    data_latex['TF Sparsity'] = data_latex['TF Sparsity'].apply(
        lambda x: f"{x:.3f}")
    data_latex['PCG Sparsity'] = data_latex['PCG Sparsity'].apply(
        lambda x: f"{x:.3f}")

    data_latex.to_latex(folder / "stats.tex", index=False)


def compute_all_genes_stats(ge, tf_list, sparsity_thresholds=None,CONFIG=None):
    """
    Computes global and per-gene statistics for the entire dataset.
    If sparsity_thresholds is None, it defaults to a range [0.05, 0.10 ... 1.00].
    """

    # Generate default threshold array if None is provided
    if sparsity_thresholds is None:
        sparsity_thresholds = [round(t, 2)
                               for t in np.arange(0.05, 1.05, 0.05)]

    # DROP DUPLICATE COLUMNS FIRST: Keep the first occurrence of each gene
    ge = ge.loc[:, ~ge.columns.duplicated(keep='first')]

    n_samples, n_genes = ge.shape
    all_genes = ge.columns.tolist()

    # 1. Gene filters
    protein_list = get_protein_coding_genes(all_genes,CONFIG)
    protein_set = set(protein_list)

    tf_found = [tf for tf in tf_list if tf in all_genes]
    tf_set = set(tf_found)

    # 2. Vectorized Statistical Computations
    gene_mean = ge.mean(axis=0)
    gene_var = ge.var(axis=0)
    gene_sum = ge.sum(axis=0)

    log_ge = np.log2(ge + 1)
    log_mean = log_ge.mean(axis=0)
    log_var = log_ge.var(axis=0)

    zero_counts = (ge == 0).sum(axis=0)
    sparsity_fraction = zero_counts / n_samples

    # 3. Compute Sparsity Threshold Distributions (Plot-Friendly Format)
    tf_sparsity = sparsity_fraction[tf_found]
    pc_sparsity = sparsity_fraction[protein_list]

    total_tf = len(tf_found)
    total_pc = len(protein_list)

    # Use lists instead of nested dicts for easy plotting
    tf_threshold_dist = {'threshold': [], 'count': [], 'percentage': []}
    pc_threshold_dist = {'threshold': [], 'count': [], 'percentage': []}

    for t in sparsity_thresholds:
        t_rounded = round(t, 2)

        # TF counts
        tf_count = int((tf_sparsity <= t).sum())
        tf_pct = round((tf_count / total_tf * 100),
                       2) if total_tf > 0 else 0.0

        tf_threshold_dist['threshold'].append(t_rounded)
        tf_threshold_dist['count'].append(tf_count)
        tf_threshold_dist['percentage'].append(tf_pct)

        # PC counts
        pc_count = int((pc_sparsity <= t).sum())
        pc_pct = round((pc_count / total_pc * 100),
                       2) if total_pc > 0 else 0.0

        pc_threshold_dist['threshold'].append(t_rounded)
        pc_threshold_dist['count'].append(pc_count)
        pc_threshold_dist['percentage'].append(pc_pct)

    # 4. Build the gene_details Dictionary
    gene_details = {}
    for gene in all_genes:
        gene_details[gene] = {
            'mean': float(gene_mean[gene]),
            'variance': float(gene_var[gene]),
            'log_mean': float(log_mean[gene]),
            'log_variance': float(log_var[gene]),
            'zero_count': int(zero_counts[gene]),
            'sparsity': float(sparsity_fraction[gene]),
            'is_protein_coding': gene in protein_set,
            'is_tf': gene in tf_set
        }

    # Global sparsity across all genes
    total_elements = n_samples * n_genes
    global_sparsity = float(
        zero_counts.sum() / total_elements) if total_elements > 0 else 0.0

    # Global sparsity across TF genes
    tf_elements = n_samples * total_tf
    sparsity_tf = float(zero_counts[tf_found].sum(
    ) / tf_elements) if tf_elements > 0 else 0.0

    # Global sparsity across protein coding genes
    pc_elements = n_samples * total_pc
    sparsity_pc = float(zero_counts[protein_list].sum(
    ) / pc_elements) if pc_elements > 0 else 0.0

    # Construct Final Dictionary
    results = {
        'n_samples': n_samples,
        'n_genes': n_genes,
        'n_tf': total_tf,
        'n_pc': total_pc,
        'sparsity': global_sparsity,
        'sparsity_tf': sparsity_tf,
        'sparsity_pc': sparsity_pc,
        'tf_sparsity_distribution': tf_threshold_dist,
        'pc_sparsity_distribution': pc_threshold_dist,
        'protein_coding_genes': protein_list,
        'tf_genes': tf_found,
        'gene_details': gene_details
    }

    return results
