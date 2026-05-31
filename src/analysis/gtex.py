"""
Stats for GTEx tissue data
"""
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
    tf = get_tflist(env)
    coreg = get_coreglist(env)

    detailed_stats = {}
    summary_stats = []

    for d in all_datasets:
        dets = input["datasets"][d]

        dets["path"] = os.path.expandvars(f"$DATA_PATH/{dets['path']}")
        # print(os.path.expandvars(dets["path"]))
        # print(dets)
        data = read_dataset(dets, env)
        # print(data)
        stats = compute_all_genes_stats(data, tf,CONFIG=env,coreg=coreg)
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


def compute_all_genes_stats(ge, tf_list, sparsity_thresholds=None,CONFIG=None,coreg=None):
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

    coreg_found = [tf for tf in coreg if tf in all_genes]
    coreg_set = set(coreg_found)

    regulators = tf_set | coreg_set
    reg_found = [tf for tf in list (regulators) if tf in all_genes]
    reg_set = set(reg_found)

    

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
    #print(pc_sparsity)
    reg_sparsity = sparsity_fraction[reg_found]

    total_tf = len(tf_found)
    total_pc = len(protein_list)
    total_coreg = len(coreg_found)
    total_reg =  len(reg_found)

    # print("---")
    # print(total_tf)
    # print(total_pc)
    # print(total_coreg)
    # print(total_reg)

    # Use lists instead of nested dicts for easy plotting
    tf_threshold_dist = {'threshold': [], 'count': [], 'percentage': []}
    pc_threshold_dist = {'threshold': [], 'count': [], 'percentage': []}
    reg_threshold_dist = {'threshold': [], 'count': [], 'percentage': []}


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

        # reg counts
        reg_count = int((reg_sparsity <= t).sum())
        reg_pct = round((reg_count / total_reg * 100),
                       2) if total_reg > 0 else 0.0

        reg_threshold_dist['threshold'].append(t_rounded)
        reg_threshold_dist['count'].append(reg_count)
        reg_threshold_dist['percentage'].append(reg_pct)

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
            'is_tf': gene in tf_set,
            'is_coreg': gene in coreg_set,
            'is_reg': gene in reg_set
        }

    # Global sparsity across all genes
    total_elements = n_samples * n_genes
    global_sparsity = float(
        zero_counts.sum() / total_elements) if total_elements > 0 else 0.0

    # Global sparsity across TF genes
    tf_elements = n_samples * total_tf
    sparsity_tf = float(zero_counts[tf_found].sum(
    ) / tf_elements) if tf_elements > 0 else 0.0

     # Global sparsity across Coreg genes
    cr_elements = n_samples * total_coreg
    sparsity_cr = float(zero_counts[coreg_found].sum(
    ) / cr_elements) if cr_elements > 0 else 0.0

    # Global sparsity across protein coding genes
    pc_elements = n_samples * total_pc
    sparsity_pc = float(zero_counts[protein_list].sum(
    ) / pc_elements) if pc_elements > 0 else 0.0

     # Global sparsity across reg genes
    re_elements = n_samples * total_reg
    sparsity_re = float(zero_counts[reg_found].sum(
    ) / re_elements) if re_elements > 0 else 0.0

    # Construct Final Dictionary
    results = {
        'n_samples': n_samples,
        'n_genes': n_genes,
        'n_tf': total_tf,
        'n_pc': total_pc,
        'n_coreg': total_coreg,
        'n_reg': total_reg,
        'sparsity': global_sparsity,
        'sparsity_tf': sparsity_tf,
        'sparsity_pc': sparsity_pc,
        'sparsity_coreg': sparsity_cr,
        'sparsity_reg': sparsity_re,
        'tf_sparsity_distribution': tf_threshold_dist,
        'pc_sparsity_distribution': pc_threshold_dist,
        'reg_sparsity_distribution': reg_threshold_dist,
        'protein_coding_genes': protein_list,
        'tf_genes': tf_found,
        'coreg_genes': coreg_found,
        'reg_genes':reg_found,
        'gene_details': gene_details
    }

    return results




def gtex_stats_diagram(input, env, options, args):
    """generate basic stats about GTEx tissue data"""
    apath = ut.get_analysis_path(env)
    folder = apath/"gtex_info"
    folder.mkdir(exist_ok=True, parents=True)
    rerun = args.rerun

    generate_gtex_reference_plots(folder)



def generate_gtex_reference_plots(folder_path):
    """
    Reads stats.pkl and generates a dynamically scaling reference diagram.
    The canvas expands automatically based on the number of records present,
    and exports directly to unconstrained vector SVG format.
    """
    folder = Path(folder_path)
    pkl_path = folder / "stats.pkl"
    
    if not pkl_path.exists():
        raise FileNotFoundError(f"Missing required stats file: {pkl_path}")
        
    detailed_stats = load(pkl_path)
    
    summary_rows = []
    for abbr, stats in detailed_stats.items():
        summary_rows.append({
            "abbr": abbr,
            "n_sample": stats.get("n_samples", 0),
            "n_pc": stats.get("n_pc", 0),
            "n_tf": stats.get("n_tf", 0),
            "n_coreg": stats.get("n_coreg", 0),
            "n_reg": stats.get("n_reg", 0)
        })
        
    df = pd.DataFrame(summary_rows)
    df = df.sort_values(by="n_sample", ascending=False).reset_index(drop=True)
    
    # --- Dynamic Canvas Sizing System ---
    # Calculates size based on population count so future rows scale gracefully.
    num_datasets = len(df)
    calc_width = max(15, num_datasets * 0.45) 
    calc_height = 12.0 # Generous height allocation for multi-panel details
    
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 13,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'figure.titlesize': 16
    })
    
    fig,(ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(calc_width, calc_height), sharex=False)
    fig.suptitle("GTEx Tissue Dataset Profile", weight="bold", y=0.97)
    
    # Top Panel: Sample Sizes
    sns.barplot(
        data=df, 
        x="abbr", 
        y="n_sample", 
        hue="abbr",
        ax=ax1, 
        palette="Blues_r", 
        edgecolor="0.2",
        linewidth=0.5,
        legend=False
    )
    ax1.set_title("Total RNA-Seq Sample Size per Dataset", loc="left", weight="semibold")
    ax1.set_xlabel("")
    ax1.set_ylabel("Number of Samples")
    ax1.tick_params(axis='x', rotation=45)
    
    # Value annotations for sample bars (Linear scale placement)
    for p in ax1.patches:
        height = p.get_height()
        if height > 0:
            ax1.annotate(f'{int(height)}',
                         (p.get_x() + p.get_width() / 2., height),
                         ha='center', va='bottom', rotation=90, 
                         fontsize=7.5, xytext=(0, 4), textcoords='offset points')

    # Bottom Panel: Gene counts (PCG vs TF vs CoReg vs Reg)
    gene_df = df.melt(
        id_vars=["abbr"], 
        value_vars=["n_pc", "n_reg", "n_tf", "n_coreg"],
        var_name="Gene_Category", 
        value_name="Gene_Count"
    )
    
    label_mapping = {
        "n_pc": "Protein Coding Genes (PCG)",
        "n_tf": "Transcription Factors (TF)",
        "n_coreg": "Co-Regulators (CoReg)",
        "n_reg": "Regulators (Reg)"
    }
    gene_df["Gene_Category"] = gene_df["Gene_Category"].map(label_mapping)
    
    sns.barplot(
        data=gene_df, 
        x="abbr", 
        y="Gene_Count", 
        hue="Gene_Category", 
        ax=ax2,
        palette="muted",
        edgecolor="0.2",
        linewidth=0.5
    )
    
    ax2.set_title("Counts by dataset", loc="left", weight="semibold")
    ax2.set_xlabel("Dataset Abbreviation")
    ax2.set_ylabel("Gene Count")
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend(title="Feature Category", frameon=True, facecolor="white", edgecolor="0.8", loc="upper right")
    
    # Logarithmic transformation scales extreme gene counts evenly
    ax2.set_yscale('log')
    
    # Value annotations for grouped categories (Log scale placement)
    for p in ax2.patches:
        height = p.get_height()
        if height > 0:
            # 1.12 multiplicative factor positions labels uniformly on log coordinates
            ax2.annotate(f'{int(height)}',
                         (p.get_x() + p.get_width() / 2., height * 1.12),
                         ha='center', va='bottom', rotation=90, 
                         fontsize=6, color='black')
    
    # Provide upper cushion limit for the log annotations to prevent clipping
    ax2.set_ylim(bottom=1, top=ax2.get_ylim()[1] * 4)


    # --- Parse Sparsity Distributions tracking the sorted tissue order ---
    thresh_records = []
    for abbr in df["abbr"]:
        stats = detailed_stats.get(abbr, {})
        categories = {
            'PCG': stats.get('pc_sparsity_distribution', {}),
            'TF': stats.get('tf_sparsity_distribution', {}),
            'Regulators': stats.get('reg_sparsity_distribution', {})
        }
        for cat_name, dist_data in categories.items():
            if dist_data and 'threshold' in dist_data:
                for t, pct in zip(dist_data['threshold'], dist_data['percentage']):
                    thresh_records.append({
                        'abbr': abbr,
                        'Category': cat_name,
                        'Threshold': t,
                        'Percentage': pct
                    })
    
    thresh_df = pd.DataFrame(thresh_records)
    color_palette = sns.color_palette("turbo", n_colors=num_datasets)

    # Panel 3: PCG Threshold Distribution
    sns.lineplot(
        data=thresh_df[thresh_df['Category'] == 'PCG'],
        x='Threshold', y='Percentage', hue='abbr', ax=ax3,
        palette=color_palette, linewidth=1.5, alpha=0.85, legend=False
    )
    ax3.set_title("Protein Coding Genes (PCG) Retention Curve", loc="left", weight="semibold")
    ax3.set_xlabel("Sparsity Cutoff Limit (Fraction of Zeros)")
    ax3.set_ylabel("Genes Retained (%)")
    ax3.set_ylim(-5, 105)

    # Panel 4: TF Threshold Distribution
    sns.lineplot(
        data=thresh_df[thresh_df['Category'] == 'TF'],
        x='Threshold', y='Percentage', hue='abbr', ax=ax4,
        palette=color_palette, linewidth=1.5, alpha=0.85, legend=False
    )
    ax4.set_title("Transcription Factors (TF) Retention Curve", loc="left", weight="semibold")
    ax4.set_xlabel("Sparsity Cutoff Limit (Fraction of Zeros)")
    ax4.set_ylabel("Genes Retained (%)")
    ax4.set_ylim(-5, 105)

    # Panel 5: Regulators Threshold Distribution
    sns.lineplot(
        data=thresh_df[thresh_df['Category'] == 'Regulators'],
        x='Threshold', y='Percentage', hue='abbr', ax=ax5,
        palette=color_palette, linewidth=1.5, alpha=0.85, legend=True
    )
    ax5.set_title("Global Regulators (Reg) Retention Curve", loc="left", weight="semibold")
    ax5.set_xlabel("Sparsity Cutoff Limit (Fraction of Zeros)")
    ax5.set_ylabel("Genes Retained (%)")
    ax5.set_ylim(-5, 105)
    ax5.legend(title="Datasets", bbox_to_anchor=(1.01, 1.0), loc="upper left", borderaxespad=0, frameon=True)

    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.94])
    
    # Fixed file writer parameters to output authentic vector SVG markup
    plt.savefig(folder / "gtex_info.svg", format="svg", bbox_inches="tight")
    plt.close()

def generate_gtex_dcor_cache(input, env, options, args):
    """generate basic stats about GTEx tissue data"""
    rerun = args.rerun
    data_folder = Path(os.path.expandvars(env["DATA_PATH"])) / "dcorr_cache"
    data_folder.mkdir(exist_ok=True, parents=True)

    all_datasets = sorted(list(input["datasets"].keys()))
    print(all_datasets)

    selected_datasets = options.get("dataset_list",all_datasets)

    tf = get_tflist(env)
    coreg = get_coreglist(env)
    all_sources = list(set(tf) | set(coreg))

    for d in selected_datasets:
        print(d)
        cache_file_path = data_folder/f"{d}_pc.parquet"
        json_path = data_folder/f"{d}_metadata.json"
        if cache_file_path.exists() and not rerun:
            print("cache exists")
            continue
        
        dets = input["datasets"][d]
        dets["path"] = os.path.expandvars(f"$DATA_PATH/{dets['path']}")
        dets["convert_gene_names"]=True
        data = read_dataset(dets, env)
        pc_targets = get_protein_coding_genes(data.columns.tolist(),env)
        n = 10500
        sub = data[pc_targets]
        stats = pd.DataFrame({"m": sub.mean(), "v": sub.var()})
        pool = stats.sort_values(
            ["m", "v"], ascending=False).head(n).index.tolist()

        flat_dcor_cache, gene_to_idx, num_genes = dc.generate_combined_dcor_cache(
            data, all_sources, pool
        )

        print(f"Saving 1D cache to {cache_file_path}...")
        # Save compressed 1D matrix
        table = pa.Table.from_arrays([pa.array(flat_dcor_cache)], names=['dcor'])
        pq.write_table(table, cache_file_path, compression='SNAPPY')
        
        print(f"Saving gene metadata to {json_path}...")
        # Save JSON metadata map
        metadata = {
            "num_genes": num_genes,
            "gene_to_idx": gene_to_idx
        }
        with open(json_path, 'w') as f:
            json.dump(metadata, f, indent=2)


