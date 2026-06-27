import json
import sys
from pathlib import Path
from dotenv import dotenv_values
import sqlite3
from joblib import  load
import os
import pandas as pd
import numpy as np

def get_env(env_path):
    path = Path(env_path)
    if not path.exists():
        print(f"[error] env file not found: {env_path}", file=sys.stderr)
        sys.exit(1)
    values = dict(dotenv_values(path))    
    os.environ.update(values)
    return values


def get_input(input_path):
    path = Path(input_path)
    if not path.exists():
        print(f"[error] input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def get_result_by_id(input, id):
    results = input.get("results", [])
    for result in results:
        if result.get("id") == id:
            return result
    print(f"[error] no result with id '{id}' found in input", file=sys.stderr)
    sys.exit(1)
 
def get_exp_path(input,env):
    """
    Use this everywhere to get the path to the experiment results and temp folder. 
    output_path = env.EXP_OUTPUT_PATH/input
    r
    """
    if input.get("path",None) is not  None:
        output_path = Path(env.get("EXP_OUTPUT_PATH"))/ input.get("path")
        temp_path = Path(env.get("EXP_TEMP_PATH"))/ input.get("path")
    else:
        output_path = Path(env.get("EXP_OUTPUT_PATH"))/ input.get("id")
        temp_path = Path(env.get("EXP_TEMP_PATH"))/ input.get("id")
    return output_path,temp_path

def get_exp_target_list(out_path):
    """
    check the status.db file in output folder.
    """
    db = out_path / "status.db"
    conn = sqlite3.connect(db, timeout=30)
    targets = conn.execute(
        "SELECT gene FROM genes WHERE status='done' ORDER BY gene"
    ).fetchall()
    conn.close()
    target_list = [r[0] for r in targets]
    return target_list

def get_temp_file(temp_path,target):
    """"""
    file_path =  temp_path/"results"/f"{target}.pkl"
    data = load(file_path)
    return data 


def get_analysis_path(env):
    """"""
    # print(env)
    p = Path(env["ANALYSIS_OUTPUT_PATH"])
    return p

def get_exp_gtex_info(dataset_id):
    """""" 
    path = Path(__file__).parent.parent.parent /"analysis"/"gtex_info"/"stats.csv"
    d = pd.read_csv(path)
    dets = d[d["abbr"]==dataset_id]
    return dets.to_dict(orient='records')[0]

def get_gtex_key(include_sample_count=False):
    """
    return key:value for dataset labels 
    """
    path = Path(__file__).parent.parent.parent /"analysis"/"gtex_info"/"stats.csv"
    d = pd.read_csv(path)
    dbkey = {}
    for _, row in d.iterrows():
        ky = str(row["abbr"])
        label = str(row["name"])
        if include_sample_count:
            label = f"{label} ( {row["n_sample"]} samples)"
        dbkey[ky] = label

    return dbkey


# 15-Color Palette curated for plotting and data visualization
color_palette = {
    # Deep/Dark tones for text, axes, and dominant baselines
    "dark_navy": "#23233B",
    "deep_slate": "#2C4268",
    "midnight_blue": "#023047",
    
    # Vibrant primary colors for main data series/lines
    "classic_blue": "#007BBA",
    "vibrant_cyan": "#00A9E2",
    "ocean_blue": "#126782",
    "medium_blue": "#219EBC",
    
    # Light tones for backgrounds, secondary bars, or grids
    "sky_blue": "#73BFDC",
    "pale_blue": "#8ECAE6",
    "ice_blue": "#BCE3FA",
    
    # High-contrast accent colors for callouts or highlights
    "bright_yellow": "#FFB703",
    "gold_orange": "#FD9E02",
    "vibrant_orange": "#FB8500",
    
    # Neutrals for benchmarks, targets, or secondary data
    "medium_gray": "#9B9C9B",
    "light_gray": "#EBEBEB"
}

default_sns_configs = {

    "axes.titlesize": 8,       # Sets all subplot titles globally
    "axes.labelsize": 6,       # Sets all x and y axis label sizes globally
    "xtick.labelsize": 6,      # Sets all x-axis tick sizes globally
    "ytick.labelsize": 6,      # Sets all y-axis tick sizes globally

    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlepad": 8,
    "svg.fonttype": "none",
    "font.family": "sans-serif"
}


index_aggregated = [
    {
        "key":"PPI(H)",
        "cols":["shared_PPI_partners_score_hippie","shortest_PPI_path_score_hippie","lcc_score_hippie","tc_score_hippie","density_score_hippie"],
        "label":"Overall PPI Score (Hippie)"
    },
    {
        "key":"PPI(B)",
        "cols":["shared_PPI_partners_score_biogrid","shortest_PPI_path_score_biogrid","lcc_score_biogrid","tc_score_biogrid","density_score_biogrid"],
        "label":"Overall PPI Score (Biogrid)"
    },
     {
        "key":"PPI(S)",
        "cols":["shared_PPI_partners_score_stringdb","shortest_PPI_path_score_stringdb","lcc_score_stringdb","tc_score_stringdb","density_score_stringdb"],
        "label":"Overall PPI Score (Stringdb)"
    },
    {
        "key":"GO",
        "cols":["goa_similarity_jc","goa_similarity_lin","goa_similarity_resnik"],
        "label":"Overall GO Score"
    },
    {
        "key":"DC",
        "cols":["DC1","DC2"],
        "label":"Overall Distance Correlation Score"
    },
    {
        "key":"TFBS",
        "cols":["TFBS_affinity_score"],
        "label":"TFBS Affinity Score"
    }
]

index_labels = {
    "PPI1(H)":"PPI Shared Partners Score (Hippie)",
    "PPI1(B)":"PPI Shared Partners Score (Biogrid)",
    "PPI1(S)":"PPI Shared Partners Score (Stringdb)",
    "PPI2(H)":"PPI Shortest Path Score (Hippie)",
    "PPI2(B)":"PPI Shortest Path Score (Biogrid)",
    "PPI2(S)":"PPI Shortest Path Score (Stringdb)",
    "PPI3(H)":"PPI Induced Subgraph Density Score (Hippie)",
    "PPI3(B)":"PPI Induced Subgraph Density Score (Biogrid)",
    "PPI3(S)":"PPI Induced Subgraph Density Score (Stringdb)",
    "PPI4(H)":"PPI Induced Subgraph Largest Connected Component Score (Hippie)",
    "PPI4(B)":"PPI Induced Subgraph Largest Connected Component Score (Biogrid)",
    "PPI4(S)":"PPI Induced Subgraph Largest Connected Component Score (Stringdb)",
    "PPI5(H)":"PPI Induced Subgraph Target Connectivity Score (Hippie)",
    "PPI5(B)":"PPI Induced Subgraph Target Connectivity Score (Biogrid)",
    "PPI5(S)":"PPI Induced Subgraph Target Connectivity Score (Stringdb)",
    "GO1":"GO Based Semantic Similarity (Jiang and Conrath)",
    "GO2":"GO Based Semantic Similarity (Lin)",
    "GO3":"GO Based Semantic Similarity (Resnik)",
    "DC1":"Source Source Distance Correlation Score",
    "DC2":"Source Target Distance Correlation Score",
    "TFBS1":"Transcription Factor Binding Affinity Score",
    "CTEC":"Common Targets Expression Coherence Score ",
    "CTPPI":"Common Targets PPI Coherence Score",
    "CTFE":"Common Targets Function Enrichment Score",
    "PPI(H)":"Overall PPI Score (Hippie)",
    "PPI(B)":"Overall PPI Score (Biogrid)",
    "PPI(S)":"Overall PPI Score (Stringdb)",
    "GO":"Overall GO Score",
    "DC":"Overall Distance Correlation Score",
    "TFBS":"TFBS Affinity Score"
}

col_index_map = {
    "shared_PPI_partners_score_hippie":"PPI1(H)",
    "shared_PPI_partners_score_biogrid":"PPI1(B)",
    "shared_PPI_partners_score_stringdb":"PPI1(S)",
    "shortest_PPI_path_score_hippie":"PPI2(H)",
    "shortest_PPI_path_score_stringdb":"PPI2(S)",
    "shortest_PPI_path_score_biogrid":"PPI2(B)",
    "density_score_hippie":"PPI3(H)",
    "lcc_score_hippie":"PPI4(H)",
    "tc_score_hippie":"PPI5(H)",
    "density_score_stringdb":"PPI3(S)",
    "lcc_score_stringdb":"PPI4(S)",
    "tc_score_stringdb":"PPI5(S)",
    "density_score_biogrid":"PPI3(B)",
    "lcc_score_biogrid":"PPI4(B)",
    "tc_score_biogrid":"PPI5(B)",
    "goa_similarity_jc":"GO1",
    "goa_similarity_lin":"GO2",
    "goa_similarity_resnik":"GO3",
    "DC1":"DC1",
    "DC2":"DC2",
    "TFBS_affinity_score":"TFBS1",
    "TFBS_affinity_sum_score":"TFBS2",
    "dCor_sources_score":"DC1",
    "dCor_target_score":"DC2"   
}

distance_map = {
        "euclidean_distance": "EU",
        "manhattan_distance": "MH",
        "canberra_distance": "CB",
        "chebyshev_distance": "CH",
        "cosine_distance": "CS",
        "weighted_jaccard": "WJ",
        "sorensen_distance": "SR",
        "jensenshannon_distance": "JS",
}

col_freq = {
    "dCor_targets_score":"CTEC",
    "ppi_shared_partners_hippie":"CTPPI",
    "ppi_shared_partners_stringdb":"CTPPI",
    "ppi_shared_partners_biogrid":"CTPPI",
    "go_score":"CTFE"
}

index_col1 = [
    "shared_PPI_partners_score_biogrid",
    "shortest_PPI_path_score_biogrid",
    "density_score_biogrid",
    "lcc_score_biogrid",
    "tc_score_biogrid",
    "goa_similarity_jc",
    "goa_similarity_lin",
    "goa_similarity_resnik",
    "TFBS_affinity_score",
    "TFBS_affinity_sum_score",
    "dCor_sources_score",
    "dCor_target_score"
]

index_col2 = [
    "shared_PPI_partners_score_stringdb",
    "shortest_PPI_path_score_stringdb",
    "density_score_stringdb",
    "lcc_score_stringdb",
    "tc_score_stringdb",
    "goa_similarity_jc",
    "goa_similarity_lin",
    "goa_similarity_resnik",
    "TFBS_affinity_score",
    "TFBS_affinity_sum_score",
    "dCor_sources_score",
    "dCor_target_score"
]

index_col3 = [
    "shared_PPI_partners_score_hippie",
    "shortest_PPI_path_score_hippie",
    "density_score_hippie",
    "lcc_score_hippie",
    "tc_score_hippie",
    "goa_similarity_jc",
    "goa_similarity_lin",
    "goa_similarity_resnik",
    "TFBS_affinity_score",
    "TFBS_affinity_sum_score",
    "dCor_sources_score",
    "dCor_target_score"
]


tool_map = {
    "coregnet": {
        "color": color_palette["bright_yellow"],
        "label": "CoRegNet",
    },
    "coregtor": {
        "color": color_palette["classic_blue"],
        "label": "CoReGToR",
    },
}


freq_index_aggregated = [
    {
        "key":"CTS",
        "cols":["dCor_targets_score","ppi_shared_partners_stringdb","go_score"],
        "label":"Overall Common Targets Score"
    }
]


def process_and_rename_dataset(df, config_df):
    """
    1. Computes the custom aggregated scores (Sums instead of Means).
    2. Renames columns using the standard col_index_map shorthand tokens.
    3. Adds short-hand 'distance_measure' column based on config mapping.
    4. Dynamically appends 3 layers of clustering method details.
    5. Returns both the processed dataframe and the updated clustering labels dictionary.
    """
    processed_df = df.copy()
    
    # Step 1: Add Aggregations first while raw data column keys exist
    for agg in index_aggregated:
        key = agg["key"]
        cols_to_sum = agg["cols"]
        valid_cols = [c for c in cols_to_sum if c in processed_df.columns]
        if valid_cols:
            processed_df[key] = processed_df[valid_cols].sum(axis=1, min_count=1)
        else:
            processed_df[key] = np.nan
            
    # Step 2: Rename all remaining raw columns to your shorthand map
    processed_df = processed_df.rename(columns=col_index_map)

    # add a total score 
    all_cols = [ i['key'] for i in index_aggregated]
    processed_df["OVERALL_SCORE"] = processed_df[all_cols].sum(axis=1, min_count=1)
    # --- Pre-processing for Step 3 & 4 ---
    config_dict = config_df.set_index("config_Id").to_dict(orient="index")
    
    # Base mapping dictionary
    dynamic_clustering_labels = {
        'HC': 'Hierarchical Clustering',
        'CD': 'Community Detection',
        'HB': 'HDBSCAN',
        'HC1': 'Hierarchical (Inconsistency)',
        'HC2': 'Hierarchical (One Cluster 10)',
        'CD_TH': 'Community Detection (Threshold)',
        'CD_KNN': 'Community Detection (KNN)',
        'HB1': 'HDBSCAN Default'
    }

    # Isolate community detection configurations to dynamically assign CL tokens
    cd_mask = config_df["clustering_method"] == "community_detection"
    unique_cd_combos = (
        config_df[cd_mask][["clustering_option1", "clustering_option2"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    
    combo_to_cl_id = {}
    for idx, row in enumerate(unique_cd_combos.itertuples(index=False), start=1):
        cl_id = f"CL{idx}"
        combo_to_cl_id[(row.clustering_option1, row.clustering_option2)] = cl_id
        
        # Clean up labels for map (e.g. "threshold_percentile_50" -> "Threshold Percentile 50")
        opt1_clean = row.clustering_option1.replace('_', ' ').title()
        dynamic_clustering_labels[cl_id] = f"{opt1_clean} (r={row.clustering_option2})"

    # Step 3: Map Distance Measures
    raw_distances = processed_df["config_name"].map(
        lambda x: config_dict.get(x, {}).get("distance_measure")
    )
    processed_df["distance_measure"] = raw_distances.map(distance_map).fillna(raw_distances)

    # --- Step 4: Map Clustering Methods ---
    def get_cfg(cfg_name):
        return config_dict.get(cfg_name, {})

    # Column 1
    method1_map = {"hierarchical": "HC", "community_detection": "CD", "hdbscan": "HB"}
    processed_df["clustering_method1"] = processed_df["config_name"].map(
        lambda x: method1_map.get(get_cfg(x).get("clustering_method"))
    )

    # Column 2
    def map_method2(cfg_name):
        cfg = get_cfg(cfg_name)
        method = cfg.get("clustering_method")
        opt1 = str(cfg.get("clustering_option1"))
        
        if method == "hierarchical":
            return "HC1" if "inconsistency" in opt1 else "HC2"
        elif method == "community_detection":
            return "CD_TH" if "threshold" in opt1 else "CD_KNN"
        elif method == "hdbscan":
            return "HB1"
        return None

    processed_df["clustering_method2"] = processed_df["config_name"].map(map_method2)

    # Column 3
    def map_method3(cfg_name):
        cfg = get_cfg(cfg_name)
        if cfg.get("clustering_method") == "community_detection":
            key = (cfg.get("clustering_option1"), cfg.get("clustering_option2"))
            return combo_to_cl_id.get(key)
        return None

    processed_df["clustering_method3"] = processed_df["config_name"].map(map_method3)

    return processed_df, dynamic_clustering_labels



def process_and_rename_freq(df, config_df):
    """
    1. Computes the custom aggregated scores (Sums instead of Means).
    2. Renames columns using the standard col_index_map shorthand tokens.
    3. Adds short-hand 'distance_measure' column based on config mapping.
    4. Dynamically appends 3 layers of clustering method details.
    5. Returns both the processed dataframe and the updated clustering labels dictionary.
    """
    processed_df = df.copy()
    
    # Step 1: Add Aggregations first while raw data column keys exist
    for agg in freq_index_aggregated:
        key = agg["key"]
        cols_to_sum = agg["cols"]
        valid_cols = [c for c in cols_to_sum if c in processed_df.columns]
        if valid_cols:
            processed_df[key] = processed_df[valid_cols].sum(axis=1, min_count=1)
        else:
            processed_df[key] = np.nan
            
    # Step 2: Rename all remaining raw columns to your shorthand map
    processed_df = processed_df.rename(columns=col_freq)

    # --- Pre-processing for Step 3 & 4 ---
    config_dict = config_df.set_index("config_Id").to_dict(orient="index")
    
    # Base mapping dictionary
    dynamic_clustering_labels = {
        'HC': 'Hierarchical Clustering',
        'CD': 'Community Detection',
        'HB': 'HDBSCAN',
        'HC1': 'Hierarchical (Inconsistency)',
        'HC2': 'Hierarchical (One Cluster 10)',
        'CD_TH': 'Community Detection (Threshold)',
        'CD_KNN': 'Community Detection (KNN)',
        'HB1': 'HDBSCAN Default'
    }

    # Isolate community detection configurations to dynamically assign CL tokens
    cd_mask = config_df["clustering_method"] == "community_detection"
    unique_cd_combos = (
        config_df[cd_mask][["clustering_option1", "clustering_option2"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    
    combo_to_cl_id = {}
    for idx, row in enumerate(unique_cd_combos.itertuples(index=False), start=1):
        cl_id = f"CL{idx}"
        combo_to_cl_id[(row.clustering_option1, row.clustering_option2)] = cl_id
        
        # Clean up labels for map (e.g. "threshold_percentile_50" -> "Threshold Percentile 50")
        opt1_clean = row.clustering_option1.replace('_', ' ').title()
        dynamic_clustering_labels[cl_id] = f"{opt1_clean} (r={row.clustering_option2})"

    # Step 3: Map Distance Measures
    raw_distances = processed_df["config_name"].map(
        lambda x: config_dict.get(x, {}).get("distance_measure")
    )
    processed_df["distance_measure"] = raw_distances.map(distance_map).fillna(raw_distances)

    # --- Step 4: Map Clustering Methods ---
    def get_cfg(cfg_name):
        return config_dict.get(cfg_name, {})

    # Column 1
    method1_map = {"hierarchical": "HC", "community_detection": "CD", "hdbscan": "HB"}
    processed_df["clustering_method1"] = processed_df["config_name"].map(
        lambda x: method1_map.get(get_cfg(x).get("clustering_method"))
    )

    # Column 2
    def map_method2(cfg_name):
        cfg = get_cfg(cfg_name)
        method = cfg.get("clustering_method")
        opt1 = str(cfg.get("clustering_option1"))
        
        if method == "hierarchical":
            return "HC1" if "inconsistency" in opt1 else "HC2"
        elif method == "community_detection":
            return "CD_TH" if "threshold" in opt1 else "CD_KNN"
        elif method == "hdbscan":
            return "HB1"
        return None

    processed_df["clustering_method2"] = processed_df["config_name"].map(map_method2)

    # Column 3
    def map_method3(cfg_name):
        cfg = get_cfg(cfg_name)
        if cfg.get("clustering_method") == "community_detection":
            key = (cfg.get("clustering_option1"), cfg.get("clustering_option2"))
            return combo_to_cl_id.get(key)
        return None

    processed_df["clustering_method3"] = processed_df["config_name"].map(map_method3)

    return processed_df, dynamic_clustering_labels



def get_cluster_list_result(input_data, result_name):
    """
    Looks for a result with "type": "generate_result_file" and a matching 
    "result_name", then returns its "cluster_list". Raises ValueError if not found.
    """
    # Safeguard in case the "results" key doesn't exist or is empty
    results = input_data.get("results", [])
    
    for result in results:
        if result.get("type") == "generate_result_file" and result.get("result_name") == result_name:
            # Return the cluster list if found
            return result.get("cluster_list", [])
            
    # If the loop finishes without returning, raise the error
    raise ValueError(f"No result found with type 'generate_result_file' and result_name '{result_name}'")


def get_ps_configs(config_name="ps_config"):
    """returns the parameter selection configs by reading the file in the experiments folder"""
    pconfig_path = Path(__file__).parent.parent.parent /"experiments"/"analysis"/f"{config_name}.csv"
    pconfig = pd.read_csv(pconfig_path)
    return pconfig


def get_input_source_list(output_path):
    """"""
    input_data = get_input(output_path/"input.json")
    source_list = input_data["source_genes"]
    target_list = input_data["target_genes"]
    return source_list,target_list
