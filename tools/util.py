import pandas as pd
import numpy as np
import json
import sys
import os
import secrets
from pathlib import Path
import argparse
from dotenv import dotenv_values
import sqlite3
from datetime import datetime

import traceback


def get_env():
    """this code is run within the container where volumes are mounted as fixed locations"""
    env = {
        "DATA_PATH": Path("/app/dataset"),
        "EXP_TEMP_PATH": Path("/app/temp"),
        "EXP_OUTPUT_PATH": Path("/app/output"),
        "EXP_INPUT_PATH": Path("/app/input")
    }
    return env


def get_input(input_path):
    path = Path(input_path)
    if not path.exists():
        print(f"[error] input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


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


def setup_experiment(exp, config):
    setup_allowed = exp.get("type") in ["run_coregnet","run_rtnduals"]
    o,t = get_exp_path(exp,config)
    if not  setup_allowed:
        print("invalid tool type")
        return
    # Construct the full path to the experiment folder
    exp_dir = o 
    temp_dir = t
    done_file = temp_dir / "init_done.txt"
    input_file = exp_dir / "input.json"

    # Check if already initialized
    if done_file.exists():
        print("Exp folder already exists and is initialized")
        return

    """
    1. load the data 
    2. compute source and target list
    3. setup the folder structure
    4. create the input.json file in output folder
    """
    print(config)
    # Create folders if it doesn't exist
    exp_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    dataset = read_dataset(exp["dataset"], config)
    tf_list = get_tflist(config)
    #print(tf_list)
    input = {**exp}
    # source_list
    
    input["source_genes"] = get_filtered_genes(
        dataset, exp.get("sources"), tf_list,config)
    input["target_genes"] = get_filtered_genes(
        dataset, exp.get("targets"), tf_list,config)

    with open(input_file, "w") as f:
        json.dump(input, f)


    Path(done_file).write_text(f"Init at {datetime.now()}")

    return



def start_experiment(exp_data,env):
    """
    check if input.json exists inside the output folder
    reads the dataset and returns both
    """
    out,temp = get_exp_path(exp_data,env)
    input_file_path = out/"input.json"
    CONFIG = get_env()
    # print(CONFIG)
    if not input_file_path.exists():
        print("Experiment not yet setup. Initializing")
        setup_experiment(exp_data,CONFIG)
    
    with open(input_file_path) as f:
        input_data = json.load(f)

    
    os.environ.setdefault("DATA_PATH", env.get("DATA_PATH") )
    #print(CONFIG)
    dataset = read_dataset(exp_data["dataset"],CONFIG)
    #print(111)
    return input_data, dataset



def read_dataset(details, config):
    """"""
    dtype = details.get("type", None)
    if dtype is None:
        raise ValueError("invalid dataset.type")

    #print(config.get("data_path"))
    if dtype == "gct":
        return read_gct(details.get("path"), config, convert_gene_names=details.get("convert_gene_names", True))


def read_gct(file_path1, CONFIG=None, convert_gene_names=False):
    """
    Read Gene Cluster Text (GCT) format file into a pandas DataFrame.

    GCT is a tab-delimited format which include
    - Line 1: Version information
    - Line 2: Dimensions (genes x samples)  
    - Line 3+: Header with Name, Description, and sample columns
    - Data rows: Gene information and expression values 

    Assuming there is  "Description" column that has the name of genes have the gene names for each gene row.

    Args:
        file_path (str or Path) : Path to the GCT file.

    Returns:
        pd.DataFrame :  DataFrame with genes as rows and samples as columns. The index is gene_name (from Description column), columns are sample identifiers. Each cell had gene expression levels

    Notes:
    ------
    - Removes the 'Name' column (gene IDs) and uses 'Description' as gene names

    """
    
    file_path = Path(os.path.expandvars(str(file_path1)))
    #print(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    try:
        # Read GCT file, skipping version and dimension lines
        df = pd.read_csv(file_path, skiprows=2, sep="\t")
        #print(df)
        gene = df["Name"].values.tolist()
        if convert_gene_names:
            if CONFIG is None:
                raise ValueError("DATA_PATH not provided")
            mps = get_mappings(CONFIG, gene_list=gene, source='gene_id',
                               target='gene_name')
            df["Name"] = df["Name"].map(mps)
        if "Name" not in df.columns or "Description" not in df.columns:
            raise ValueError(
                "GCT file must contain 'Name' and 'Description' columns")
        # remove Name column, rename Description to gene_name, set as index
        df = df.drop(columns=["Name"]).rename(
            columns={"Description": "gene_name"})
        df = df.set_index("gene_name")
        df = df.transpose().rename_axis("sample_name")
        print(1111)
        return df
    except Exception as e:
        print(e)
        traceback.print_exc()
        raise ValueError(f"Error reading GCT file {file_path}: {str(e)}")




def get_mappings(CONFIG, gene_list, source, target, batch_size=900):
    """
    Map gene identifiers from source format to target format.
    Assuming tfitpy setup was run and the source db file is available.
    Args:
        gene_list: List of identifiers to map
        source: Source attribute name in database (e.g., 'gene_name', 'gene_id')
        target: Target attribute name in database (e.g., 'gene_id', 'gene_name')
        batch_size: Number of items to process per query (default 900, under SQLite's 999 limit)
    Returns:
        Dictionary mapping {source_value: target_value, ...}
    """
    db_path = Path(CONFIG["DATA_PATH"])/"gencode"/"gene_name_mapping.db"
    #print(db_path.exists())
    con = sqlite3.connect(db_path)

    try:
        # Remove None/NaN and get unique values
        unique_values = [x for x in set(
            gene_list) if x is not None and pd.notna(x)]
        if len(unique_values) == 0:
            print("No values to map")
            return {}

        # Process in batches to avoid SQLite variable limit
        mapping_dict = {}
        for i in range(0, len(unique_values), batch_size):
            batch = unique_values[i:i + batch_size]
            placeholders = ','.join(['?'] * len(batch))
            query = f"""
                SELECT DISTINCT {source}, {target}
                FROM mappings
                WHERE {source} IN ({placeholders})
                AND {target} IS NOT NULL
            """

            # Execute query and update mapping dict
            mapping_df = pd.read_sql_query(query, con, params=batch)
            mapping_dict.update(
                dict(zip(mapping_df[source], mapping_df[target])))

        # Report statistics
        total = len(gene_list)
        unique_count = len(unique_values)
        mapped_count = len(mapping_dict)
        failed_count = unique_count - mapped_count
        null_in_results = sum(
            1 for gene in gene_list if gene not in mapping_dict)

        print(f"Mapping from {source} to {target}:")
        print(f"  Total input values: {total}")
        print(f"  Unique values: {unique_count}")
        print(f"  Successfully mapped: {mapped_count}")
        print(f"  Failed to map: {failed_count}")
        print(
            f"  Null results: {null_in_results} ({null_in_results/total*100:.2f}%)")
        return mapping_dict

    finally:
        con.close()

def get_tflist(CONFIG):
    """
    """
    tf_path = CONFIG["DATA_PATH"] / "tflist"/ "allTFs_hg38.txt"
    df = pd.read_csv(tf_path, names=["gene_name"], header=None)
    return  df["gene_name"].tolist()


def get_protein_coding_genes(gene_list, CONFIG):
    """
    Fetches all protein-coding genes once and intersects with input list in memory.
    Best for one-time calls with large input lists.
    """
    db_path = CONFIG["DATA_PATH"] / "gencode" / "gene_name_mapping.db"
    
    # 1. Fetch all protein-coding genes from the DB
    query = "SELECT DISTINCT gene_name FROM mappings WHERE gene_type = 'protein_coding'"
    
    try:
        with sqlite3.connect(db_path) as con:
            cursor = con.execute(query)
            pc_reference_set = {row[0] for row in cursor.fetchall()}
    except sqlite3.OperationalError as e:
        print(f"Database error: {e}")
        return []

    # 2. Use set intersection to find matches
    gene_list_set = set(gene_list)
    protein_coding = list(pc_reference_set.intersection(gene_list_set))
    
    return protein_coding

def get_filtered_genes(df, options, tflist, config):
    """
    Selects and filters genes from df based on type and stats.
    """
    #print(options)
    g_type = options.get("type", "all")

    # 1. Selection
    if g_type == "tf":
        pool = [g for g in tflist if g in df.columns]
    elif g_type == "no_tf":
        pool = [g for g in df.columns if g not in tflist]
    elif g_type == "pc_only":
        pool = get_protein_coding_genes(df.columns.tolist(),config)
    elif g_type == "custom":
        pool = [g for g in options.get("items", []) if g in df.columns]
    else:
        raise ValueError("Invalid selection type provided")

    if not pool:
        return []

    # 2. Thresholding (Non-zero occupancy)
    nz_thresh = options.get("non_zero_threshold")
    if nz_thresh is not None:
        occupancy = (df[pool] > 0).mean()
        pool = occupancy[occupancy >= nz_thresh].index.tolist()

    # 3. Ranking (Top N by Mean/Var)
    n = options.get("top_expressed_n")
    if n is not None and len(pool) > n:
        sub = df[pool]
        stats = pd.DataFrame({"m": sub.mean(), "v": sub.var()})
        pool = stats.sort_values(
            ["m", "v"], ascending=False).head(n).index.tolist()

    return pool
