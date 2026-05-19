import json
import sys
from pathlib import Path
from dotenv import dotenv_values
import sqlite3
from joblib import  load
import os
import pandas as pd

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
    output_path = env.EXP_OUTPUT_PATH/input
    r
    """
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
