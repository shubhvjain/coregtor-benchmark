#new code, metadata part is missing 
#!/usr/bin/env python3
"""
Download datasets defined in datasets.json.

Usage:
    python dataset.py --all
    python dataset.py --id amygdala
    python dataset.py --list
    python dataset.py --all --rerun
"""

import json
import os
import argparse
import sys
from pathlib import Path
import src.results.util as ut
import pooch

from tfitpy import get_gene_products


def load_list() -> list[dict]:
    path = Path(__file__).parent.parent.parent /"experiments"/"analysis"/"gtex.json"
    with open(path) as f:
        return json.load(f)["datasets"]


def get_data_root(env) -> Path:
    data_path = env["DATA_PATH"]
    if not data_path:
        raise EnvironmentError("DATA_PATH env var not set")
    return Path(data_path)


def is_downloaded(data_root: Path,dataset_path: str) -> bool:
    d = data_root / dataset_path
    return d.exists()



def download_dataset(dataset_id, dataset, data_root, rerun=False):
    # 1. Setup directories
    dataset_dir = data_root / dataset_id
    dataset_dir.mkdir(parents=True, exist_ok=True)
    
    # Target file path (e.g., data_root/SKL/gene_tpm_v11_skin_sun_exposed_lower_leg.gct)
    # Using Path(dataset["path"]).name ensures we just get the filename itself
    final_filename = Path(dataset["path"]).name
    final_file_path = dataset_dir / final_filename
    
    # 2. Check if already downloaded
    if not rerun and final_file_path.exists():
        print(f"skip {dataset_id} - already exists at {final_file_path}")
        return

    url = dataset["url"]
    
    # 3. Setup the decompression processor
    # We supply the explicit final name we want (without .gz)
    if url.endswith(".gz"):
        processor = pooch.Decompress(name=final_filename)
    elif url.endswith(".zip"):
        processor = pooch.Unzip()
    elif ".tar" in url:
        processor = pooch.Untar()
    else:
        processor = None

    # 4. Download and unpack
    print(f"Downloading and unpacking {dataset_id}...")
    pooch.retrieve(
        url=url,
        fname=Path(url).name,  # Saves temporary archive as the original URL name
        path=dataset_dir,
        progressbar=True,
        processor=processor,
    )
    
    # 5. Clean up the compressed file
    # Pooch downloads the .gz to dataset_dir / Path(url).name
    compressed_file_path = dataset_dir / Path(url).name
    if compressed_file_path.exists() and url.endswith(".gz"):
        compressed_file_path.unlink()
        print(f"Removed compressed archive: {compressed_file_path.name}")


def download_all(env,rerun: bool = False):
    data_root = get_data_root(env)
    datasets = load_list()
    ok, failed = 0, []
    for ds in datasets.keys():
        try:
            id = ds
            dets = datasets[ds]
            download_dataset(id,dets, data_root, rerun=rerun)
            ok += 1
        except Exception as e:
            print(f"failed {ds['id']}: {e}")
            failed.append(ds["id"])
    print(f"\n{ok}/{len(datasets)} successful")
    if failed:
        print(f"failed: {', '.join(failed)}")


def download_one(env,dataset_id: str, rerun: bool = False):
    data_root = get_data_root(env)
    datasets = load_list()
    ds =  datasets[dataset_id] #next((d for d in datasets if d["id"] == dataset_id), None)
    if ds is None:
        print(f"'{dataset_id}' not found. available: {', '.join(d['id'] for d in datasets)}")
        sys.exit(1)
    download_dataset(dataset_id,ds, data_root, rerun=rerun)


def list_datasets(env):
    data_root = get_data_root(env)
    datasets = load_list()
    for ds in datasets:
        status = "ok" if is_downloaded(ds["id"], data_root) else "missing"
        print(f"{ds['id']:<25} {status}")


def generate_gene_lists(env):
    """"""
    path = Path(__file__).parent.parent.parent/"analysis"/"gene_lists.json"
    with open(path) as f:
        data = json.load(f)
    for d in data["list"]:
        dataset_dir = Path(env["DATA_PATH"]) / d["title"]
        dataset_dir.mkdir(parents=True, exist_ok=True)
        list_file = dataset_dir/"list.csv"
        df = get_gene_products(env["DATA_PATH"],d["term"])
        df.to_csv(list_file,index=False)



def main():
    parser = argparse.ArgumentParser(description="Download datasets")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all",  action="store_true")
    group.add_argument("--id",   metavar="ID")
    group.add_argument("--list", action="store_true")
    group.add_argument("--gene_list", action="store_true")
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--env",   required=True,
                        metavar="PATH", help="Path to the .env file.")
    args = parser.parse_args()

    env = ut.get_env(args.env)
    if args.all:
        download_all(env,rerun=args.rerun)
    elif args.id:
        download_one(env,args.id, rerun=args.rerun)
    elif args.list:
        list_datasets()
    elif args.gene_list:
        generate_gene_lists(env)


if __name__ == "__main__":
    main()