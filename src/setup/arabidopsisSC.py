"""
Code to setup the Arabidopsis Single Cell dataset
Denyer et al , Spatiotemporal Developmental Trajectories in the Arabidopsis Root Revealed Using High-Throughput Single-Cell RNA Sequencing  (https://doi.org/10.1016/j.devcel.2019.02.022)
"""
from pathlib import Path
import argparse
import os

import pooch
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

import src.results.util as ut


MARKER_GENES = {
    "Columella": [
        "AT2G33230","AT3G10320","AT4G18350","AT5G22550","AT1G26680",
        "AT2G18230","AT3G18230","AT3G28540","AT3G52890","AT3G55180",
        "AT5G07550","AT3G29810","AT5G48300","AT5G64860","AT3G44240",
        "AT3G62320","AT5G48130","AT1G04600","AT5G58580"
    ],
    "QC": [
        "AT1G68640","AT3G26120","AT5G23780","AT5G17430","AT2G03830",
        "AT2G23050","AT1G30530","AT5G25900","AT5G45780","AT1G11420",
        "AT3G13510","AT3G30350","AT1G16070","AT1G75580"
    ],
    "Cortex": [
        "AT1G05570","AT1G15460","AT1G62510","AT3G12700","AT3G14850",
        "AT5G03570","AT5G44130","AT5G53370","AT5G55250","AT1G64380",
        "AT5G48000","AT3G21670","AT1G30380","AT1G52410","AT3G26300",
        "AT4G00700","AT3G56080","AT5G27350"
    ],
    "Endodermis": [
        "AT1G60690","AT1G61590","AT1G71740","AT2G21100","AT2G27370",
        "AT2G30210","AT2G39430","AT2G40160","AT2G48130","AT3G11550",
        "AT3G22620","AT4G02090","AT4G17215","AT5G06200","AT5G42180",
        "AT5G66390","AT1G01120","AT2G14900"
    ],
    "Pericycle": [
        "AT1G74430","AT3G24310","AT4G28250","AT5G25920","AT4G29770",
        "AT1G58025","AT4G34650","AT5G06080","AT1G52400","AT2G19900",
        "AT5G24240","AT3G06020","AT3G29635","AT4G17760"
    ],
    "Protophloem": [
        "AT1G14730","AT5G40260","AT5G57130","AT5G62940","AT1G15080",
        "AT3G01680","AT4G13600","AT5G67460","AT5G22870","AT1G54330",
        "AT1G60030","AT3G17730","AT4G17910","AT4G36620","AT4G03270",
        "AT1G06490","AT1G11915","AT1G13050","AT1G31335","AT1G33480"
    ],
    "Phloem": [
        "AT1G11570","AT1G69580","AT1G69970","AT2G41200","AT3G47180",
        "AT4G04620","AT4G16400","AT5G01480","AT5G04680","AT5G16610",
        "AT5G53030","AT1G02950","AT1G05610","AT1G05770","AT1G10360",
        "AT1G13920"
    ],
    "Protoxylem": [
        "AT1G04540","AT1G04760","AT1G06020","AT1G07120","AT1G09610",
        "AT1G09890","AT1G10800","AT1G12450","AT1G14190","AT1G16520",
        "AT1G17950","AT1G19230","AT1G20720","AT1G21780","AT1G26570",
        "AT1G27920","AT1G29200","AT1G48280","AT1G51640"
    ],
    "Meri. Xylem": [
        "AT1G80100","AT4G37650","AT1G29950","AT4G02750","AT4G12620",
        "AT1G74150","AT5G08260","AT5G63140","AT2G25060","AT3G14730",
        "AT3G22790","AT1G29940","AT5G06210","AT3G61100","AT2G13540",
        "AT3G53190","AT5G37010","AT5G25475","AT3G54090","AT3G20430"
    ],
    "Atrichoblast": [
        "AT4G00730","AT5G66800","AT1G68560","AT1G31950","AT4G25250",
        "AT1G79840","AT5G19750","AT4G00360"
    ],
    "Trichoblast": [
        "AT1G03550","AT1G45545","AT1G48930","AT1G69240","AT2G46380",
        "AT3G51410","AT3G60630","AT4G04930","AT1G13950","AT2G03220",
        "AT4G07960","AT4G25160","AT5G07450","AT5G65160","AT4G09500",
        "AT3G46270","AT4G16920","AT1G07795"
    ],
}


def get_dataset_dir(data_path):
    return Path(data_path).expanduser().resolve() / "ArabidopsisSC"


def download_datasets(data_path):
    dataset_dir = get_dataset_dir(data_path)
    dataset_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "GSE123818_Root_single_cell_wt_datamatrix.csv.gz": (
            "https://www.ncbi.nlm.nih.gov/geo/download/"
            "?acc=GSE123818&format=file&file=GSE123818%5FRoot%5Fsingle%5Fcell%5Fwt%5Fdatamatrix%2Ecsv%2Egz"
        ),
        "GSE123818_Root_single_cell_shr_datamatrix.csv.gz": (
            "https://www.ncbi.nlm.nih.gov/geo/download/"
            "?acc=GSE123818&format=file&file=GSE123818%5FRoot%5Fsingle%5Fcell%5Fshr%5Fdatamatrix%2Ecsv%2Egz"
        ),
    }

    downloaded = {}
    for fname, url in files.items():
        local_path = pooch.retrieve(
            url=url,
            known_hash=None,
            fname=fname,
            path=dataset_dir,
            progressbar=True,
        )
        downloaded[fname] = Path(local_path)
    print("done")
    return downloaded


def get_ge_raw(data_path):
    ge_path = get_dataset_dir(data_path) / "GSE123818_Root_single_cell_wt_datamatrix.csv.gz"
    return pd.read_csv(ge_path, index_col=0)


def generate_clusters(
    ge_data,
    min_cells=2,
    min_genes=100,
    max_genes=12000,
    max_counts=120000,
    n_hvg=2000,
    n_neighbors=10,
    n_pcs=40,
    resolution=0.35,
):
    ge_data = ge_data.T

    adata = ad.AnnData(ge_data.values)
    adata.obs_names = ge_data.index.astype(str)
    adata.var_names = ge_data.columns.astype(str)

    adata.obs["n_genes_by_counts"] = (adata.X > 0).sum(axis=1)
    adata.obs["total_counts"] = adata.X.sum(axis=1)

    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes].copy()
    adata = adata[adata.obs["total_counts"] < max_counts].copy()

    adata.raw = adata.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=False)

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)

    sc.tl.pca(adata_hvg, svd_solver="arpack")
    sc.pp.neighbors(adata_hvg, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata_hvg, resolution=resolution)

    adata.obs["leiden"] = adata_hvg.obs["leiden"]

    return adata


def annotate_clusters_from_markers(adata, marker_dict, cluster_key="leiden"):
    cluster_to_label = {}
    cluster_scores = {}

    clusters = adata.obs[cluster_key].astype(str).unique()

    for cl in sorted(clusters):
        cells = adata[adata.obs[cluster_key].astype(str) == cl]

        scores = {}
        for label, genes in marker_dict.items():
            found = [g for g in genes if g in cells.var_names]

            if not found:
                scores[label] = 0.0
                continue

            mat = cells[:, found].X
            if hasattr(mat, "toarray"):
                mat = mat.toarray()

            scores[label] = mat.mean()

        best_label = max(scores, key=scores.get)
        cluster_to_label[cl] = best_label
        cluster_scores[cl] = scores

    return cluster_to_label, cluster_scores


def save_clustered_annotated_h5ad(data_path, marker_dict):
    ge_data = get_ge_raw(data_path)
    adata = generate_clusters(ge_data)
    print("clusters generated")
    cluster_labels, cluster_scores = annotate_clusters_from_markers(adata, marker_dict)
    adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(cluster_labels).astype("category")

    for cluster_id, label in cluster_labels.items():
        adata.uns[f"cluster_{cluster_id}_label"] = label
        adata.uns[f"cluster_{cluster_id}_scores"] = cluster_scores[cluster_id]

    out_path = get_dataset_dir(data_path) / "wt_clustered_annotated.h5ad"
    adata.write_h5ad(out_path, compression="gzip")
    return out_path


def get_cluster_expression(adata, clusters, mode="raw"):
    if isinstance(clusters, (str, int)):
        clusters = [clusters]

    clusters = [str(c) for c in clusters]
    mask = adata.obs["leiden"].astype(str).isin(clusters)
    adata_cluster = adata[mask].copy()

    if adata_cluster.n_obs == 0:
        raise ValueError(f"No cells found in clusters {clusters!r}")

    if mode == "lognorm":
        mat = adata_cluster.X
        var_names = adata_cluster.var_names
    elif mode == "raw":
        if adata_cluster.raw is None:
            raise ValueError("adata.raw is None; raw counts not available")
        mat = adata_cluster.raw.X
        var_names = adata_cluster.raw.var_names
    else:
        raise ValueError("mode must be 'lognorm' or 'raw'")

    if hasattr(mat, "toarray"):
        mat = mat.toarray()

    return pd.DataFrame(mat, index=adata_cluster.obs_names, columns=var_names)


def normalize_log_ge(ge_data, target_sum=1e4):
    cell_sums = ge_data.sum(axis=1)
    cell_sums[cell_sums == 0] = 1.0
    norm_counts = ge_data.div(cell_sums, axis=0) * float(target_sum)
    return np.log1p(norm_counts + 1.0)


def find_clusters_for_cell_type(adata, cell_type):
    mask = adata.obs["cell_type"].astype(str) == str(cell_type)
    clusters = adata.obs.loc[mask, "leiden"].astype(str).unique().tolist()

    if not clusters:
        raise ValueError(f"No clusters found for cell type {cell_type!r}")

    return sorted(clusters)


def save_gene_expression_for_cell_type(data_path, cell_type, output_name=None):
    h5ad_path = get_dataset_dir(data_path) / "wt_clustered_annotated.h5ad"
    adata = sc.read_h5ad(h5ad_path)

    clusters = find_clusters_for_cell_type(adata, cell_type)
    raw_expr = get_cluster_expression(adata, clusters, mode="raw")
    norm_expr = normalize_log_ge(raw_expr)

    if output_name is None:
        safe_name = str(cell_type).replace(" ", "_").replace(".", "")
        output_name = f"{safe_name}_lognorm.csv.gz"

    out_path = get_dataset_dir(data_path) / output_name
    norm_expr.to_csv(out_path, index=True, compression="gzip")
    return out_path


def build_parser():
    parser = argparse.ArgumentParser(
        description="Arabidopsis single-cell dataset setup",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--env", required=True, metavar="PATH", help="Path to the .env file.")
    parser.add_argument("--step", required=True, choices=["download", "cluster", "expression"])
    parser.add_argument("--cell-type", help="Cell type label for step=expression, e.g. Endodermis")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    env = ut.get_env(args.env)
    data_path = Path(env["DATA_PATH"]).expanduser().resolve()

    if args.step == "download":
        download_datasets(data_path)

    elif args.step == "cluster":
        save_clustered_annotated_h5ad(data_path, MARKER_GENES)

    elif args.step == "expression":
        if not args.cell_type:
            raise ValueError("--cell-type is required for step=expression")
        save_gene_expression_for_cell_type(data_path, args.cell_type)


if __name__ == "__main__":
    main()