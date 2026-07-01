import pandas as pd
from tfitpy.datasets import load
from tfitpy.utils import generate_tf_pairs
from tfitpy import validate
import src.results.util as ut
import src.pipeline.util as put
from joblib import Parallel, delayed, dump
from pathlib import Path
import os
import json

def add_cluster_indices(cluster_df, data_path, ge_data, organism="human"):
    """
    cluster_df must have at least a 'target' column.
    Computes and returns score columns for cluster_df and
    a dict mapping each uid/cluster_uid to its evidence dict.
    """

    # load the datasets
    target_list = cluster_df["target"].tolist()
    cache = load(
        data_path,
        gene_expression_data=ge_data,
        targets=target_list,
        organism=organism
    )

    evidence_objects = {}
    score_rows = []

    def get_row_uid(row):
        if "cluster_uid" in row.index:
            return row["cluster_uid"]
        elif "uid" in row.index:
            return row["uid"]
        else:
            raise KeyError("Neither 'cluster_uid' nor 'uid' present in row")

    use_pairwise_cache = True
    go_ea = False
    if organism == "arabidopsis":
        use_pairwise_cache = False
        go_ea = True


    for _, row in cluster_df.iterrows():
        sources = row["sources"].split(";")
        source_pairs = generate_tf_pairs(sources)
        target = row["target"]
        uid = get_row_uid(row)
        # print(uid)
        scores, evidence = validate(
            sources,
            target,
            cache,
            use_pairwise_cache=use_pairwise_cache,
            data_path=data_path,
            go_ea=go_ea,
            pairs=source_pairs,
            organism=organism
        )

        

        # store evidence dict per uid
        evidence_objects[uid] = evidence

        # store scores dict to become new columns
        score_rows.append(scores)

    # turn list of score dicts into a DataFrame and join
    new_cols = pd.DataFrame(score_rows, index=cluster_df.index)
    df = cluster_df.join(new_cols)

    return df, evidence_objects


def add_cluster_indices_main(
    data,
    data_path,
    ge_data,
    organism="human",
    parallel=True,
    n_jobs=4,
    batch_size=1000,
):
    """
    Main wrapper for index generation.
    Returns:
        combined_df, combined_evidence
    """

    if not parallel:
        return add_cluster_indices(data, data_path, ge_data, organism=organism)

    batches = [
        data.iloc[i:i + batch_size].copy()
        for i in range(0, len(data), batch_size)
    ]
    print(len(batches))
    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(add_cluster_indices)(
            batch,
            data_path,
            ge_data,
            organism
        )
        for batch in batches
    )

    df_parts = [result[0] for result in results]
    combined_df = pd.concat(df_parts, axis=0)

    evidence_parts = [result[1] for result in results]
    combined_evidence = {}
    for evidence_dict in evidence_parts:
        combined_evidence.update(evidence_dict)

    return combined_df, combined_evidence


def read_results_file(input, env, options, args):
    """"""
    out_path, temp_path = ut.get_exp_path(input, env)
    bio_path = Path(os.path.expandvars(env["DATA_PATH"]))
    if input["type"] == "run_coregnet":
        # print()
        f_path = out_path / f"{input['result_file_name']}.csv"
        # print(f_path)
        data = pd.read_csv(f_path)
        return data
    elif input["type"] == "run_coregtor":
        # print()
        cname = options.get("cluster_name",None)
        if cname is None:
            raise ValueError("no cluster_name in the options")
        f_path = temp_path / "clusters"/ f"{cname}.csv"
        data = pd.read_csv(f_path)
        return data
    else:
        raise ValueError("invalid exp type")


def save_results_file(input, env, options, args,result_df,evidence):
    """
    """
    out_path, temp_path = ut.get_exp_path(input, env)
    bio_path = Path(os.path.expandvars(env["DATA_PATH"]))
    if input["type"] == "run_coregnet":
        f_path = out_path / f"{input['result_file_name']}_indices.csv"
        result_df.to_csv(f_path,index=False)

        f_evid = out_path / f"{input['result_file_name']}_evidence.pkl"
        dump(evidence, f_evid, compress=3)
    elif input["type"] == "run_coregtor":
        cname = options.get("cluster_name",None)
        f_path = out_path / f"{cname}_indices.csv"
        result_df.to_csv(f_path,index=False)

        f_evid = out_path / f"{cname}_evidence.pkl"
        dump(evidence, f_evid, compress=3)
    else:
        raise ValueError("invalid exp type")




def compute_validation_indices(input, env, options, args):
    """
    """
    out_path, temp_path = ut.get_exp_path(input, env)
    bio_path = Path(os.path.expandvars(env["DATA_PATH"]))
    rerun = args.rerun
    batch = args.batch
    njobs = args.njobs if args.njobs is not None else 4
    if njobs == -1:
        njobs = 4
    organism = input.get('organism',"human")
    #   print("hi")
    # read the gene expression data
    ge_data = put.read_dataset(input["dataset"], env)

    # read the file with clusters
    cluster_df = read_results_file(input, env, options, args)
    #   print(cluster_df)

    results, evidence = add_cluster_indices_main(cluster_df, data_path=bio_path,ge_data=ge_data,organism=organism,parallel=True,n_jobs=njobs,batch_size=batch)

    # print(results)
    save_results_file(input, env, options, args,results, evidence)
    # first create and save indices file
