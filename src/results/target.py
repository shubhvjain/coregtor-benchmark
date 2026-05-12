import src.results.util as ut
from joblib import load, Parallel, delayed

from src.results.explore import _get_gene_stats, _get_dist_stats
import json


def get_target_info(target, input, data):
    """
    generate a dict with all info related to a target from its pkl run file
    """
    gf_status = _get_gene_stats(data["results"]["gene_frequency"])
    results = {
        "target": target,
        "gf": gf_status,   
    }
    dmeasures = input["context_comparison"]["methods"]
    for dms in dmeasures:
        dmatrix = data["results"]["distance"][dms]
        stat = _get_dist_stats(dmatrix)
        results[f"{dms}"] = stat
    # print(data)
    results["stats"] = {
        "n_cores_used": data["stats"]["n_cores"],
        "n_paths":data["stats"]["quality"]["n_paths"],
        "n_unique_roots":data["stats"]["quality"]["n_unique_roots"],
    
        "model_train_time": round(data["stats"]["timing"]["model_train"],5),
        "tree_path_time" : round(data["stats"]["timing"]["context_create"],5),
        "paths_extract_time": round(data["stats"]["timing"]["paths_extract"],5),
        "total_time": round(data["stats"]["timing"]["model_train"]+ data["stats"]["timing"]["context_create"]+ data["stats"]["timing"]["paths_extract"] ,5) ,
    }
    
    return results


def target_info(input, env, options, args):
    """
    Create stats about targets and store them in a json.
    """
    out_path, temp_path = ut.get_exp_path(input, env)
    rerun = args.rerun
    n_jobs = args.njobs if args.njobs is not None else 4
    out_file = out_path / "target_stats.json"

    if out_file.exists() and not rerun:
        print("File already exists")
        return

    all_targets = ut.get_exp_target_list(out_path)

    def process_target(target):
        try:
            data = ut.get_temp_file(temp_path, target)
        except FileNotFoundError:
            return {"target": target, "success": False}

        result = get_target_info(target, input, data)
        result["success"] = True
        return result

    results = Parallel(n_jobs=n_jobs)(
        delayed(process_target)(target) for target in all_targets
    )

    with open(out_file, "w") as f:
        json.dump(results, f, indent=1)  # results is already a list

    return results


def get_target_crm(input, env, options, args):
    """
    Returns the set of predicted co-regulatory modules for the target in an experiment as a dataframe 
    """
    result_type = options.get("result_name", "default")
    out_path, temp_path = ut.get_exp_path(input, env)


"""
To interactively explore an experiment
"""


class Explore:
    def __init__(self, env_path, exp_file_path):
        self.env = ut.get_env(env_path)
        self.input = ut.get_input(exp_file_path)
        output_path, temp_path = ut.get_exp_path(self.input, self.env)
        self.output_path = output_path
        self.temp_path = temp_path
        self.default_args = {"interactive": True}
        self.target_cache = {}
        print("done")

    def _target_cache(self, target):
        if target in self.target_cache:
            return self.target_cache[target]
        else:
            tfile_path = self.temp_path / "results" / f"{target}.pkl"
            tfile = load(tfile_path)
            self.target_cache[target] = tfile
            return tfile

    def target_gf(self, target):
        """"""
        data = self._target_cache(target)
        return data["results"]["gene_frequency"]

    def target_gf_stats(self, target):
        """"""
        gf = self.target_gf(target)
        return _get_gene_stats(gf)

    def target_crm(self, target, result_name="default"):
        """"""
        options = {"result_type": result_name, "target": target}
        args = {"interactive": True}
        return get_target_crm(self.input, self.env, options, self.default_args)
