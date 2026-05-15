"""
This is to setup ppi network analysis 
"""
import src.results.util as ut
import argparse
from pathlib import Path

from tfitpy.network_analysis_ppi.generate_null import generate_null_networks


def generate_files(env,args):
    """"""
    data_path = Path(env.get("DATA_PATH"))
    ppis = ["hippie","stringdb","biogrid"]
    options = {
        "njobs":  args.njobs if args.njobs is not None else -1 ,
        "nmodels":  200,
        "rerun": args.rerun if args.rerun is not None else False,
    }
    for p in ppis:
        opt = {**options}
        opt["ppi_key"] = p
        generate_null_networks(data_path,opt)



def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate setup files for ppi network analysis",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--env",   required=True,
                        metavar="PATH", help="Path to the .env file.")
    

    # Optional arguments — add more as needed
    parser.add_argument("--njobs", type=int, default=-1,
                        help="Number of parallel jobs. Defaults to -1 (all cores).")
    parser.add_argument("--rerun", action="store_true", default=False,
                        help="Force rerun even if output already exists.")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    env = ut.get_env(args.env)
    generate_files(env,args)


if __name__ == "__main__":
    main()