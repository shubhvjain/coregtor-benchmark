import argparse
import json
import sys
from pathlib import Path
from dotenv import dotenv_values
import src.results.util as ut
from joblib import Parallel, delayed, load
import pandas as pd

from src.analysis.gtex import generate_gtex_stats

def hello(input, env, options, args):
    target = input.get("target", "world")
    print(f"Hello, {target}!")
    print(f"  env keys : {list(env.keys())}")
    print(f"  verbose  : {args.verbose}")



METHOD_REGISTRY = {
    "hello": hello,
    "generate_gtex_stats":generate_gtex_stats
}

# CLI entry point


def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate analysis",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--env",   required=True,
                        metavar="PATH", help="Path to the .env file.")
    parser.add_argument("--input", required=True,
                        metavar="PATH", help="Path to the JSON input file.")
    parser.add_argument(
        "--id",
        required=True,
        metavar="ID",
        help="ID of the result object in input['results'] to run.",
    )

    # Optional arguments — add more as needed
    parser.add_argument("--njobs", type=int, default=-1,
                        help="Number of parallel jobs. Defaults to -1 (all cores).")
    parser.add_argument("--rerun", action="store_true", default=False,
                        help="Force rerun even if output already exists.")
    parser.add_argument("--verbose", action="store_true",
                        default=False, help="Enable verbose output.")
    parser.add_argument("--dry-run", action="store_true",
                        default=False, help="Load files but skip execution.")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    env = ut.get_env(args.env)
    input_data = ut.get_input(args.input)

    options = ut.get_result_by_id(input_data, args.id)

    result_type = options.get("type")
    if result_type not in METHOD_REGISTRY:
        print(
            f"[error] unknown type '{result_type}'. Registered types: {', '.join(METHOD_REGISTRY.keys())}", file=sys.stderr)
        sys.exit(1)

    fn = METHOD_REGISTRY[result_type]

    if args.dry_run:
        print(
            f"[dry-run] would call '{result_type}' for id '{args.id}' with {len(env)} env key(s).")
        return

    fn(input_data, env, options, args)


if __name__ == "__main__":
    main()
