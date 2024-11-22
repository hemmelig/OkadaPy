"""
A command-line interface to the `OkadaPy` package.

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import argparse
from pathlib import Path

from okada import read_model
from okada import evaluate


def main(args=None):
    """Entry point for the `okada` command-line utility."""

    parser = argparse.ArgumentParser()
    sub_parser = parser.add_subparsers(
        title="commands",
        dest="command",
        required=True,
        help="Specify the mode to run in.",
    )

    evaluate_parser = sub_parser.add_parser(
        "evaluate", help="Evaluate the analytical solution for a given model."
    )
    evaluate_parser.add_argument(
        "-i",
        "--input",
        help="Path to an okada-compatible model file.",
        required=True,
    )
    evaluate_parser.add_argument(
        "-m",
        "--mode",
        help="Specify the mode in which to evaluate the model.",
        choices=["displacement", "stress", "strain"],
        required=True,
    )
    evaluate_parser.add_argument(
        "-d",
        "--depth",
        help="Depth at which to compute the analytical solution.",
        required=True,
        type=float,
    )
    evaluate_parser.add_argument(
        "-t",
        "--threads",
        help="Specify the number of threads to use for computation.",
        required=False,
        type=int,
        default=1,
    )
    evaluate_parser.add_argument(
        "-o",
        "--output-dir",
        help="Specify an output directory.",
        required=False,
        default=Path.cwd(),
    )
    args = parser.parse_args(args)

    print("=" * 79)
    print("\tOkadaPy - analytical solutions for deformation, strain, and stress")
    print("=" * 79)

    model = read_model(args.input)
    solution = evaluate(model, args.depth, args.mode, args.threads)

    # solution.write()
    # print(solution.strain)

    print("=" * 79)
