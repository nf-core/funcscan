#!/usr/bin/env python


"""Provide a command line tool to merge AMP prediction tables."""


import argparse
import logging
from pathlib import Path

import pandas as pd


logger = logging.getLogger()


def join_probabilities(tables):
    """
    Join all tables together.

    Tables are joined using an outer join which ensures that all sequence predictions
    are in the resulting table but this may introduce missing values.

    """
    result = tables[0].join(tables[1:], how="outer", sort=True)
    return result.reset_index()


def parse_macrel(filename):
    """Parse tab-delimited macrel AMP predictions and return a homogenous table."""
    df = pd.read_table(filename, sep="\t", skiprows=1)
    df.rename(columns={"Access": "sequence_id", "AMP_probability": "macrel"})
    return df[["sequence_id", "macrel"]].set_index("sequence_id")


def parse_ampir(filename):
    """Parse tab-delimited ampir AMP predictions and return a homogenous table."""
    # TODO: Define column names if the table doesn't have any.
    df = pd.read_table(filename, sep="\t", header=None, names=[])
    # df.rename(columns={"?": "sequence_id", "?": "ampir"})
    return df[["sequence_id", "ampir"]].set_index("sequence_id")


def parse_arguments(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge different AMP prediction tables into one.",
    )
    parser.add_argument(
        "--macrel",
        metavar="PATH",
        type=Path,
        help="Input file in tab-delimited format containing macrel predictions.",
    )
    parser.add_argument(
        "--ampir",
        metavar="PATH",
        type=Path,
        help="Input file in tab-delimited format containing ampir predictions.",
    )
    parser.add_argument(
        "--output",
        metavar="PATH",
        type=Path,
        default=Path("predictions.tsv.gz"),
        help="Output file in tab-delimited format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """"""
    args = parse_arguments(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    # Ensure that directories for the output path exist.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tables = []
    if args.macrel:
        assert args.macrel.is_file(), f"File not found: {args.macrel}"
        logger.info("Parse macrel predictions.")
        macrel_df = parse_macrel(args.macrel)
        tables.append(macrel_df)
    if args.ampir:
        assert args.ampir.is_file(), f"File not found: {args.ampir}"
        logger.info("Parse ampir predictions.")
        ampir_df = parse_ampir(args.ampir)
        tables.append(ampir_df)
    # TODO: add other prediction tools
    result = join_probabilities(tables)
    result.to_csv(args.output, header=True, index=True, sep="\t")


if __name__ == "__main__":
    main()
