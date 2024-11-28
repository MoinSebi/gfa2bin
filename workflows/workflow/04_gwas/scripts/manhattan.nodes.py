#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s    [%(levelname)s] - %(filename)s: %(message)s',
    datefmt='%d/%m/%Y %H:%M:%S',  # 24-hour format
    handlers=[logging.StreamHandler(stream=sys.stderr)]
)

def read_assoc(filename: str) -> pd.DataFrame:
    """
    Parse associations files from GEMMA
    :param filename: GEMMA output file
    :return: Data in pandas Dataframe
    """
    df = pd.read_csv(filename, sep = "\t")
    last_col = df.columns[-1]
    df["log"] = - np.log10(df[last_col])

    return df


def plotting_manhattan(df: pd.DataFrame, filename: str):
    df_len = len(df)
    df = df.loc[df["log"] > 2]

    plt.figure(figsize=(10,6))
    plt.ylim(0, max(df["log"])+1)

    plt.scatter([int(x) for x in df["ps"]],df["log"], s = 3, color = "blue")
    if df_len > 0:
        plt.axhline(y=-np.log10(0.05/df_len), color='r', linestyle='-')
    plt.xlabel("Node id")
    plt.ylabel("-log(p)")
    plt.savefig(filename, dpi=800)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manhattan plotter from GEMMA input (nodes).")
    parser.add_argument('-i', '--input', type=str, help='Path to the input file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to the output file (pdf[default]/png)', required=True)
    args = parser.parse_args()

    logging.info(f"Reading input file: {args.input}")
    df = read_assoc(args.input)

    logging.info(f"Plotting manhattan plot")

    # Check file name extension
    if not args.output.endswith(".pdf") and not args.output.endswith(".png"):
        args.output += ".pdf"
    plotting_manhattan(df, args.output)