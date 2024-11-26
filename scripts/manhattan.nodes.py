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
    df["log"] = - np.log10(df["p_lrt"])

    return df


def plotting_manhattan(df: pd.DataFrame, filename: str):
    df = df.loc[df["log"] > 2]

    plt.figure(figsize=(10,6))


    plt.scatter([int(x) for x in df["pos"]],df["log"], s = 3, color = "blue")
    plt.axhline(y=-np.log10(0.05/(len(df))), color='r', linestyle='-')
    plt.xlabel("Node id")
    plt.ylabel("-log(p)")
    plt.savefig(filename + ".manhatten.pdf")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manhattan plotter from GEMMA input (nodes).")
    parser.add_argument('-i', '--input', type=str, help='Path to the input file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to the output file', required=True)
    args = parser.parse_args()

    logging.info(f"Reading input file: {args.input}")
    df = read_assoc(args.input)

    logging.info(f"Plotting manhattan plot")
    plotting_manhattan(df, args.output)