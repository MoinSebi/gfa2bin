#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import logging
import sys

from matplotlib.cm import ScalarMappable, get_cmap
from matplotlib.colors import Normalize

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s    [%(levelname)s] - %(filename)s: %(message)s',
    datefmt='%d/%m/%Y %H:%M:%S',  # 24-hour format
    handlers=[logging.StreamHandler(stream=sys.stderr)]
)


def read_assoc(f: str) -> pd.DataFrame:
    """
    Read the gemme assoc file

    :param f:Filename
    :return:
    """
    df = pd.read_csv(f, sep = "\t")
    df["log"] = - np.log10(df["p_lrt"])
    df["node"] = df["ps"]
    return df

def read_distance(filename: str) -> pd.DataFrame:
    """
    Read gfa2bin nearest distance table

    :param filename: File name of gretl table
    :return: Data in pandas Dataframe
    """
    df_dis = pd.read_csv(filename, sep = "\t")
    df_dis = df_dis.sort_values("node")
    df_dis["node"] = df_dis["node"].astype(int)
    return df_dis

def merge_df(dfpval, dfdis) -> pd.DataFrame:
    """
    Merge two pandas Dataframe by nodes

    :param dfpval: Dataframe of the plvalues (GEMMA)
    :param dfdis: Dataframe of the gfa2bin distance
    :return: Merged dataframe
    """
    result = pd.merge(dfdis, dfpval, on="node", how="left")
    result2 = result.sort_values("node")
    result3 = result2.loc[True != result2["los"].isna()]
    return result3

def filter_df(df) -> pd.DataFrame:
    """
    Filter the dataframe by pval

    :param df: Merged pandas Dataframe
    :return: Returned filtered Dataframe
    """
    result = df.loc[df["log"] > 2]
    result2 = result.sort_values("path")
    return result2


def plot(result4: pd.DataFrame, filename: str):
    """
    Plot the merged Dataframe as a manhattan plot. Color by distance, position by reference location
    X label by reference entry

    :param result4: Merged pandas Dataframe with all the information
    :param filename: Output plot (pdf)
    :return:
    """
    rollingm = 0
    cmap = get_cmap("plasma_r")  # Use the "oranges" colormap
    norm = Normalize(vmin=-1, vmax=1000)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    paths = list(set(result4["path"]))
    ticks = []#
    tick_labels = paths
    plt.figure(figsize = (10,5))
    for x in paths:
        dfbt1 = result4.loc[result4["path"] == x]
        dfbt1 = dfbt1.loc[dfbt1["distance"] < 0]
        # Map "dist" values to colors using the "oranges" colormap
        colors1 = sm.to_rgba(dfbt1["distance"])

        # Plot with a single marker style (e.g., 'o')
        plt.scatter(dfbt1["position"] + rollingm, dfbt1["log"], s=12, marker="o", c=colors1, alpha = 0.5)
        ticks.append(int(np.mean([rollingm, rollingm + max(dfbt1["position"])])))
        dfbt1 = result4.loc[result4["path"] == x]
        dfbt1 = dfbt1.loc[dfbt1["distance"] > 0]
        # Map "dist" values to colors using the "oranges" colormap
        colors1 = sm.to_rgba(dfbt1["distance"])

        # Plot with a single marker style (e.g., 'o')
        plt.scatter(dfbt1["position"] + rollingm, dfbt1["log"], s=12, marker="o", c=colors1, alpha = 0.5)
        print(dfbt1["position"])
        #ticks.append(int(np.mean([ro'llingm, rollingm + max(dfbt1["position"])])))


        rollingm += max(dfbt1["position"]) + 10000
        print(rollingm)

    plt.ylabel("-log$_{10}$ ($\it{P}$ value)")
    plt.xlabel("Reference name")
    plt.xticks(ticks=ticks, labels=tick_labels)
    plt.savefig(filename + ".manhattan.pdf")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manhattan plotter from GEMMA input (nodes).")
    parser.add_argument('-i', '--input', type=str, help='Path to the input file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to the output file', required=True)
    parser.add_argument('-d', '--distance', type=str, help="Distance file", required=True)
    args = parser.parse_args()

    logging.info(f"Reading input file: {args.input}")
    dfpval = read_assoc(args.input)


    logging.info(f"Read gfa2bin distance file: {args.distance}")
    dfdistance = read_distance(args.distance)

    logging.info("Merge the two Dataframes")
    df = merge_df(dfpval, dfdistance)

    logging.info("Filter Dataframe")
    dffilter = filter_df(df)

    logging.info("Manhattan plot")
    plot(dffilter, args.output)