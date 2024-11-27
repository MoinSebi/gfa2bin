#!/usr/bin/env python3
import logging

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib

def getLog(file: str, n: int) -> (list, int):
    df = pd.read_csv(file, sep = "\t")
    last_col = df.columns[-1]
    df["log"] = - np.log10(df[last_col])
    p = sorted(df["log"])
    o = len(p)
    p = p[:int(o*0.9)][::n] + p[o:len(p)]
    return p, o


def get_uniform_dist(o: int, n: int) -> list:
    uni = sorted(-np.log10(np.linspace(1/o,1,o)))
    uni = uni[:int(o*0.9)][::n] + uni[o:len(uni)]
    return uni


def qqplot(p: list, uni: list, filename: str):
    plt.figure(figsize=(7,7))
    maxx = max(max(p), max(uni))
    plt.scatter(uni, p, s = 3, color = "blue")
    x1 = [-2, maxx]
    x2 = [-2, 20]
    plt.plot(x1, x2, marker = 'o', c = "black")
    plt.xlim([-0.5,maxx + 0.5])
    plt.ylim([-0.5,maxx + 0.5])
    plt.xlabel("Expected -log$_{10}$ ($\it{P}$ value)")
    plt.ylabel("Observed -log$_{10}$ ($\it{P}$ value)")
    plt.legend()
    plt.savefig(filename + ".qq.thresh.pdf", dpi = 800)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="QQ plotter from GEMMA input.")
    parser.add_argument('-i', '--input', type=str, help='Path to the GEMMA file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file name (pdf)', required=True)
    parser.add_argument('-s', '--skip', type=int, help='Skip every N entry', required=False)
    args = parser.parse_args()

    logging.info(f"Reading input file: {args.input}")
    p, o = getLog(args.input, args.skip)

    print(p)

    logging.info(f"Calculating uniform distribution")
    uni = get_uniform_dist(o, args.skip)

    logging.info(f"Plotting QQ plot")
    qqplot(p, uni, args.output)
