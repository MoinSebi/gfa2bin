#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: moinSebi

"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')



def getLog(file): 
    df = pd.read_csv(file, sep = "\t")
    df["log"] = - np.log10(df["p_lrt"])
    p = sorted(df["log"])
    #uni = sorted(np.random.normal(size = len(p)))
    uni = sorted(-np.log10(np.linspace(1/len(p),1,len(p))))

    return (p,uni)


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file (FAM)", required=True)

    parser.add_argument("-o", "--output", help="Output file (FAM)", required=True)

    args = parser.parse_args()

    print("Reading")
    h1 = getLog(args.input)
    print("Plotting")
    plt.figure(figsize=(7,7))
    print("Plotting")

    plt.scatter(h1[1], h1[0], s = 1)
    print("Plotting")
    x1 = [0, 10]
    x2 = [0, 10]
    plt.plot(x1, x2, marker = 'o', c = "black")
    print("Plotting")


    a = max(h1[0])
    plt.xlim([0,a+1])
    plt.ylim([0,a+1])
    plt.xlabel("Expected (-logP)")
    plt.ylabel("Observed (-logP)")
    print("Plotting")

    plt.savefig(args.output, dpi = 400)
    plt.close()