#!/usr/bin/env python2.7
import sys
import cPickle as pik
import math
from collections import defaultdict
import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt


def main(args):
    data = defaultdict(lambda: np.full(len(args.input) + 1, np.nan))
    for i, fname in enumerate(args.input):
        allele = pik.load(open(fname, 'rb'))
        wtcnt = 0.0
        total = 0.0
        for k in allele:
        # allele[k] == [position, codon, residue, count]
            total += allele[k][3]
            if allele[k][0] == 0:
                # wild type
                wtcnt += allele[k][3]

        for k in allele:
            fwt = wtcnt / total
            fmut = allele[k][3] / total
            if fmut > 0:
                fitness = math.log(fmut / fwt, 2)
            else:
                fitness = 0
            data[k][i] = fitness

    for k in data:
        # data[k] # numpy array of 3 values
        y = data[k][0], data[k][1], data[k][2]
        x = 0, 1.87, 3.82
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        data[k][3] = slope
    slopes = np.array([data[k][3] for k in data])
    plt.hist(slopes, 25)
    plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Calculate fitness scores.")
    parser.add_argument("input", nargs=3)
    sys.exit(main(parser.parse_args()))
