#!/usr/bin/env python2.7
import sys
import cPickle as pik


def main(args):
    allele = pik.load(open(args.input, 'rb'))
    for k in allele:
        if allele[k][0] == 0:
            # wild type
            wt = None


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Calculate fitness scores.")
    parser.add_argument("input", nargs="?")
    sys.exit(parser.parse_args())
