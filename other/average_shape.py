from __future__ import division
import pandas as pd
import itertools
import argparse
import numpy as np
import os.path as op

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Name of the directory to be processed', type=str)
    args = parser.parse_args()

    # Get data and store in dataframe
    shapeDf = pd.read_csv(args.file,sep='\t',header=None,na_values='NA')

    means = shapeDf.mean(axis=0)

    # Print a header row
    print op.basename(args.file), '\t'.join(map(str,means))

    # for i,j in enumerate(means):
    #     print i, '\t', j
