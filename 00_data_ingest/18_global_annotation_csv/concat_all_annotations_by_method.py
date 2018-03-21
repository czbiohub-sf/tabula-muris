#!/usr/bin/env python3.6

import argparse
import csv
import glob
import os.path

import pandas as pd
import numpy as np

METHODS = ("facs", "droplet")
INPUT_ANNOTATIONS_DIR = "../03_tissue_annotation_csv"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", type=str, choices=METHODS, required=True)
    return parser.parse_args()


def get_header(csv_filename):
    with open(csv_filename) as csv_file:
        header = next(csv.reader(csv_file))
    return header


subtissues = "Aorta", "Diaphragm"

def main():
    args = parse_args()
    method = args.method

    all_annotations_csv_name = f"annotations_{method}.csv"
    if os.path.isfile(all_annotations_csv_name):
        raise RuntimeError(f"{all_annotations_csv_name} already exists. Please remove before running.")
    input_annotations = glob.iglob(f"{INPUT_ANNOTATIONS_DIR}/*_{method.lower()}_annotation.csv")

    dfs = []
    for annotation_file in input_annotations:
        df = pd.read_csv(annotation_file)

        # These subtissues were separated out for analysis and should just be
        # called the "tissue"
        if (method == 'facs') and (df['subtissue'][0] in subtissues):
            df['tissue'] = df['subtissue']
            df['subtissue'] = np.nan

        dfs.append(df)
    all_annotations = pd.concat(dfs)

    # index=False means don't write row numbers in the file
    all_annotations.to_csv(all_annotations_csv_name, index=False)


if __name__ == "__main__":
    main()
