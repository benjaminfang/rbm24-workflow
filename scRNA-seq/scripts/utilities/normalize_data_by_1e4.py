"""
1) read in ambient corrected data.
2) remove gene which less than three cells have it(gene_n < 3 droped).
3) normalize the count data to 1e4 by cells.
"""

import numpy as np
import argparse
from io_10x import read_10x_as_dense_mtx, write_dense_mtx_10x
import os


def filter_gene_have_cell(mtx, N=3):
    dt_out = []
    for row in mtx:
        if np.sum(row > 0) >= 3:
            dt_out.append(True)
        else:
            dt_out.append(False)
    return dt_out
            

def norm_mtx(mtx, N=1e4):
    mtx_sum = np.sum(mtx, axis=0)
    mtx = mtx * N / mtx_sum
    return mtx


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("count_dire", help="directory where the matrix data is stored")
    parser.add_argument("--out_dire", help="output directory", default="out")
    args = parser.parse_args()
    return args.count_dire, args.out_dire


def main_norm_data():
    count_dire, out_dire = getargs()
    bar_f = os.path.join(count_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(count_dire, "features.tsv.gz")
    mtx_f = os.path.join(count_dire, "matrix.mtx.gz")
    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")

    keep = filter_gene_have_cell(mtx)
    mtx = mtx[keep]
    features = list(np.asarray(features)[keep])
    mtx = norm_mtx(mtx)
    if not os.path.exists(out_dire):
        os.mkdir(out_dire)
    write_dense_mtx_10x(barcodes, features, mtx, out_dire)


if __name__ == "__main__":
    main_norm_data()