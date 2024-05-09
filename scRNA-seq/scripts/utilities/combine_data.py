import sys
import gzip
from io_10x import read_10x_as_dense_mtx, write_dense_mtx_10x
import os
import numpy as np


def load_data(dd):

    bar_f = os.path.join(dd, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    fea_f = os.path.join(dd, "filtered_feature_bc_matrix", "features.tsv.gz")
    mtx_f = os.path.join(dd, "filtered_feature_bc_matrix", "matrix.mtx.gz")
    bar, fea, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "int32")
    return bar, fea, mtx


def filter_bar(bar, mtx):
    mtx_sum = np.sum(mtx, axis=0)
    f_bool = mtx_sum > 1000
    mtx = mtx.T[f_bool]
    mtx = mtx.T
    bar = np.asarray(bar)[f_bool]

    return bar, mtx


def check_consistence(fea1, fea2):
    fea1 = [e[0] for e in fea1]
    fea2 = [e[0] for e in fea2]
    for e in zip(fea1, fea2):
        if e[0] != e[1]:
            print("error")
            break


if __name__ == "__main__":
    bar1, fea1, mtx1 = load_data(sys.argv[1])
    bar2, fea2, mtx2 = load_data(sys.argv[2])
    germ_cell = [l.rstrip().split()[0] for l in open(sys.argv[3])]
    germ_marker = [l.rstrip() for l in open(sys.argv[4])]
    print(mtx1.shape)
    print(mtx2.shape)
    bar1, mtx1 = filter_bar(bar1, mtx1)
    bar2, mtx2 = filter_bar(bar2, mtx2)
    print(mtx1.shape)
    print(mtx2.shape)
    check_consistence(fea1, fea2)

    bar = np.hstack([bar1, bar2])
    mtx = np.hstack([mtx1, mtx2])

    germ_marker_idx = []
    germ_idx = []
    none_idx = []

    marker_align = []
    for i in range(len(fea1)):
        if fea1[i][0] in germ_marker:
            germ_marker_idx.append(i)
            marker_align.append(fea1[i][0])
    for i, bb in enumerate(bar):
        if bb in germ_cell:
            germ_idx.append(i)
        else:
            none_idx.append(i)
    
    mtx = mtx[germ_marker_idx]
    germ_mtx = mtx.T[germ_idx]
    none_mtx = mtx.T[none_idx]

    print(marker_align)
    print(np.sum(germ_mtx, axis=0))
    print(np.sum(none_mtx, axis=0))
