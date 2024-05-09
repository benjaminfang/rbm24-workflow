import sys
import gzip
from io_10x import read_10x_as_dense_mtx, write_dense_mtx_10x
import os
import numpy as np

def load_data(dires):
    dt_out = []
    for dd in dires:
        print(dd)
        bar_f = os.path.join(dd, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
        fea_f = os.path.join(dd, "filtered_feature_bc_matrix", "features.tsv.gz")
        mtx_f = os.path.join(dd, "filtered_feature_bc_matrix", "matrix.mtx.gz")
        bar, fea, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "int32")
        dt_out.append([bar, fea, mtx])
    return dt_out


def pick_barcode(bars, barcodes, mtx):
    pick_idx = []
    for idx in range(len(barcodes)):
        if barcodes[idx] in bars:
            pick_idx.append(idx)
    mtx = mtx.T[pick_idx]
    return mtx.T


if __name__ == "__main__":
    barcodes = [l.decode().rstrip() for l in gzip.open(sys.argv[1])]
    
    sid_group = [[], [], [], []]
    for bar in barcodes:
        sid = bar.split("_")[0]
        if sid == "s95":
            sid_group[0].append(bar)
        elif sid == "s96":
            sid_group[1].append(bar)
        elif sid == "s97":
            sid_group[2].append(bar)
        elif sid == "s98":
            sid_group[3].append(bar)
        else:
            print("error")
            exit()

    out_dire = sys.argv[2]
    if not os.path.exists(out_dire):
        os.mkdir(out_dire)

    data = load_data(sys.argv[3:])

    mtx_clo = []
    for i in range(4):
        mtx_1 = pick_barcode(sid_group[i], data[i][0], data[i][2])
        print(mtx_1.shape)
        mtx_clo.append(mtx_1)
    mtx_all = np.hstack(mtx_clo)
    sid_group_all = []
    for gg in sid_group:
        sid_group_all += gg
    
    write_dense_mtx_10x(sid_group_all, data[0][1], mtx_all, out_dire)