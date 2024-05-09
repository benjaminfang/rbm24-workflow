import argparse
import os
import numpy as np
import gzip
from io_10x import write_dense_mtx_10x


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_prefix", help="sample prefix", nargs=4)
    parser.add_argument("data_path", help="10x data directorys", nargs=4)
    parser.add_argument("--mt_genes", help="mitochodira gene list", required=True)
    parser.add_argument("--germ_marker", help="germline gene markers", required=True)
    parser.add_argument("--out", default="out", help="out directory")
    args = parser.parse_args()
    return args.sample_prefix, args.data_path, args.mt_genes, args.germ_marker, args.out


def get_file_name(sample_pre, data_path):
    decontx_file = sample_pre + "_decontX.mtx"
    doublet_file = sample_pre + "_scrublet.txt"
    filter_bool_file = sample_pre + "_count_correct_filter_bool.txt"
    barcodes_file = os.path.join(data_path, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    features_file = os.path.join(data_path, "outs", "filtered_feature_bc_matrix", "features.tsv.gz")
    return decontx_file, doublet_file, filter_bool_file, barcodes_file, features_file


def read_decontX_mtx(mtx_file):
    fin = open(mtx_file, "r")
    fin.readline()
    nrows, ncols, nvals = fin.readline().rstrip().split()
    nrows, ncols, nvals  = int(nrows), int(ncols), int(nvals)
    mtx = np.zeros((nrows, ncols), dtype="float32")
    vals_count = 0 
    for line in fin:
        row, col, val = line.rstrip().split()
        row, col, val = int(row), int(col), float(val)
        mtx[row - 1, col - 1] = val 
        vals_count += 1
    assert vals_count == nvals
    fin.close()
    return mtx


def read_doublet_file(doublet_file):
    dt_out_bool = []
    fin = open(doublet_file, "r")
    for line in fin:
        line = line.rstrip().split(",")[1]
        #Fasle means not doublets, some append true to keep the barcodes
        if line == "False":
            dt_out_bool.append(1)
        else:
            dt_out_bool.append(0)
    dt_out_bool = np.bool_(dt_out_bool)
    return dt_out_bool


def read_filter_bool_file(filter_bool_file):
    dt_out_bool = []
    fin = open(filter_bool_file, "r")
    for line in fin:
        line = line.rstrip()
        if line == "0":
            dt_out_bool.append(0)
        else:
            dt_out_bool.append(1)
    dt_out_bool = np.bool_(dt_out_bool)
    return dt_out_bool


def read_barcodes_features(barcodes_f, features_f):
    barcodes = []
    features = []
    fin = gzip.open(barcodes_f, "r")
    for line in fin:
        barcodes.append(line.decode().rstrip())
    fin.close()

    fin = gzip.open(features_f, "r")
    for line in fin:
        features.append(line.decode().rstrip().split("\t"))
    fin.close()

    return barcodes, features


def combin_file_name(sample_prefixs):
    return "_".join([os.path.basename(e) for e in sample_prefixs])


def check_features_consistence(features_s):
    print(len(features_s))
    print(len(features_s[0]))

    for idx in range(len(features_s[0])):
        tmp = set()
        for jdx in range(len(features_s)):
            tmp.add(features_s[jdx][idx][0])
        if len(tmp) > 1:
            print("features not consistent")
            exit()


if __name__ == "__main__":
    sample_prefixs, data_paths, mt_genes, germ_markers, out = getargs()
    mt_genes = [line.rstrip() for line in open(mt_genes)]
    germ_markers = [line.rstrip() for line in open(germ_markers)]
    print(mt_genes)
    print(germ_markers)
    barcodes_s = []
    features_s = []
    mtx_s = []
    com_dire = combin_file_name(sample_prefixs)
    com_dire = os.path.join(out, com_dire)
    if not os.path.exists(com_dire):
        os.mkdir(com_dire)
    #print(com_dire)

    for sample, data in zip(sample_prefixs, data_paths):
        decontx_file, doublet_file, filter_bool_file, barcodes_file,\
            features_file = get_file_name(sample, data)
        mtx = read_decontX_mtx(decontx_file)
        doublet_bool = read_doublet_file(doublet_file)
        filter_bool = read_filter_bool_file(filter_bool_file)
        barcodes, features = read_barcodes_features(barcodes_file, features_file)
        keep_bool = doublet_bool & filter_bool
        print(len(barcodes))
        print(len(keep_bool))
        print(sum(keep_bool))
        
        barcodes_keep = np.asarray(barcodes)[keep_bool]
        features_keep = features
        mtx_keep = mtx.T[keep_bool].T
        barcodes_s.append(barcodes_keep)
        features_s.append(features_keep)
        mtx_s.append(mtx_keep)
    
    com_barcodes = []
    com_features = []
    com_mtx = []
    for idx, bar in enumerate(barcodes_s):
        sid = os.path.basename(sample_prefixs[idx])
        bar = [sid + "_" + ele.split("-")[0] for ele in bar]
        com_barcodes += list(bar)

    check_features_consistence(features_s)
    com_features = features_s[0]
    mt_idx = []
    germ_markers_idx = []
    germ_markers_align = []
    for idx in range(len(com_features)):
        gene_id = com_features[idx][0]
        if gene_id in mt_genes:
            mt_idx.append(idx)
        if gene_id in germ_markers:
            germ_markers_idx.append(idx)
            germ_markers_align.append(gene_id)

    com_mtx = np.hstack(mtx_s)
    
    stat_file = open(com_dire + "_barcodes_stat.txt", "w")
    print("\t".join(["barcodes", "count_sum", "gene_num", "top20_pct", "mt_pct"] + germ_markers_align), file=stat_file)
    for idx in range(len(com_barcodes)):
        bar = com_barcodes[idx]
        row = com_mtx[:, idx]
        count_num = sum(row)
        gene_num = sum(row > 0)
        top20_pct = sum(np.sort(row)[-20:]) / count_num
        mt_pct = sum(row[mt_idx]) / count_num
        germ_markers_count = row[germ_markers_idx]
        print("\t".join([bar] + [str(count_num), str(gene_num)] + \
                        [format(top20_pct, ".5f"), format(mt_pct, ".5f")] + \
                        [str(e) for e in germ_markers_count]), \
                        file=stat_file)
    stat_file.close()

    write_dense_mtx_10x(com_barcodes, com_features, com_mtx, com_dire)


