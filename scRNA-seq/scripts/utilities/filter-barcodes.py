"""
1) Correction of ambient RNA
Using decontX to correct count data.

2) Doublet Detection
Using scrublet to detect  doublets, and the doublets was removed in subsequent
process.

3) Filtering cells
Cells whose count number great than 500 was used to calculate MAD of count,
gene number, top 20 gene proportion. The cell out of the range of 5 time MAD is
regards as outlayer cells. and cells whose mitochondira count percent great than
 8% is count as outlayer cells. outlayer cells are filtered out. 

"""
import argparse
import os
import subprocess as sp
from io_10x import read_10x_as_dense_mtx
import scrublet as scr
import numpy as np


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("cellranger_out_dire", help="directory of cellranger out")
    parser.add_argument("--mt_file", help="mitochodira gene list file")
    parser.add_argument("--decontx_script", help="file of decontx script")
    parser.add_argument("--out", help="prefix of out files")
    args = parser.parse_args()
    return args.cellranger_out_dire, args.mt_file, args.decontx_script, args.out


def get_mtx_file_path(cellranger_dire, subdir="filtered"):
    if subdir == "filtered":
        barcodes_file = os.path.join(cellranger_dire, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
        features_file = os.path.join(cellranger_dire, "outs", "filtered_feature_bc_matrix", "features.tsv.gz")
        matrix_file = os.path.join(cellranger_dire, "outs", "filtered_feature_bc_matrix", "matrix.mtx.gz")
    else:
        barcodes_file = os.path.join(cellranger_dire, "outs", "raw_feature_bc_matrix", "barcodes.tsv.gz")
        features_file = os.path.join(cellranger_dire, "outs", "raw_feature_bc_matrix", "features.tsv.gz")
        matrix_file = os.path.join(cellranger_dire, "outs", "raw_feature_bc_matrix", "matrix.mtx.gz")

    return barcodes_file, features_file, matrix_file


def correct_ambient(cellranger_dire, decontx_script, out):
    """
    correct ambient using decontX

    Parameters
    ---------------
    cellranger_dire: directory where cellranger result data lives.
    out: prefix of corrected count mtx data file.

    Returns
    --------------
    out: the file name of result file.
    """

    out = out + "_decontX.mtx"
    sp.run(["Rscript", decontx_script, cellranger_dire, out], capture_output=True)
    return out


def detect_doublet(mtx, out):
    mtx_t = mtx.T
    scrub = scr.Scrublet(mtx_t)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    out = out  + "_scrublet.txt"
    fout = open(out, "w")
    for ele in zip(doublet_scores, predicted_doublets):
        print(",".join([str(e) for e in ele]), file=fout)
    fout.close()
    return out


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


def calculate_MAD_boundaries(dt_array, n_mad):
    mad = np.median(np.abs(dt_array - np.median(dt_array)))
    array_median = np.median(dt_array)
    mad_up = array_median + n_mad * mad
    mad_down = array_median - n_mad * mad
    return mad_up, mad_down


def stats_mtx(mtx, mt_idx):
    mtx_sum = np.sum(mtx, axis=0)
    mtx_gene = np.sum(mtx > 0, axis=0)
    mtx_top20_pct = np.sum(np.sort(mtx, axis=0)[-20:, :], axis=0) / mtx_sum
    mtx_mt_num_pct = np.sum(mtx[mt_idx, :], axis=0) / mtx_sum

    print("total count", sum(mtx_sum))
    print("mtx_sum quantiles", np.quantile(mtx_sum, [0, .25, .5, .75, 1]))
    print("mtx_gene_quantile", np.quantile(mtx_gene, [0, .25, .5, .75, 1]))
    print("mtx_top20_pct", np.quantile(mtx_top20_pct, [0, .25, .5, .75, 1]))
    print("mtx_mt_num_pct_quantile", np.quantile(mtx_mt_num_pct, [0, .25, .5, .75, 1]))


def filter_count_mtx(mtx, mt_idx, min_count, out):
    """
    filter count matrix.

    Parameters
    ---------------
    mtx: matrix, row is features, colomn is barcodes

    Returns
    ---------------
    out: file store filtering result. cell passes filter marked as true.

    """
    n_mad = 5
    mtx_sum = np.sum(mtx, axis=0)
    mtx_sum_bool = mtx_sum > min_count
    count_sum_log1p = np.log1p(mtx_sum)
    mtx_gene = np.sum(mtx > 0, axis=0)
    gene_num_log1p = np.log1p(mtx_gene)
    top20_pct = np.sum(np.sort(mtx, axis=0)[-20:, :], axis=0) / mtx_sum
    mt_pct = np.sum(mtx[mt_idx, :], axis=0) / mtx_sum

    count_mad_up, count_mad_down = calculate_MAD_boundaries(count_sum_log1p[mtx_sum_bool], n_mad)
    gene_mad_up, gene_mad_down = calculate_MAD_boundaries(gene_num_log1p[mtx_sum_bool], n_mad)
    top20_pct_mad_up, top20_pct_mad_down = calculate_MAD_boundaries(top20_pct[mtx_sum_bool], n_mad)

    count_bool = (count_sum_log1p < count_mad_up) & (count_sum_log1p >  count_mad_down)
    gene_bool = (gene_num_log1p < gene_mad_up) & (gene_num_log1p > gene_mad_down)
    top20_bool = (top20_pct < top20_pct_mad_up) & (top20_pct > top20_pct_mad_down)
    #here do not filter the mt-pct 20240508
    mt_bool = mt_pct < 1
    filter_bool = mtx_sum_bool & count_bool & gene_bool & top20_bool & mt_bool
    print("not count the doublets")
    print("total count", sum(mtx_sum))
    print("mtx_sum quantiles", np.quantile(mtx_sum[filter_bool], [0, .25, .5, .75, 1]))
    print("mtx_gene_quantile", np.quantile(mtx_gene[filter_bool], [0, .25, .5, .75, 1]))
    print("mtx_top20_pct", np.quantile(top20_pct[filter_bool], [0, .25, .5, .75, 1]))
    print("mtx_mt_num_pct_quantile", np.quantile(mt_pct[filter_bool], [0, .25, .5, .75, 1]))

    out = out + "_filter_bool.txt"
    fout = open(out, "w")
    for ele in filter_bool:
        if ele:
            print("1", file=fout)
        else:
            print("0", file=fout)
    return out


if __name__ == "__main__":
    cellranger_dire, mt_file, decontx_script, out = getargs()
    print(cellranger_dire)

    barcodes_f, features_f, matrix_f = get_mtx_file_path(cellranger_dire, "filtered")
    barcodes, features, mtx = read_10x_as_dense_mtx(barcodes_f, features_f, matrix_f, dtype="int32")
    features = [e[0] for e in features]
    mt_list = [l.rstrip() for l in open(mt_file)]
    mt_idx = []
    for gene in mt_list:
        mt_idx.append(features.index(gene))
    stats_mtx(mtx, mt_idx)

    correct_mtx_f = correct_ambient(cellranger_dire, decontx_script, out)
    doublet_res_f = detect_doublet(mtx, out)

    mtx_docontx = read_decontX_mtx(correct_mtx_f)
    stats_mtx(mtx_docontx, mt_idx)

    #this result is used in next step
    filter_bool_file = filter_count_mtx(mtx_docontx, mt_idx, 500, out + "_count_correct")

