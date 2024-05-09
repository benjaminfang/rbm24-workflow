from io_10x import read_10x_as_dense_mtx
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
import sys
import os
import numpy as np
import argparse


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


def keep_corr_gene(mtx):
    keep = np.bool_(np.zeros(len(mtx)))
    idx = np.arange(len(mtx))
    left = mtx
    idx_left = idx
    while True:
        if len(left) < 2:
            break

        item = left[0]
        i = 1
        left_tmp = []
        idx_tmp = []
        for raw in left[1:]:
            if abs(np.corrcoef(raw, item)[0, 1]) > .2:
                keep[idx_left[i]] = True
                flg = 1
            else:
                left_tmp.append(raw)
                idx_tmp.append(idx_left[i])
            i += 1
        if flg:
            keep[idx[0]] = True

        left = left_tmp
        idx_left = idx_tmp
    return keep


def keep_corr_gene2(mtx):
    corr_bool = np.bool_(np.zeros(len(mtx)))
    left = list(range(len(mtx)))
    while True:
        if len(left) <= 1:
            break
        g1_idx = left[0]
        left = left[1:]
        left_2 = []
        g1_vec = mtx[g1_idx]
        flg = 0 
        for idx in left:
            g2_vec = mtx[idx]
            if abs(np.corrcoef(g1_vec, g2_vec)[0][1]) > 0.2:
                corr_bool[idx] = True
                flg = 1 
            else:
                left_2.append(idx)
            if flg:
                corr_bool[g1_idx] = True
        left = left_2

    return corr_bool


def remove_cell_cycle(mtx, cycle_mtx):
    keep = np.bool_(np.ones(len(mtx)))
    idx = np.arange(len(mtx))
    mtx_left = mtx
    while True:
        if len(cycle_mtx) < 1:
            break
        tmp_cycle = []
        for raw in cycle_mtx:
            tmp_left = []
            tmp_idx = []
            i = 0
            for raw2 in mtx_left:
                if abs(np.corrcoef(raw, raw2)[0, 1]) > .4:
                    keep[idx[i]] = False
                    tmp_cycle.append(raw2)
                else:
                    tmp_left.append(raw2)
                    tmp_idx.append(idx[i])
                i += 1
            mtx_left = tmp_left
            idx = tmp_idx
        cycle_mtx = tmp_cycle
    return keep


def remove_cell_cycle2(mtx, cycle_mtx):
    coor_bool = np.bool_(np.ones(len(mtx)))
    left = np.arange(len(mtx))
    while True:
        if len(cycle_mtx) < 1:
            break
        cycle_tmp = []
        for raw in cycle_mtx:
            left_tmp = []
            for idx in left:
                if abs(np.corrcoef(raw, mtx[idx])[0, 1]) > .4:
                    coor_bool[idx] = False
                    cycle_tmp.append(mtx[idx])
                else:
                    left_tmp.append(idx)
            left = left_tmp
        cycle_mtx = cycle_tmp
    return coor_bool


def mt_gene_cho(var_gene, mt_gene):
    dt_out = []
    for g in var_gene:
        if g in mt_gene:
            dt_out.append(True)
        else:
            dt_out.append(False)
    return dt_out


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("mtx_dire", help="directory of data")
    parser.add_argument("--germ_marker", help="germcel marker genes", required=True)
    parser.add_argument("--mt_genes", help="miticondria gene list", required=True)
    parser.add_argument("--cell_cycle", help="cell cycle and house keeping genes list", required=True)
    parser.add_argument("--out", help="file name prefix of result file", default="out")

    args = parser.parse_args()
    return args.mtx_dire, args.germ_marker, args.mt_genes, args.cell_cycle, args.out


if __name__ == "__main__":
    mtx_dire, germ_marker, mt_genes, cell_cycle, out = getargs()
    bar_f = os.path.join(mtx_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(mtx_dire, "features.tsv.gz")
    mtx_f = os.path.join(mtx_dire, "matrix.mtx.gz")
    germ_marker = [l.rstrip() for l in open(germ_marker)]
    cell_cycle = [l.rstrip() for l in open(cell_cycle)]
    mt_genes = [l.rstrip() for l in open(mt_genes)]
    gene_fout = out

    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")

    gene_by_cell_lg = filter_gene_have_cell(mtx)
    mtx = mtx[gene_by_cell_lg]
    features = list(np.asarray(features)[gene_by_cell_lg])
    mtx = norm_mtx(mtx)

    
    cell_cycle_idx = []
    var_names = [e[0] for e in features]
    for gene in var_names:
        if gene in cell_cycle:
            cell_cycle_idx.append(True)
        else:
            cell_cycle_idx.append(False)

    mtx_cycle = mtx[cell_cycle_idx]
    mtx = mtx[np.logical_not(cell_cycle_idx)]
    features = list(np.asarray(features)[np.logical_not(cell_cycle_idx)])
    var_names = [e[0] for e in features]

    count = csr_matrix(mtx.T)
    adata = ad.AnnData(count)
    adata.obs_name = barcodes
    adata.var_names = var_names
    sc.pp.log1p(adata)
    #adata.obs["sample"] = [e.split("_")[0] for e in barcodes]
    #print("here")
    sc.pp.highly_variable_genes(adata, min_mean=0.0001, max_mean=3, min_disp=0.05, n_top_genes=2000, flavor="seurat")
    #print("here")
    hvg_lg = np.asarray([e for e in adata.var.highly_variable])
    gene_sel = np.asarray([e for e in adata.var_names])
    gene_sel = gene_sel[hvg_lg]
    for g in germ_marker:
        if g not in gene_sel:
            print(g)
    mtx_hvg = mtx[hvg_lg]
    mt_gene_bool = mt_gene_cho(gene_sel, mt_genes)
    print(sum(mt_gene_bool))
    gene_sel = gene_sel[np.logical_not(mt_gene_bool)]
    mtx_hvg = mtx_hvg[np.logical_not(mt_gene_bool)]

    keep = keep_corr_gene2(mtx_hvg)
    print("corr .2", np.sum(keep))
    mtx_hvg = mtx_hvg[keep]
    gene_sel = gene_sel[keep]
    keep = remove_cell_cycle2(mtx_hvg, mtx_cycle)
    print("not cycle", np.sum(keep))
    mtx_hvg = mtx_hvg[keep]
    gene_sel = gene_sel[keep]
    fout = open(out + "_hvg.txt", "w")
    for gene in gene_sel:
        print(gene, file=fout)
    fout.close()


