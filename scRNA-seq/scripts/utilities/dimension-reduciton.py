import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
import sys 
import os
import numpy as np
from io_10x import read_10x_as_dense_mtx
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


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("mtx_dire", help="data directory")
    parser.add_argument("hvg", help="high variable genes")
    parser.add_argument("--germ_marker", help="germcel marker genes", required=True)
    parser.add_argument("--mt_genes", help="miticondria gene list", required=True)
    parser.add_argument("--max_scale_value", help="max value of pp.scale", default=10, type=int)
    parser.add_argument("--out", help="file name prefix of result file", default="out")

    args = parser.parse_args()
    return args.mtx_dire, args.hvg, args.germ_marker, args.mt_genes, args.out, args.max_scale_value

if __name__ == "__main__":
    mtx_dire, hvg, germ_marker, mt_genes, out, max_scale_value = getargs()
    germ_gene_f = germ_marker
    hvg_f = hvg
    mt_gene_f = mt_genes
    out_h5ad = out + ".h5ad"
    max_v = max_scale_value
    bar_f = os.path.join(mtx_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(mtx_dire, "features.tsv.gz")
    mtx_f = os.path.join(mtx_dire, "matrix.mtx.gz")

    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")
    hvg = [l.rstrip() for l in open(hvg_f)]
    germ_gene = [l.rstrip() for l in open(germ_gene_f)]
    mt_gene = [l.rstrip() for l in open(mt_gene_f)]

    keep = filter_gene_have_cell(mtx)
    mtx = mtx[keep]
    features = list(np.asarray(features)[keep])
    feat_g = [e[0] for e in features]
    for g in germ_gene:
        if g not in feat_g:
            print(g)
    
    count = csr_matrix(mtx.T)
    adata = ad.AnnData(count)
    adata.var_names = [e[0] for e in features]
    adata.obs_names = barcodes
    adata.obs["sample"] = [e.split("_")[0] for e in barcodes]

    mt_bool = []
    for gene in adata.var_names:
        if gene in mt_gene:
            mt_bool.append(True)
        else:
            mt_bool.append(False)
    adata.var["mt"] = np.array(mt_bool)

    sel_gene  = hvg + germ_gene
    #sel_gene  = hvg
    print(len(set(sel_gene)))
    sel_bool = []
    for gene in adata.var_names:
        if gene in sel_gene:
            sel_bool.append(True)
        else:
            sel_bool.append(False)
    print("select gene", sum(sel_bool))

    adata.var["highly_variable"] = np.array(sel_bool)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.raw = adata

    adata.obs["nanos3_nm"] = adata.obs_vector("nanos3")
    adata.obs["ddx4_nm"] = adata.obs_vector("ddx4")
    adata.obs["dnd1_nm"] = adata.obs_vector("dnd1")
    adata.obs["tdrd7a_nm"] = adata.obs_vector("tdrd7a")
    adata.obs["ca15b_nm"] = adata.obs_vector("ca15b")

    sc.pp.log1p(adata)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=max_v)

    sc.tl.pca(adata, svd_solver='arpack', use_highly_variable="highly_variable")
    print(adata)
    adata.write(out_h5ad)
    
