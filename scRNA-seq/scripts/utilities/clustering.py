import warnings
import sys
import scanpy as sc
#import harmonypy as hm
import argparse
import os


def getargs():
    pass
    parser = argparse.ArgumentParser()
    parser.add_argument("adata", help="file name of adata")
    parser.add_argument("--out", help="result file name", default="out")
    parser.add_argument("--revolution", help="revolution of leiden.", default=.5, type=float)

    args = parser.parse_args()
    return args.adata, args.out, args.revolution


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    adata, out, revolution = getargs()
    fout_name = out
    resl = revolution
    adata = sc.read_h5ad(adata)
    #pcs = adata.varm["PCs"]
    sc.external.pp.bbknn(adata, batch_key='sample', n_pcs=40)

    #6hpf using n_neighbors 10, n_pcs 50
    #24hpf using n_neighbors 10, n_pcs 50
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata, resolution=resl)
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False)
    sc.tl.umap(adata, init_pos='paga')
    sc.tl.umap(adata)
    sc.tl.tsne(adata)

    cwd = os.getcwd()
    os.chdir(os.path.dirname(fout_name))
    sc.pl.umap(adata, color=["leiden", "ddx4_nm", "dnd1_nm", "nanos3_nm",\
                            "tdrd7a_nm", "ca15b_nm"], color_map="Reds", save="_" + os.path.basename(fout_name) + "_rv_" + str(revolution) + "_umap.svg", show=False, legend_loc="on data")
    sc.pl.tsne(adata, color=["leiden", "ddx4_nm", "dnd1_nm", "nanos3_nm",\
                            "tdrd7a_nm", "ca15b_nm"], color_map="Reds", save="_" + os.path.basename(fout_name) + "_rv_" + str(revolution) + "_tsne.svg", show=False, legend_loc="on data")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True)
    os.chdir(cwd)
    #sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

    cellids = adata.obs_names
    cluster_info = adata.obs.leiden
    umap_xy=adata.obsm["X_umap"]
    cell_n = len(cellids)
    fout_cell_cluster = open(fout_name + "_rv_" + str(revolution) + "_cell_cluster.tsv", "w")
    for idx in range(cell_n):
        print("\t".join([cellids[idx], str(int(cluster_info.iloc[idx])), str(umap_xy[idx, 0]), str(umap_xy[idx, 1])]), file=fout_cell_cluster)
    fout_cell_cluster.close()


    sig_genes = adata.uns["rank_genes_groups"]
    names_gene = sig_genes["names"]
    scores = sig_genes["scores"]
    pvals = sig_genes["pvals"]
    pvals_adj = sig_genes["pvals_adj"]
    logfoldchanges = sig_genes["logfoldchanges"]
    gene_n = len(names_gene)
    cluster_n = len(adata.uns["leiden_sizes"])
    fout_sig_genes = open(fout_name + "_rv_" + str(revolution) + "_DE_genes.tsv", "w")
    head_line = []
    for i in range(cluster_n):
        head_line += ["cluster_" + str(i) + "_gene_name", "cluster_" + str(i) + "_score",
            "cluster_" + str(i) + "_pval", "cluster_" + str(i) + "_pval_adj", "cluster_" + str(i) + "_log2foldchange"]
    print("\t".join(head_line), file=fout_sig_genes)

    for idx in range(gene_n):
        sig_gene_n = 0
        zip_itm = zip(names_gene[idx], scores[idx], pvals[idx], pvals_adj[idx], logfoldchanges[idx])
        
        line_info = []
        for itm in zip_itm:
            #print(itm)
            gene_name_, score_, pval_, p_adj_, logch_ = itm
            #print(gene_name_, score_, pval_, p_adj_, logch_)
            if p_adj_ < 0.01:
                sig_gene_n += 1
                line_info += [gene_name_, str(score_), str(pval_), str(p_adj_), str(logch_)]
            else:
                line_info += ["", "", "", "", ""]
        print("\t".join(line_info), file=fout_sig_genes)
        if not sig_gene_n:
            break

    fout_sig_genes.close()
    
