import argparse
import matplotlib.pyplot as plt
import numpy as np
import copy
from io_10x import read_10x_as_dense_mtx
import os
import matplotlib as mpl


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("cluster_data", help="cluster data file, contained barcodes, cluster id, umapx, umapy columns")
    parser.add_argument("cluster_info", help="cluster information, annotated the cell type for each cluster id")
    parser.add_argument("--marker_genes", help="cluster marker genes for cluster expression")
    parser.add_argument("--count_dire", help="count directory")
    parser.add_argument("--heatmap_genes", help="list of heat map genes")
    parser.add_argument("--germ_marker", help="germcell marker gene list")
    parser.add_argument("--scatter_genes", help="marker genes for scatter plot")
    parser.add_argument("--violin_genes", help="gene list for violin plot")
    parser.add_argument("--ctl", help="control group ids", required=True, type=str)
    parser.add_argument("--exp", help="experiment group ids", required=True, type=str)
    parser.add_argument("--out", help="file name of the results", default="out")

    args = parser.parse_args()
    return (args.cluster_data, args.cluster_info, args.ctl, args.exp,
            args.out, args.marker_genes, args.count_dire,
            args.heatmap_genes, args.germ_marker, args.scatter_genes, args.violin_genes)


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


def read_cluster_data(cluster_data_file):
    cluster_data = {}
    fin = open(cluster_data_file)
    for line in fin:
        cellid, clusterid, umapx, umapy = line.rstrip().split("\t")
        umapx = float(umapx)
        umapy = float(umapy)
        sid = cellid.split("_")[0]
        if sid in cluster_data:
            if clusterid in cluster_data[sid]:
                cluster_data[sid][clusterid].append([cellid, umapx, umapy])
            else:
                cluster_data[sid][clusterid] = [[cellid, umapx, umapy]]
        else:
            cluster_data[sid] = {clusterid: [[cellid, umapx, umapy]]}

    return cluster_data


def read_cluster_info(cluster_info_file):
    cluster_info = {}
    fin = open(cluster_info_file, "r")
    for line in fin:
        clu, clu_name, color, alpha = line.rstrip().split(",")
        cluster_info[clu] = [clu_name, color, float(alpha)]
    return cluster_info


def plot_cluster(ax, x, y, color, alpha, label):
    size = 10
    if label == "Unknown":
        ax.scatter(x, y, c=color, alpha=alpha, label=label, s=size, linewidths=0, zorder=0)
    else:
        ax.scatter(x, y, c=color, alpha=alpha, label=label, s=size, linewidths=0)


def plot_scatter(cluster_data_file, cluster_info_file, ctl, exp, figout):
    """
    Plot barcode umap plot and annotated cell types.
    """
    cluster_data = read_cluster_data(cluster_data_file)
    cluster_info = read_cluster_info(cluster_info_file)
    ctl = ctl.split(",")
    exp = exp.split(",")

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18.5, 5), layout="constrained")

    cluster_plot_all = {}
    cluster_plot_ctl = {}
    cluster_plot_exp = {}
    for sid in cluster_data:
        for cluster in cluster_data[sid]:
            cluster_data_this = cluster_data[sid][cluster]
            if cluster not in cluster_plot_all:
                cluster_plot_all[cluster] = [[], []]
            for cell in cluster_data_this:
                cluster_plot_all[cluster][0].append(cell[1])
                cluster_plot_all[cluster][1].append(cell[2])

        if sid in ctl:
            for cluster in cluster_data[sid]:
                cluster_data_this = cluster_data[sid][cluster]
                if cluster not in cluster_plot_ctl:
                    cluster_plot_ctl[cluster] = [[], []]
                for cell in cluster_data_this:
                    cluster_plot_ctl[cluster][0].append(cell[1])
                    cluster_plot_ctl[cluster][1].append(cell[2])

        if sid in exp:
            for cluster in cluster_data[sid]:
                cluster_data_this = cluster_data[sid][cluster]
                if cluster not in cluster_plot_exp:
                    cluster_plot_exp[cluster] = [[], []]
                for cell in cluster_data_this:
                    cluster_plot_exp[cluster][0].append(cell[1])
                    cluster_plot_exp[cluster][1].append(cell[2])


    for cluster in cluster_plot_all:
        label, color, alpha = cluster_info[cluster]
        plot_cluster(axs[0], cluster_plot_all[cluster][0], cluster_plot_all[cluster][1], color, alpha, label)

        
    for cluster in cluster_plot_ctl:
        label, color, alpha = cluster_info[cluster]
        plot_cluster(axs[1], cluster_plot_ctl[cluster][0], cluster_plot_ctl[cluster][1], color, alpha, label)

    for cluster in cluster_plot_exp:
        label, color, alpha = cluster_info[cluster]
        plot_cluster(axs[2], cluster_plot_exp[cluster][0], cluster_plot_exp[cluster][1], color, alpha, label)
    

    xlim = axs[0].get_xlim()
    ylim = axs[0].get_ylim()
    axs[0].set_title("Control and Mrbm24a")
    axs[1].set_title("Control")
    axs[2].set_title("Mrbm24a")
    for ax in axs:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        #ax.set_xticks([], [])
        #ax.set_yticks([], [])
        ax.set_xlabel("UMAP_1")
        ax.set_ylabel("UMAP_2")
    handles, labels = axs[0].get_legend_handles_labels()
    
    lab_my = ["Ectoderm", "Lateral and Ventral Mesoderm", "Dorsal Organizer", "Endoderm", "EVL", "Dorsal Forerunner Cell", "PGC", "Unknown"]
    han_my = []
    h_l_dic = {}
    for h, l in zip(handles, labels):
        h_l_dic[l] = h
    for e in lab_my:
        han_my.append(h_l_dic[e])
    leg = fig.legend(handles=han_my, labels=lab_my, loc="outside right upper")
    for ll in leg.legend_handles:
        ll.set_sizes([40])
        ll.set_alpha(1)
    
    fig.savefig(figout + ".svg")


def output_cell_num(cell_num, cluster_info, figout):
    fout_num = open(figout + "_cell_num.csv", "w")
    fout_pct = open(figout + "_cell_pct.csv", "w")
    cell_sum = {}
    for sid in cell_num:
        cell_sum[sid] = sum([cell_num[sid][cid] for cid in cell_num[sid]])
    head = ["cluster_id", "cluster_name", "s95", "s96", "s97", "s98"]
    print(",".join(head), file=fout_num)
    print(",".join(head), file=fout_pct)
    sid_ord = ["s95", "s96", "s97", "s98"]
    cluster_ord = [0, 1, 4, 6, 7, 8, 9, 2, 3, 5]
    cluster_ord = [str(e) for e in cluster_ord]

    for cid in cluster_ord:
        print(",".join([cid, cluster_info[cid][0]] + [format(cell_num[sid][cid], ".5f") for sid in sid_ord]), file=fout_num)
        print(",".join([cid, cluster_info[cid][0]] + [format(cell_num[sid][cid] / cell_sum[sid], ".5f") for sid in sid_ord]), file=fout_pct)
        
    fout_num.close()
    fout_pct.close()


def bar_cluster_celln(cell_num, cluster_info, figout):
    cell_num = copy.deepcopy(cell_num)

    gid_order = ["s95", "s96", "s97", "s98"]
    color_dic = {"s95": "#1f77b4", "s96": "#1f77b4", "s97": "#ff7f0e", "s98": "#ff7f0e"}
    label_dic = {"s95": "Control", "s96": "Control", "s97": "Mrbm24a", "s98": "Mrbm24a"}
    cluster_ord = [0, 1, 4, 6, 7, 8, 9, 2, 3, 5]
    cluster_ord = [str(e) for e in cluster_ord]
    fig, ax = plt.subplots(layout="constrained")
    width = .2
    offset = .2
    offset_t = 0
    sid_sum = {}
    for sid in cell_num:
        sid_sum[sid] = sum([cell_num[sid][k] for k in cell_num[sid]])
    for sid in gid_order:
        cn = []
        for cid in cluster_ord:
            cn.append(cell_num[sid][cid])
        #combine unknowns
        last_2 = sum(cn[-2:])
        cn = cn[:-2]
        cn[-1] += last_2
        x = np.linspace(-.3, 6.7, 8)
        cn = np.asarray(cn) / sid_sum[sid]
        ax.bar(x + offset * offset_t, cn, width=width, label=label_dic[sid], color=color_dic[sid])
        offset_t += 1
    cname = []
    for cid in cluster_ord:
        cname.append(cluster_info[cid][0])
    cname = cname[:-2]
    print(cname)
    loc = list(np.arange(8))
    print(loc)
    xticks = ax.set_xticks(loc, cname, rotation=90)
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[0], handles[2]]
    labels = [labels[0], labels[2]]
    fig.legend(handles=handles, labels=labels)
    fig.savefig(figout + "_num_bar.svg")


def bar_stack_pct(cell_num, cluster_info, figout):
    cell_num = copy.deepcopy(cell_num)
    
    gid_order = ["s95", "s96", "s97", "s98"]
    cluster_ord = [0, 1, 4, 6, 7, 8, 9, 2, 3, 5]
    cluster_ord = [str(e) for e in cluster_ord]

    cell_num_samp = []
    for sid in cell_num:
        tmp = []
        for cid in cluster_ord:
            tmp.append(cell_num[sid][cid])
        tmp[-3] += sum(tmp[-2:])
        cell_num_samp.append(tmp[:-2])
    print(cell_num_samp)
    #experiment, control
    tmp = [None, None]
    tmp[0] = list(np.asarray(cell_num_samp[0]) + np.asarray(cell_num_samp[1]))
    tmp[1] = list(np.asarray(cell_num_samp[2]) + np.asarray(cell_num_samp[3]))

    print(tmp)
    cell_pct = np.asarray(tmp)
    print(cell_pct)
    cell_pct = cell_pct.T / np.sum(cell_pct, axis=1)
    print(cell_pct)

    fig, ax = plt.subplots(figsize=(4, 16), layout="constrained")
    x = ["Control", "Mrbm24a"]
    bottom = np.zeros(2)
    idx = 0
    for row in cell_pct:
        ax.bar(x, row, bottom=bottom, width=.5, label=cluster_info[cluster_ord[idx]][0])
        bottom += row
        idx += 1
    xtk = ax.get_xticklabels()
    for xx in xtk:
        xx.set_rotation(90)
    hands, labels = ax.get_legend_handles_labels()

    leg = fig.legend(handles=hands, labels=labels, loc="outside right upper")

    fig.savefig(figout + "_bar_stack.svg") 


def plot_cloumn(cluster_data_file, cluster_info_file, ctl, exp, figout):
    """
    Plot bar figure for cell number in each clusters.
    """
    cluster_data = read_cluster_data(cluster_data_file)
    cluster_info = read_cluster_info(cluster_info_file)
    ctl = ctl.split(",")
    exp = exp.split(",")

    cell_num = {}
    cluster_id = set()
    for sid in cluster_data:
        cell_num[sid] = {}
        for cid in cluster_data[sid]:
            cell_num[sid][cid] = len(cluster_data[sid][cid])
            cluster_id.add(cid)
    for cid in cluster_id:
        for gid in cell_num:
            if cid not in cell_num[gid]:
                cell_num[gid][cid] = 0
    print(cell_num)
    bar_cluster_celln(cell_num, cluster_info, figout)
    output_cell_num(cell_num, cluster_info, figout)
    bar_stack_pct(cell_num, cluster_info, figout)


def read_marker_genes(marker_genes):
    dt = {}
    fin = open(marker_genes, "r")
    key = None
    for line in fin:
        line = line.rstrip()
        if line[0] == ">":
            key = line[1:]
            dt[key] = []
        else:
            dt[key].append(line)
    fin.close()
    return dt


def get_trim_mtx(features, marker_genes_dt):
    marker_genes = []
    for cid in marker_genes_dt:
        marker_genes += marker_genes_dt[cid]
    trim_bool = []
    for fea in features:
        if fea[0] in marker_genes:
            trim_bool.append(True)
        else:
            trim_bool.append(False)
    return trim_bool


def struct_count_dt(barcodes, features, mtx):
    dt = {}
    genes = [fea[0] for fea in features]
    mtx_t = mtx.T
    idx = 0
    for row in mtx_t:
        dt[barcodes[idx]] = {ele[0]: ele[1] for ele in zip(genes, row)}
        idx += 1
    return dt


def combin_sample(cluster_data, ctl, exp):
    dt_out = {"ctl": {}, "exp": {}}
    for sid in cluster_data:
        for cid in cluster_data[sid]:
            if sid in ctl:
                if cid not in dt_out["ctl"]:
                    dt_out["ctl"][cid] = cluster_data[sid][cid]
                else:
                    dt_out["ctl"][cid] += cluster_data[sid][cid]
            elif sid in exp:
                if cid not in dt_out["exp"]:
                    dt_out["exp"][cid] = cluster_data[sid][cid]
                else:
                    dt_out["exp"][cid] += cluster_data[sid][cid]
    return dt_out


def query_gene_count_data(data, geneid, barcodes_gene_count):
    dt_out = {}
    for cid in data:
        dt_out[cid] = []
        for item in data[cid]:
            barcodes, umapx, umapy = item
            dt_out[cid].append([barcodes, umapx, umapy, barcodes_gene_count[barcodes][geneid]])
    return dt_out


def plot_marker_genes_cluster(fig, ax, x, y, c, geneid, vmax):
    cb = ax.scatter(x, y, c=c, s=3, cmap="Reds", vmin=0, vmax=vmax)
    ax.set_title(geneid)
    fig.colorbar(cb, ax=ax)


def plot_marker_genes(cluster_data_file, cluster_info_file, ctl, exp, figout, marker_genes, count_dire):
    """
    Read count data (cellrange gz format)which may is original count
    or normalized.
    plot scatter the color is the expression level of marker genes.
    """
    bar_f = os.path.join(count_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(count_dire, "features.tsv.gz")
    mtx_f = os.path.join(count_dire, "matrix.mtx.gz")
    ctl = ctl.split(",")
    exp = exp.split(",")
    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")
    filter_bool = filter_gene_have_cell(mtx)
    mtx = mtx[filter_bool]
    features = list(np.asarray(features)[filter_bool])
    mtx = norm_mtx(mtx)

    marker_genes_dt = read_marker_genes(marker_genes)
    print(marker_genes_dt)
    trim_bool = get_trim_mtx(features, marker_genes_dt)
    mtx = mtx[trim_bool]
    features = list(np.asarray(features)[trim_bool])

    barcodes_gene_count = struct_count_dt(barcodes, features, mtx)
    del mtx
    cluster_data = read_cluster_data(cluster_data_file)
    cluster_info = read_cluster_info(cluster_info_file)
    print(cluster_info)

    cluster_data_com = combin_sample(cluster_data, ctl, exp)

    ctl_dt = cluster_data_com["ctl"]
    exp_dt = cluster_data_com["exp"]
    for cid in marker_genes_dt:
        ncols = 2
        nrows = len(marker_genes_dt[cid])
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4.5 * ncols, 4 * nrows), layout="constrained")
        idx = 0
        fig_name = figout + "_" + cid + ".svg"
        for geneid in marker_genes_dt[cid]:
            ctl_dt_gene = query_gene_count_data(ctl_dt, geneid, barcodes_gene_count)
            exp_dt_gene = query_gene_count_data(exp_dt, geneid, barcodes_gene_count)
            x_ctl = []
            y_ctl = []
            c_ctl = []
            for cid_2 in ctl_dt_gene:
                for item in ctl_dt_gene[cid_2]:
                    x_ctl.append(item[1])
                    y_ctl.append(item[2])
                    c_ctl.append(item[3])
            x_exp = []
            y_exp = []
            c_exp = []
            for cid_3 in exp_dt_gene:
                for item in exp_dt_gene[cid_3]:
                    x_exp.append(item[1])
                    y_exp.append(item[2])
                    c_exp.append(item[3])
            vmax = max(c_ctl + c_exp)
            plot_marker_genes_cluster(fig, axs[idx, 0], x_ctl, y_ctl, c_ctl, geneid, vmax)
            plot_marker_genes_cluster(fig, axs[idx, 1], x_exp, y_exp, c_exp, geneid, vmax)
            idx += 1
        print(cluster_info[cid])
        fig.suptitle(cluster_info[cid][0] + " sibling vs Mrbm24a(right column)")
        fig.savefig(fig_name)


def select_genes_idx(features, headmap_genes_list):
    idx_list = []
    gene_list = [fea[0] for fea in features]
    for gene in headmap_genes_list:
        idx_list.append(gene_list.index(gene))    
    return idx_list


def choose_barcodes(barcodes, barcodes_list, mtx):
    dt_out = []
    mtx_t = mtx.T
    for block in barcodes_list:
        idx = []
        for bar in block:
            idx.append(barcodes.index(bar))
        dt_out.append(mtx_t[idx].T)
    return dt_out


def scale_mtx_by_gene(mtx):
    #zero to max scaled to 0 to one
    mtx = mtx.T / np.max(mtx, axis=1)
    mtx = mtx.T
    mtx = np.log1p(mtx * 100)

    return mtx


def plot_heatmap_fuc(cluster_ids, heat_dt, figout, surfix, cluster_info):
    width_r = []
    for block in heat_dt:
        width_r.append(block.shape[1])
    fig, axs = plt.subplots(ncols=len(cluster_ids), figsize=(2 * len(cluster_ids), 8), layout="constrained")
    idx = 0
    for cid in cluster_ids:
        print(heat_dt[idx].shape)
        cb = axs[idx].pcolormesh(heat_dt[idx][::-1], cmap="viridis")
        axs[idx].set_xticks([], [])
        axs[idx].set_yticks([], [])
        axs[idx].set_title(str(cid))
        idx += 1
    cmap = mpl.cm.viridis
    norm= mpl.colors.Normalize(vmin=np.log1p(1), vmax=np.log1p(100))
    fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm), ax=axs)
    fig.savefig(figout + "_" + surfix + "_heatmap.tiff")


def plot_heatmap(cluster_data_file, cluster_info_file, ctl, exp, figout, count_dire, headmap_genes):
    """
    plot heatmap the marker gene expression in barcodes.
    """
    cluster_data = read_cluster_data(cluster_data_file)
    cluster_info = read_cluster_info(cluster_info_file)
    ctl = ctl.split(",")
    exp = exp.split(",")
    heatmap_genes = read_marker_genes(headmap_genes)
    bar_f = os.path.join(count_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(count_dire, "features.tsv.gz")
    mtx_f = os.path.join(count_dire, "matrix.mtx.gz")
    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")
    filter_bool = filter_gene_have_cell(mtx)
    mtx = mtx[filter_bool]
    features = list(np.asarray(features)[filter_bool])
    mtx = norm_mtx(mtx)

    heatmap_cluster = [0, 1, 4, 6, 7, 8, 9]
    heatmap_cluster = [str(ele) for ele in heatmap_cluster]
    heatmap_genes_list = []
    
    for cid in heatmap_cluster:
        heatmap_genes_list += heatmap_genes[cid]
    cluster_data_com = combin_sample(cluster_data, ctl, exp)
    barcodes_list_ctl = []
    barcodes_list_exp = []
    barcodes_list_ctl_exp = []
    ctl_cid = []
    exp_cid = []
    ctl_exp_cid = []
    for cid in heatmap_cluster:
        ctl_cid_bar = cluster_data_com["ctl"].get(cid)
        exp_cid_bar = cluster_data_com["exp"].get(cid)
        tmp = []
        ctl_exp_cid.append(cid)
        if ctl_cid_bar:
            ctl_cid.append(cid)
            print("ctl cluster", cid)
            barcodes_list_ctl.append([e[0] for e in ctl_cid_bar])
            tmp += [e[0] for e in ctl_cid_bar]
        if exp_cid_bar:
            exp_cid.append(cid)
            print("exp cluster", cid)
            barcodes_list_exp.append([e[0] for e in exp_cid_bar])
            tmp += [e[0] for e in exp_cid_bar]
        barcodes_list_ctl_exp.append(tmp)


    gene_idx = select_genes_idx(features, heatmap_genes_list)
    mtx = mtx[gene_idx]
    mtx = scale_mtx_by_gene(mtx)
    features = np.asarray(features)[gene_idx]
    heat_ctl = choose_barcodes(barcodes, barcodes_list_ctl, mtx)
    heat_exp = choose_barcodes(barcodes, barcodes_list_exp, mtx)
    heat_ctl_exp = choose_barcodes(barcodes, barcodes_list_ctl_exp, mtx)
    plot_heatmap_fuc(ctl_cid, heat_ctl, figout, "ctl", cluster_info)
    plot_heatmap_fuc(exp_cid, heat_exp, figout, "exp", cluster_info)
    plot_heatmap_fuc(ctl_exp_cid, heat_ctl_exp, figout, "ctl_exp", cluster_info)


def split_none_germ_germ(ctl_cid_bar):
    none_germ_cid = [str(e) for e in range(9)]
    germ_cid = ["9"]
    none_germ_bars = []
    germ_bars = []
    for cid in none_germ_cid:
        none_germ_bars += ctl_cid_bar[cid]
    for cid in germ_cid:
        germ_bars += ctl_cid_bar[cid]
    return [e[0] for e in none_germ_bars], [e[0] for e in germ_bars]    


def gene_count_stat(bars, barcodes, mtx):
    choose_bool = []
    for bb in barcodes:
        if bb in bars:
            choose_bool.append(True)
        else:
            choose_bool.append(False)
    mtx = mtx.T[choose_bool]
    sums = np.sum(mtx, axis=0)
    return sums


def compare_marker_gene(cluster_data_file, cluster_info_file, ctl, exp, figout, count_dire, germ_marker):
    """
    compare germcell marker count between non-PGCs and PGCs
    """
    cluster_data = read_cluster_data(cluster_data_file)
    cluster_info = read_cluster_info(cluster_info_file)
    ctl = ctl.split(",")
    exp = exp.split(",")
    bar_f = os.path.join(count_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(count_dire, "features.tsv.gz")
    mtx_f = os.path.join(count_dire, "matrix.mtx.gz")
    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")
    filter_bool = filter_gene_have_cell(mtx)
    mtx = mtx[filter_bool]
    features = list(np.asarray(features)[filter_bool])
    mtx = norm_mtx(mtx)
    germ_marker = [l.rstrip() for l in open(germ_marker)]
    print(germ_marker)
    
    """    
    germ_marker_idx = []
    germ_marker_align = []
    for i in range(len(features)):
        if features[i][0] in germ_marker:
            germ_marker_idx.append(i)
            germ_marker_align.append(features[i][0])
    print(germ_marker_align)
    mtx = mtx[germ_marker_idx]

    mtx = mtx.T
    i = 0
    for b in barcodes:
        print(b, mtx[i])
        i += 1
    """
    
    cluster_data_com = combin_sample(cluster_data, ctl, exp)
    ctl_cid_bar = cluster_data_com["ctl"]
    print(germ_marker)
      
    gene_idx = select_genes_idx(features, germ_marker)
    mtx = mtx[gene_idx]
    print(mtx.shape)
    features = np.asarray(features)[gene_idx]
    none_germ_bars, germ_bars = split_none_germ_germ(ctl_cid_bar)
    
    print(germ_bars)
    none_germ_count = gene_count_stat(none_germ_bars, barcodes, mtx)
    germ_count = gene_count_stat(germ_bars, barcodes, mtx)
    print(list(germ_count))
    print(list(none_germ_count))
    print(germ_marker)
    print("\t".join(["geneid", "PGC", "non-PGC"]))
    for i in range(len(germ_count)):
        print("\t".join([germ_marker[i], str(germ_count[i]), str(none_germ_count[i])]))


    fig, ax = plt.subplots(layout="constrained")
    x = np.arange(len(germ_marker))
    ax.bar(x - .2, germ_count, width=.4, label="PGCs")
    ax.bar(x + .2, none_germ_count, width=.4, label="none-PGCs")

    pos = np.arange(len(germ_marker))
    lab = germ_marker
    print(pos)
    print(lab)
    ax.set_xticks(pos, lab, rotation=90)
    ax.set_yscale("log", base=10)
    ax.legend()
    fig.savefig(figout + "_germ_marker.svg")


def split_mtx_by_cid(barcodes, mtx, barcodes_cid_dic):
    dt_out = {}
    for cid in barcodes_cid_dic:
        idx = []
        for bar in barcodes_cid_dic[cid]:
            idx.append(barcodes.index(bar))
        dt_out[cid] = mtx.T[idx]
    return dt_out


def choose_top_genes(mtx_splited, features, tops=20):
    dt_out = {}
    for cid in mtx_splited:
        mtx  = mtx_splited[cid]
        mtx_sum = np.sum(mtx, axis=0)
        sort_list = list(zip(mtx_sum, range(len(features))))
        sort_list.sort(key=lambda x:x[0], reverse=True)
        idx = [e[1] for e in sort_list[:tops]]
        dt_out[cid] = [e[0] for e in np.asarray(features)[idx]]
    return dt_out


def get_top_genes(cluster_data_file, count_dire, out, tops=20):
    cluster_data = read_cluster_data(cluster_data_file)
    bar_f = os.path.join(count_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(count_dire, "features.tsv.gz")
    mtx_f = os.path.join(count_dire, "matrix.mtx.gz")
    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")
    filter_bool = filter_gene_have_cell(mtx)
    mtx = mtx[filter_bool]
    features = list(np.asarray(features)[filter_bool])
    mtx = norm_mtx(mtx)

    barcodes_cid_dic = {}
    for sid in cluster_data:
        for cid in cluster_data[sid]:
            if cid not in barcodes_cid_dic:
                barcodes_cid_dic[cid] = [e[0] for e in cluster_data[sid][cid]]
            else:
                barcodes_cid_dic[cid] += [e[0] for e in cluster_data[sid][cid]]
    mtx_splited = split_mtx_by_cid(barcodes, mtx, barcodes_cid_dic)
    top_genes = choose_top_genes(mtx_splited, features)
    print(top_genes)
    keys = list(top_genes.keys())
    keys.sort()
    fout = open(out + "_top_gene.txt", "w")
    for key in keys:
        print(">" + key, file=fout)
        for gene in top_genes[key]:
            print(gene, file=fout)
    fout.close()


def calculate_scatter_genes_data(mtx_split, features, clusterids, scatter_genes_align):
    dt_out = {}    
    features = [ele[0] for ele in features]
    gene_idx = []
    for gg in scatter_genes_align:
        gene_idx.append(features.index(gg))

    for cid in clusterids:
        dt_out[cid] = {}
        mtx = mtx_split[cid]
        mtx = mtx.T[gene_idx]
        idx = 0
        for row in mtx:
            dt_out[cid][scatter_genes_align[idx]] = [np.sum(row > 0) / len(row), np.average(row)]
            idx += 1
    return dt_out


def plot_the_scatter(scatter_gene_data, clusterids, scatter_genes_align, out, cluster_info):
    fig, axsd = plt.subplot_mosaic([["left", "rup"], ["left", "rbot"]], figsize=(13, 3), layout="constrained", width_ratios=[5, 1])
    ax = axsd["left"]
    xx = 0
    clusterids.reverse()
    for gg in scatter_genes_align:
        pct = []
        mexp = []
        x = []
        y = []
        yy = 0
        for cid in clusterids:
            x.append(xx)
            y.append(yy)
            yy += 1
            pct.append(scatter_gene_data[cid][gg][0] * 100)
            mexp.append(scatter_gene_data[cid][gg][1])
        xx += 1
        mexp = np.asarray(mexp) / max(mexp)
        ax.scatter(x, y, s=pct, c=mexp, cmap="Reds", edgecolors="#777777")
    ax.margins(.02, .1)
    ax.set_xticks(list(range(len(scatter_genes_align))), scatter_genes_align, rotation=90)
    y_lables = [cluster_info[c][0] for c in clusterids]
    ax.set_yticks(list(range(len(clusterids))), y_lables)

    axpct = axsd["rup"]
    title = axpct.set_title("Frection of cells\n in group")
    title.set_fontsize(8)
    x = np.linspace(1, 5, 5)
    y = np.ones(5)
    size = np.linspace(20, 100, 5)
    axpct.scatter(x, y, s=size, color="#aaaaaa")
    spines = axpct.spines
    for bod in spines:
        spines[bod].set_visible(False)
    axpct.set_xticks(x, ["20%", "40%", "60%", "80%", "100%"], rotation=90)
    axpct.set_yticks([], [])
    axpct.margins(.1, .1)
    axpct.set_ylim([.9, 1.1])

    axmexp = axsd["rbot"]
    title = axmexp.set_title("Relative Mean Expresion")
    title.set_fontsize(8)
    spines = axmexp.spines
    for bod in spines:
        spines[bod].set_visible(False)
    cmap = mpl.cm.Reds
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=axmexp, orientation='horizontal')
    fig.savefig(out + "_pct_mean_exp.svg")


def plot_gene_scatter(cluster_data_file, out, scatter_genes, cluster_info_file):
    cluster_data = read_cluster_data(cluster_data_file)
    cluster_info = read_cluster_info(cluster_info_file)
    bar_f = os.path.join(count_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(count_dire, "features.tsv.gz")
    mtx_f = os.path.join(count_dire, "matrix.mtx.gz")
    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")
    filter_bool = filter_gene_have_cell(mtx)
    mtx = mtx[filter_bool]
    features = list(np.asarray(features)[filter_bool])
    mtx = norm_mtx(mtx)

    barcodes_cid_dic = {}
    for sid in cluster_data:
        for cid in cluster_data[sid]:
            if cid not in barcodes_cid_dic:
                barcodes_cid_dic[cid] = [e[0] for e in cluster_data[sid][cid]]
            else:
                barcodes_cid_dic[cid] += [e[0] for e in cluster_data[sid][cid]]

    scatter_genes_dic = {}
    key = None
    for l in open(scatter_genes):
        l = l.rstrip()
        if l[0] == ">":
            key = l[1:]
            scatter_genes_dic[key] = []
        else:
            scatter_genes_dic[key].append(l)

    clusterids = list(scatter_genes_dic.keys())
    clusterids.sort()
    scatter_genes_align = []
    for cid in clusterids:
        scatter_genes_align += scatter_genes_dic[cid]

    mtx_split = split_mtx_by_cid(barcodes, mtx, barcodes_cid_dic)
    scatter_gene_data = calculate_scatter_genes_data(mtx_split, features, clusterids, scatter_genes_align)
    plot_the_scatter(scatter_gene_data, clusterids, scatter_genes_align, out, cluster_info)


def plot_germ_marker_violin(cluster_data_file, cluster_info_file, ctl, exp, out, count_dire, germ_marker):
    cluster_data = read_cluster_data(cluster_data_file)
    cluster_info = read_cluster_info(cluster_info_file)
    ctl = ctl.split(",")
    exp = exp.split(",")
    bar_f = os.path.join(count_dire, "barcodes.tsv.gz")
    fea_f = os.path.join(count_dire, "features.tsv.gz")
    mtx_f = os.path.join(count_dire, "matrix.mtx.gz")
    barcodes, features, mtx = read_10x_as_dense_mtx(bar_f, fea_f, mtx_f, "float32")
    filter_bool = filter_gene_have_cell(mtx)
    mtx = mtx[filter_bool]
    features = list(np.asarray(features)[filter_bool])
    mtx = norm_mtx(mtx)
    germ_marker = [l.rstrip() for l in open(germ_marker)]
    print(germ_marker)
    cluster_data = combin_sample(cluster_data, ctl, exp)
    ctl_bars_idx = []
    exp_bars_idx = []
    for cid in cluster_data["ctl"]:
        ctl_bars_idx += [barcodes.index(e[0]) for e in cluster_data["ctl"][cid]]
    for cid in cluster_data["exp"]:
        exp_bars_idx += [barcodes.index(e[0]) for e in cluster_data["exp"][cid]]

    features_gene = [e[0] for e in features]
    germ_idx = []
    for gg in germ_marker:
        try:
            germ_idx.append(features_gene.index(gg))
        except:
            print(gg, "not found")

    mtx = mtx[germ_idx]
    ctl_mtx = mtx.T[ctl_bars_idx]
    exp_mtx = mtx.T[exp_bars_idx]

    ctl_mtx = np.sum(ctl_mtx, axis=0)
    exp_mtx = np.sum(exp_mtx, axis=0)

    v1 = (ctl_mtx / (ctl_mtx + exp_mtx)) * 2
    v2 = (exp_mtx / (ctl_mtx + exp_mtx)) * 2
    print(v1)
    print(v2)
    v1_x = np.random.normal(1, .1, len(v1))
    v2_x = np.random.normal(2, .1, len(v2))
    fig, ax = plt.subplots()
    ax.violinplot([v1, v2], bw_method=.3)
    ax.scatter(v1_x, v1, color="#333333", s=6)
    ax.scatter(v2_x, v2, color="#333333", s=6)
    ax.set_xticks([1, 2], ["Control", "Mrbm24a"])
    fig.savefig(out + "_violinplot.svg")


if __name__ == "__main__":
    cluster_data_file, cluster_info_file, ctl, exp, out, marker_genes, count_dire, \
        headmap_genes, germ_marker, scatter_genes, violin_genes = getargs()
    plot_scatter(cluster_data_file, cluster_info_file, ctl, exp, out)
    plot_cloumn(cluster_data_file, cluster_info_file, ctl, exp, out)
    plot_marker_genes(cluster_data_file, cluster_info_file, ctl, exp, out, marker_genes, count_dire)
    plot_heatmap(cluster_data_file, cluster_info_file, ctl, exp, out, count_dire, headmap_genes)
    compare_marker_gene(cluster_data_file, cluster_info_file, ctl, exp, out, count_dire, germ_marker)
    get_top_genes(cluster_data_file, count_dire, out)
    plot_gene_scatter(cluster_data_file, out, scatter_genes, cluster_info_file)
    plot_germ_marker_violin(cluster_data_file, cluster_info_file, ctl, exp, out, count_dire, violin_genes)
