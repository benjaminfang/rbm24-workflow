#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", help="out file name")
    parser.add_argument("--norm", choices={"fpkm", "tpm", "deseq_norm"}, help="normalization method")
    parser.add_argument("data_file", nargs=3, help="unioned data tables")
    parser.add_argument("exp_name", nargs=3, help="experimental time point name")

    args = parser.parse_args()
    return args.out, args.norm, args.data_file, args.exp_name


def read_data(filename):
    head = []
    genes = []
    data = []
    fin = open(filename)
    head = fin.readline().rstrip().split("\t")[1:]
    for line in fin:
        line = line.rstrip().split("\t")
        genes.append(line[0])
        tmp = []
        for ele in line[1:]:
            if ele == "NA":
                tmp.append(np.nan)
            else:
                tmp.append(float(ele))
        data.append(np.asarray(tmp, dtype=np.float32))
    data = np.asarray(data)
    data_sum = np.sum(data, axis=0)
    row_idx = []
    i = 0
    for e in head:
        if e.startswith("M_") or e.startswith("Sibling_"):
            row_idx.append(i)
        i += 1
    print(data_sum[row_idx])
    return head, genes, data


def pick_data(head, data, norm_method):
    exp_ord = []
    ctl_ord = []
    idx = 0
    for field in head:
        if field.startswith("M_"):
            exp_ord.append(idx)
            idx += 1    
        elif field.startswith("Sibling_"):
            ctl_ord.append(idx)
            idx += 1
    data_idx = []
    idx = 0
    for field in head:
        if field == norm_method:
            data_idx.append(idx)
        idx += 1
    
    mean_padj = []
    idx = 0
    for field in head:
        if field == "baseMean" or field == "log2FoldChange" or field == "padj" or field == "pvalue":
            mean_padj.append(idx)
        idx += 1


    exp_idx = np.asarray(data_idx)[exp_ord]
    ctl_idx = np.asarray(data_idx)[ctl_ord]
    
    print(exp_idx)
    print(ctl_idx)
    print(mean_padj)
    print("--------")
    
    data = data.T
    exp_data = data[exp_idx]
    ctl_data = data[ctl_idx]
    mean_padj_data = data[mean_padj]

    return exp_data.T, ctl_data.T, mean_padj_data.T


def plot_ma(data, ax, ax_tilte):
    l2fc_max = 5
    data = data.T
    meanexp = data[0]
    l2fc = data[1]
    padj = data[2]
    """
    keep_bool = []
    for idx, ele in enumerate(padj):
        if not np.isnan(ele):
            keep_bool.append(idx)
    meanexp = meanexp[keep_bool]
    l2fc = l2fc[keep_bool]
    padj = padj[keep_bool]
    """

    out_mean = []
    out_l2fc = []
    out_color = []
    out_alpha = []
    out_edc = []
    out_edw = []

    in_mean = []
    in_l2fc = []
    color = []
    alpha = []
    edgecolor = []
    edgewith = []
    for me, fc, pj in zip(meanexp, l2fc, padj):
        if abs(fc) > l2fc_max:
            out_mean.append(me)
            if fc > 0:
                out_l2fc.append(l2fc_max)
            else:
                out_l2fc.append(-l2fc_max)
            if np.isnan(pj):
                out_color.append("white")
                out_alpha.append(.5)
                out_edc.append("black")
                out_edw.append(.1)

            else:
                if pj >= .01:
                    out_color.append("gray")
                    out_alpha.append(.5)
                    out_edc.append("gray")
                    out_edw.append(0)
                else:
                    out_color.append("red")
                    out_alpha.append(1)
                    out_edc.append("gray")
                    out_edw.append(0)
        else:
            in_mean.append(me)
            in_l2fc.append(fc)
            if np.isnan(pj):
                color.append("white")
                alpha.append(.5)
                edgecolor.append("black")
                edgewith.append(.1)
            else:
                if pj < .01 and abs(fc) > 1:
                    color.append("red")
                    alpha.append(1)
                    edgecolor.append("gray")
                    edgewith.append(0)
                else:
                    color.append("gray")
                    alpha.append(.5)
                    edgecolor.append("gray")
                    edgewith.append(0)
    if len(out_mean) > 0:
        ax.scatter(out_mean, out_l2fc, color=out_color, alpha=out_alpha, s=5, linewidths=out_edw, edgecolor=out_edc, marker="^")
    ax.scatter(in_mean, in_l2fc, color=color, alpha=alpha, linewidths=edgewith, edgecolor=edgecolor, marker="o", s=5)
    ax.set_title(ax_tilte)
    ax.set_xscale("log", base=10)
    ax.set_ylim([-5.5, 5.5])



def plot_ma2(geneid, data, ax, ax_tilte):
    """
    Same as function plot_ma except the nan value will shown as a not significant
    point.
    """
    print(ax_tilte)

    germ_marker = ["nanos3", "ddx4", "dnd1", "tdrd7a", "ca15b", "celf1", "rgs14a", "dazl", "hook2", "tdrd6", "gra", "h1m"]
    l2fc_max = 5
    data = data.T
    meanexp = data[0]
    l2fc = data[1]
    pval = data[2]


    out_ss = []
    out_mean = []
    out_l2fc = []
    out_color = []
    out_alpha = []
    in_ss = []
    in_mean = []
    in_l2fc = []
    color = []
    alpha = []
    idx = 0
    t_x = []
    t_y = []
    tt = []
    for me, fc, pv in zip(meanexp, l2fc, pval):
        if fc > l2fc_max:
            out_mean.append(me)
            out_l2fc.append(l2fc_max)
            if np.isnan(pv):
                out_color.append("gray")
                out_alpha.append(.5)
                out_ss.append(5)
            else:
                if pv >= .01:
                    out_color.append("gray")
                    out_alpha.append(.5)
                    out_ss.append(5)
                elif pv < .01 and (geneid[idx] not in germ_marker):
                    out_color.append("red")
                    out_alpha.append(1)
                    out_ss.append(5)
                else:
                    t_x.append(me)
                    t_y.append(l2fc_max)
                    tt.append(geneid[idx])
                    out_color.append("#00ff01")
                    out_alpha.append(1)
                    out_ss.append(10)
        elif fc < -l2fc_max:
            out_mean.append(me)
            out_l2fc.append(-l2fc_max)
            if np.isnan(pv):
                out_color.append("gray")
                out_alpha.append(.5)
                out_ss.append(5)
            else:
                if pv >= .01:
                    out_color.append("gray")
                    out_alpha.append(.5)
                    out_ss.append(5)
                elif pv < .01 and (geneid[idx] not in germ_marker):
                    out_color.append("blue")
                    out_alpha.append(1)
                    out_ss.append(5)
                else:
                    t_x.append(me)
                    t_y.append(-l2fc_max)
                    tt.append(geneid[idx])
                    out_color.append("#00ff01")
                    out_alpha.append(1)
                    out_ss.append(10)
        else:
            in_mean.append(me)
            in_l2fc.append(fc)
            if np.isnan(pv):
                color.append("gray")
                alpha.append(.5)
                in_ss.append(5)
            else:
                if pv < .01 and fc > 1 and (geneid[idx] not in germ_marker):
                    color.append("red")
                    alpha.append(1)
                    in_ss.append(5)
                elif pv < .01 and fc < -1 and (geneid[idx] not in germ_marker):
                    color.append("blue")
                    alpha.append(1)
                    in_ss.append(5)
                elif pv < .01 and abs(fc) > 1 and (geneid[idx] in germ_marker):
                    t_x.append(me)
                    t_y.append(fc)
                    tt.append(geneid[idx])
                    color.append("#00ff01")
                    alpha.append(1)
                    in_ss.append(10)
                else:
                    color.append("gray")
                    alpha.append(.5)
                    in_ss.append(5)
        idx += 1
    
    if len(out_mean) > 0:
        ax.scatter(out_mean, out_l2fc, color=out_color, alpha=out_alpha, linewidths=0, s=out_ss, marker="^")
    ax.scatter(in_mean, in_l2fc, color=color, alpha=alpha, linewidths=0, marker="o", s=in_ss)
    for x, y, te in zip(t_x, t_y, tt):
        ax.text(x, y, te)
    ax.set_title(ax_tilte)
    ax.set_xscale("log", base=10)
    ax.set_ylim([-5.5, 5.5])
    ax.set_xlabel("mean expression ($log_{10}$)")
    ax.set_ylabel("log fold change ($log_{2}$)")


def plot_marker_exp(all_data, ax, marker, exp_name):
    exp_mean = []
    ctl_mean = []
    exp_std = []
    ctl_std = []
    for point in all_data:
        genes, exp, ctl, mean_padj = point
        marker_idx = genes.index(marker)
        exp = exp[marker_idx]
        ctl = ctl[marker_idx]
        exp_mean.append(np.mean(exp))
        exp_std.append(np.std(exp))
        ctl_mean.append(np.mean(ctl))
        ctl_std.append(np.std(ctl))
    
    ax.plot(exp_name, exp_mean, marker="o", label="exp")
    ax.plot(exp_name, ctl_mean, marker="x", label="ctl")
    ax.errorbar(exp_name, exp_mean, yerr=exp_std, linestyle="")
    ax.errorbar(exp_name, ctl_mean, yerr=ctl_std, linestyle="")
    ax.set_title(marker)
    ax.legend()


def plot_marker_exp_column(exp_data, germ_marker, ax, exp_name, norm_method, fout_col):
    genes, exp, ctl, mean_padj = exp_data
    marker_idx = []
    for gg in germ_marker:
        marker_idx.append(genes.index(gg))

    exp_item = exp[marker_idx]
    ctl_item = ctl[marker_idx]
    print(exp_name, file=fout_col)
    print("\t".join(["geneid"] + ["exp"] * exp_item.shape[1] + ["ctl"] * ctl_item.shape[1]), file=fout_col)
    dt = np.hstack((exp_item, ctl_item))

    idx = 0
    for row in dt:
        print("\t".join([germ_marker[idx]] + [format(e, ".5f") for e in row]), file=fout_col)
        idx += 1

    exp_item_avg = np.average(exp_item, axis=1)
    ctl_item_avg = np.average(ctl_item, axis=1)
    exp_item_sd = np.std(exp_item, axis=1)
    ctl_item_sd = np.std(ctl_item, axis=1)
    x = np.arange(len(exp_item))
    padj = mean_padj[marker_idx][:, 2]
    fc = mean_padj[marker_idx][:, 1]
    avg_max = np.max(np.vstack((exp_item_avg, ctl_item_avg)), axis=0)
    ax.bar(x - .15, ctl_item_avg, yerr=ctl_item_sd, width=.3, label="Ctl")
    ax.bar(x + .15, exp_item_avg, yerr=exp_item_sd, width=.3, label="Exp")
    for tx, ty, pp, ff in zip(x, avg_max, padj, fc):
        if pp < .01 and abs(ff) > 1: 
            ax.text(tx, ty + 5, "*")
        
    ax.legend()
    ax.set_xticks(x, germ_marker, rotation=90)
    ax.set_title(exp_name + f"({norm_method.upper()})")
    
    return


if __name__ == "__main__":
    print(">>>")
    out, norm_method, data_file, exp_name = getargs()
    print(exp_name)
    #mir202 not has reads mapped.
    germ_marker = ["nanos3", "ddx4", "dnd1", "tdrd7a", "ca15b", "celf1", "rgs14a", "dazl", "hook2", "tdrd6", "gra", "h1m"]
    all_data = []
    for ff in data_file:
        print(ff)
        head, genes, data = read_data(ff)
        exp, ctl, mean_padj = pick_data(head, data, norm_method)
        all_data.append([genes, exp, ctl, mean_padj])
    
    
    fig, axs = plt.subplots(nrows=1, ncols=len(exp_name), figsize=(5 * len(exp_name), 4.8), layout="constrained")
    for i, point in enumerate(exp_name):
        plot_ma2(all_data[i][0], all_data[i][3], axs[i], point)
    fig.savefig(out + "_ma.pdf")

    """
    fig, axs = plt.subplots(nrows=1, ncols=len(exp), figsize=(5 * len(germ_marker), 4.8), layout="constrained")
    for i, marker in enumerate(germ_marker):
        plot_marker_exp(all_data, axs[i], marker, exp_name)
    fig.savefig(out + "_germmarker.pdf")
    """
    """
    fout_col = open(out + "-column-data.tsv", "w")
    fig, axs = plt.subplots(nrows=1, ncols=len(exp_name), figsize=(5 * len(exp_name), 4.8), layout="constrained")
    for i, exp_data in enumerate(all_data):
        plot_marker_exp_column(exp_data, germ_marker, axs[i], exp_name[i], norm_method, fout_col)
    fig.savefig(out + "-germmarker-bar-" + norm_method + ".pdf")
    fout_col.close()
    """
    

