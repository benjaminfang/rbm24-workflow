import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats_file", required=True, 
                        help="barcode statistic file, generated in data intergrate step.")
    parser.add_argument("--out", default="out", help="out file name")
    parser.add_argument("sample_prefix", nargs=4, help="name of sample.")
    args = parser.parse_args()
    return args.stats_file, args.sample_prefix, args.out


def plot_data(data, sample_pre, out):
    dt_index = np.asarray(list([e.split("_")[0] for e in data.index]))
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8), layout="constrained")
    count_sum_all = []
    gene_num_all = []
    top20_pct_all = []
    mt_pct_all = []
    for sam in sample_pre:
        print(">>>", sam)
        dt_sam = data[dt_index == sam]
        count_sum = np.int32(np.round(dt_sam["count_sum"].to_numpy()))
        gene_num = np.int32(np.round(dt_sam["gene_num"].to_numpy()))
        top20_pct = dt_sam["top20_pct"].to_numpy()
        mt_pct = dt_sam["mt_pct"].to_numpy()
        print("count_sum", np.quantile(count_sum, [0, .25, .5, .75, 1]))
        print("gene_sum", np.quantile(gene_num, [0, .25, .5, .75, 1]))
        print("top20_pct", np.quantile(top20_pct, [0, .25, .5, .75, 1]))
        print("mt_pct", np.quantile(mt_pct, [0, .25, .5, .75, 1]))
        count_sum_all.append(count_sum)
        gene_num_all.append(gene_num)
        top20_pct_all.append(top20_pct)
        mt_pct_all.append(mt_pct)
    for idx, sam in enumerate(sample_pre):
        axs[0, 0].hist(count_sum_all[idx], bins=100, density=True, cumulative=True, histtype="step", label=sam)
        axs[0, 1].hist(gene_num_all[idx], bins=100, density=True, cumulative=True, histtype="step", label=sam)
        axs[1, 0].hist(top20_pct_all[idx], bins=100, density=True, cumulative=True, histtype="step", label=sam)
        axs[1, 1].hist(mt_pct_all[idx], bins=100, density=True, cumulative=True, histtype="step", label=sam)
    axs[0, 0].legend()
    axs[0, 0].set_title("sum of count")
    axs[0, 1].legend()
    axs[0, 1].set_title("sum of gene number")
    axs[1, 0].legend()
    axs[1, 0].set_title("percent of top20 genes")
    axs[1, 1].legend()
    axs[1, 1].set_title("percent of mt genes")
    fig.savefig(out + ".svg")


def main():
    stats_file, sample_pre, out = getargs()
    data = pd.read_csv(stats_file, sep="\t", header=0, index_col=0)
    plot_data(data, sample_pre, out)


if __name__ == "__main__":
    main()
