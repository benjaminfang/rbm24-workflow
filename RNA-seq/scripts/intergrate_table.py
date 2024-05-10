#!

import argparse
import numpy as np


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", help="out file name")
    parser.add_argument("--gene_len", help="length of gene")
    parser.add_argument("--raw_count", help="original count mtx")
    parser.add_argument("--fpkm_count", help="FPKM normalize mtx")
    parser.add_argument("--tpm_count", help="TPM count mtx")
    parser.add_argument("--deseq_count", help="DESeq2 normalized mtx")
    parser.add_argument("--deseq_res", help="DESeq2 DE data")
    args = parser.parse_args()
    return args.out, args.gene_len, args.raw_count, args.fpkm_count, args.tpm_count, args.deseq_count, args.deseq_res


def read_gene_len(gene_len):
    data = {}
    fin = open(gene_len, "r")
    for line in fin:
        line = line.rstrip().split("\t")
        kk, val = line
        data[kk] = val
    fin.close()
    return data    


def read_countmtx_tsv(filename):
    data = []
    genes = []
    head = []

    fin = open(filename)
    head = fin.readline().rstrip().split("\t")[1:]
    for line in fin:
        line = line.rstrip().split("\t")
        genes.append(line[0])
        data.append(line[1:])
    fin.close()
    return head, genes, data


def read_data_deseq(filename):
    data = []
    genes = []
    head = []

    fin = open(filename)
    head = [e.strip('"') for e in fin.readline().rstrip().split(",")[1:]]
    for line in fin:
        line = line.rstrip().split(",")
        genes.append(line[0].strip('"'))
        data.append(line[1:])
    return head, genes, data


def check_genes_consistence(raw_genes, fpkm_genes, tpm_genes, deseq_genes, res_genes):
    for ele in zip(raw_genes, fpkm_genes, tpm_genes, deseq_genes, res_genes):
        if len(set(ele)) != 1:
            print("gene not consistent")
            return
    return


def check_head_consistence(raw_head, fpkm_head, tpm_head, deseq_head):
    for ele in zip(raw_head, fpkm_head, tpm_head, deseq_head):
        if len(set(ele)) != 1:
            print("head not consistence")
            return
    return


def combine_data(raw_data, fpkm_data, tpm_data, deseq_data, res_data):
    dt = []
    row_n = len(raw_data)
    fn = len(raw_data[0])

    for i in range(row_n):
        ll = []
        for j in range(fn):
            ll += [raw_data[i][j], fpkm_data[i][j], tpm_data[i][j], deseq_data[i][j]]
        ll += res_data[i]
        dt.append(ll)
    
    return dt


def output_dt(raw_head, res_head, raw_genes, genelen, data, out):
    fout = open(out, "w")
    head = ["geneid", "genelen"]
    for ele in raw_head:
        head += [ele, "fpkm", "tpm", "deseq_norm"]
    head += res_head
    print("\t".join(head), file=fout)

    data_len = len(data)
    for i in range(data_len):
        print("\t".join([raw_genes[i], genelen[i]] + data[i]), file=fout)

    fout.close()


def main():
    out, gene_len, raw_count, fpkm_count, tpm_count, deseq_count, deseq_res = getargs()
    
    genelen_dic = read_gene_len(gene_len)
    raw_head, raw_genes, raw_data = read_countmtx_tsv(raw_count)
    gene_len = []
    for gg  in raw_genes:
        gene_len.append(genelen_dic[gg])

    fpkm_head, fpkm_genes, fpkm_data = read_countmtx_tsv(fpkm_count)
    tpm_head, tpm_genes, tpm_data = read_countmtx_tsv(tpm_count)
    deseq_head, deseq_genes, deseq_data = read_data_deseq(deseq_count)
    res_head, res_genes, res_data = read_data_deseq(deseq_res)

    check_genes_consistence(raw_genes, fpkm_genes, tpm_genes, deseq_genes, res_genes)
    check_head_consistence(raw_head, fpkm_head, tpm_head, deseq_head)

    data = combine_data(raw_data, fpkm_data, tpm_data, deseq_data, res_data)
    output_dt(raw_head, res_head, raw_genes, gene_len, data, out)

    return 0


if __name__ == "__main__":
    main()