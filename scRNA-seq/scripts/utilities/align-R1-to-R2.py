import sys
import gzip

r1_f = sys.argv[1]
r2_f = sys.argv[2]

r1_f_out = sys.argv[3]

fin_r1 = gzip.open(r1_f, "r")
fin_r2 = gzip.open(r2_f, "r")
fout_r1 = gzip.open(r1_f_out, "w")

r2_tmp = []
r1_tmp = []

while True:
    r2_tmp = []
    for i in range(10000):
        line = fin_r2.readline()
        if line:
            r2_tmp.append(line.decode().split()[0])
        else:
            r2_tmp.append("")
        fin_r2.readline()
        fin_r2.readline()
        fin_r2.readline()

    if "" in r2_tmp:
        break

    eat_len = 0
    while True:
        r1_tmp = []
        if eat_len == 10000:
            break
        r1_tmp.append(fin_r1.readline())
        r1_tmp.append(fin_r1.readline())
        r1_tmp.append(fin_r1.readline())
        r1_tmp.append(fin_r1.readline())
        if r1_tmp[0].decode().split()[0] in r2_tmp:
            eat_len += 1
            for l in r1_tmp:
                fout_r1.write(l)

while True:
    r1_tmp = []

    r1_tmp.append(fin_r1.readline())
    r1_tmp.append(fin_r1.readline())
    r1_tmp.append(fin_r1.readline())
    r1_tmp.append(fin_r1.readline())
    if b"" in r1_tmp:
        break

    if r1_tmp[0].decode().split()[0] in r2_tmp:
        for l in r1_tmp:
            fout_r1.write(l)


fout_r1.close()
fin_r1.close()
fin_r2.close()

