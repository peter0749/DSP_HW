import sys
import os

if __name__ == "__main__":
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    D = dict()
    with open(input_path, "r", encoding='cp950') as fr:
        for line in fr.readlines():
            chars = line.strip().split()
            chars = list(filter(lambda x : len(x)>0, chars))
            root = chars[0] # Big5
            for poly in chars[1:]: # ZhuYin
                cc = list(filter(lambda x : len(x)>0, poly.split('/')))
                for c in cc:
                    if not c[0] in D:
                        D[c[0]] = set()
                    D[c[0]].add(root)
            D[root] = root
    with open(output_path, "w", encoding='cp950') as fw:
        for fc, tcs in D.items():
            fw.write(fc)
            for c in tcs:
                fw.write(" " + c)
            fw.write("\n")

