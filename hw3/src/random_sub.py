import sys
import os
import random

if __name__ == "__main__":
    input_path = sys.argv[1]
    map_path = sys.argv[2]
    output_path = sys.argv[3]
    D = dict()
    with open(map_path, "r", encoding='cp950') as fr:
        for line in fr.readlines():
            chars = line.strip().split()
            chars = list(filter(lambda x : len(x)>0, chars))
            root = chars[0] # Big5
            key = chars[0]
            values = []
            for c in chars[1:]:
                values += c.split('/')
            D[key] = values
    with open(input_path, "r", encoding='cp950') as fr:
        with open(output_path, "w", encoding='cp950') as fw:
            for in_line in fr:
                chars = in_line.strip().split()
                char_len = len(chars)
                n_sub = int(char_len * 0.3)
                sub_pos = random.sample(range(char_len), n_sub)
                for i in sub_pos:
                    if not chars[i] in D:
                        continue
                    c_sub = random.choice(D[chars[i]])[0] # Take the first pinin
                    chars[i] = c_sub
                out_line = ' '.join(chars)
                fw.write(out_line+'\n')
