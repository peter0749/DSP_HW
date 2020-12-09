import sys

if __name__ == '__main__':
    output_file = sys.argv[1]
    reference_file = sys.argv[2]
    n_total = 0
    n_correct = 0
    with open(output_file, encoding='cp950') as of:
        with open(reference_file, encoding='cp950') as rf:
            for line_o, line_r in zip(of, rf):
                out = line_o.strip()[4:-5].split()
                gt = line_r.strip()[4:-5].split()
                assert (len(gt) == len(out))
                for i, j in zip(out, gt):
                    n_correct += 1 if (i==j) else 0
                n_total += len(gt)
    print("Accuracy: %.4f (#wrong: %d, #total: %d)"%(n_correct / n_total, (n_total-n_correct), n_total))
