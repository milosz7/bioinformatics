def load_input(filename):
    with open(filename) as f:
        [dna, subseq] = f.read().split("\n")
        return dna, subseq


def find_subsequences(dna, seq):
    indices = []
    start = 0
    while True:
        start = dna.find(seq, start)
        if start == -1:
            break
        indices.append(start + 1)
        start += 1
    return indices


if __name__ == "__main__":
    dna, subseq = load_input("ex7.txt")
    subsequences = find_subsequences(dna, subseq)
    subsequences = [str(idx) for idx in subsequences]
    print(" ".join(subsequences))
