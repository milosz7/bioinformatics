def load_input(filename):
    with open(filename, "r") as f:
        data = f.read().strip()
    return data

def sequence_hamming_distance(seq1, seq2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))


if __name__ == '__main__':
    [seq1, seq2] = load_input("ex1.txt").split("\n")
    print(sequence_hamming_distance(seq1, seq2))
