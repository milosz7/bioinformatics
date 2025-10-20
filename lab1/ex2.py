def load_input(filename):
    with open(filename, "r") as f:
        data = f.read().strip()
    return data


def reverse_complement(sequence):
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    reverse_sequence = sequence[::-1]
    return "".join(complement_map[nucleotide] for nucleotide in reverse_sequence)


if __name__ == "__main__":
    sequence = load_input("ex2.txt")
    reverse_complement_sequence = reverse_complement(sequence)
    print(reverse_complement_sequence)
