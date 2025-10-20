from collections import Counter


def count_nucleotides(sequence):
    nucleotides = Counter(sequence)
    alphabetically_sorted = sorted(nucleotides.items(), key=lambda x: x[0])
    return alphabetically_sorted


def load_input(filename):
    with open(filename, "r") as f:
        data = f.read().strip()
    return data


def parse_sequence_counts(sequence_counts):
    counts_values = [count for _nucleotide, count in sequence_counts]
    parsed_output = " ".join(map(str, counts_values))
    return parsed_output


if __name__ == "__main__":
    sequence = load_input("ex1.txt")
    counts = count_nucleotides(sequence)
    parsed_output = parse_sequence_counts(counts)
    print(parsed_output)
