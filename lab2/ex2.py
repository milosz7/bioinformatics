from Bio import SeqIO


def load_fasta(filename):
    sequences = SeqIO.parse(filename, 'fasta')
    return [str(sequence.seq) for sequence in sequences]


def longest_common_subsequence(sequences):
    shortest_seq = min(sequences, key=len)
    for i in range(len(shortest_seq), 0, -1):
        subsequences = []
        for sequence in sequences:
            subsequences.append({sequence[j:j+i] for j in range(len(sequence) - i)})
        intersection = subsequences[0]
        for subsequence in subsequences[1:]:
            intersection &= subsequence
        if intersection:
            return max(intersection, key=len)
    return ''


if __name__ == "__main__":
    data = load_fasta("ex2.txt")
    print(longest_common_subsequence(data))