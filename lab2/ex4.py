from Bio import SeqIO


def load_fasta_with_index(filename):
    sequences = SeqIO.parse(filename, 'fasta')
    return [(sequence.id, str(sequence.seq)) for sequence in sequences]


def calculate_gc_content(sequence):
    gc_count = sum(nuc == 'G' or nuc == 'C' for nuc in sequence)
    return gc_count * 100 / len(sequence)


def find_highest_gc_count(sequences):
    gc_contents = [(seq_id, calculate_gc_content(seq)) for seq_id, seq in sequences]
    return max(gc_contents, key=lambda gc_content: gc_content[1])


if __name__ == "__main__":
    sequences = load_fasta_with_index("ex4.txt")
    seq_id, gc_content = find_highest_gc_count(sequences)
    print(f"{seq_id}\n{gc_content}")
