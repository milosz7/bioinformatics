from Bio import SeqIO
from collections import defaultdict


def load_fasta_with_index(filename):
    sequences = SeqIO.parse(filename, 'fasta')
    return [(sequence.id, str(sequence.seq)) for sequence in sequences]


def map_sequences(sequences):
    triplets_map = defaultdict(list)
    for sequence_id, sequence in sequences:
        triplets_map[sequence[:3]].append((sequence_id, sequence))
    return triplets_map


def find_adjacent_sequences(current_pair, triplets_map):
    found_pairs = set()
    current_id, current_seq = current_pair
    for sequence_id, sequence in triplets_map[current_seq[-3:]]:
        if sequence_id != current_id:
            found_pairs.add((current_id, sequence_id))
    return found_pairs


def create_adjacency_list(sequences):
    start_triplets_map = map_sequences(sequences)
    all_pairs = set()
    for sequence in sequences:
        all_pairs |= find_adjacent_sequences(sequence, start_triplets_map)
    return all_pairs


def format_output(all_pairs):
    pairs_fmt = [" ".join(pair) for pair in all_pairs]
    return "\n".join(pairs_fmt)


if __name__ == '__main__':
    sequences = load_fasta_with_index("ex5.txt")
    all_pairs = create_adjacency_list(sequences)
    print(format_output(all_pairs))

