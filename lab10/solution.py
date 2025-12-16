import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def load_fasta(filename):
    sequences = SeqIO.parse(filename, "fasta")
    return [str(sequence.seq) for sequence in sequences]


def sample_subsequences(sequence, length, coverage):
    idx_max = len(sequence) - length
    n_subsequences = int(len(sequence) * coverage / length)

    indices = random.choices(range(idx_max), k=n_subsequences)
    return [sequence[idx:idx + length] for idx in indices]


def merge_subsequences(sequences_filename):
    """
    Function to merge subsequences into a single sequence.
    Implemented using iteration as python does not support tail recursion.
    Iterative approach avoids RecursionError.
    """
    sequences = load_fasta(sequences_filename)
    print(f"Constructing contigs using {len(sequences)} fragments...")
    best_left, best_right = None, None
    while True:
        max_overlap = 0
        for left in sequences:
            for right in sequences:
                if left != right:
                    overlap = 0
                    length = min(len(left), len(right))
                    for i in range(1, length):
                        if right[:i] == left[-i:]:
                            overlap = i
                    if overlap > max_overlap:
                        best_left, best_right = left, right
                        max_overlap = overlap
        if max_overlap != 0:
            new_sequence = best_left + best_right[max_overlap:]
            sequences.remove(best_right)
            sequences.remove(best_left)
            sequences.append(new_sequence)
        else:
            break
    return sequences


def save_subsequences(fragments, filename):
    records = [SeqRecord(Seq(frag), id=f"read_{i}", description="")
               for i, frag in enumerate(fragments, 1)]
    SeqIO.write(records, filename, "fasta")


def print_summary(contigs, subsequence=None):
    print(f"Found {len(contigs)} contig(s)")
    for n, contig in enumerate(contigs, start=1):
        print(f"{n}. | length: {len(contig)} | contig: {contig}")
        if subsequence is not None:
            print(f"Is found contig a part of the sequence? - {contig in subsequence}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_filename",
                        help="Input sequence filename. Expected format is fasta. Loads the first sequence in the file."
                        )
    parser.add_argument("--contig_only", help="Treats input file as sequences to merge.", action="store_true")
    parser.add_argument("--output", help="Output sequence filename.", default="out.fa")
    parser.add_argument("--bp", help="Number of bp's to take from analysed sequence.", type=int)
    parser.add_argument("--coverage",
                        help="Redundancy factor of the sampling process", type=int, default=5
                        )
    parser.add_argument("--length", help="Length of the sequence", type=int, default=200)

    args = parser.parse_args()
    if args.contig_only:
        contigs = merge_subsequences(args.input_filename)
        print_summary(contigs)
    else:
        sequences = load_fasta(args.input_filename)
        first_sequence = sequences[0]
        subsequence = first_sequence[:args.bp]
        fragments = sample_subsequences(subsequence, args.length, args.coverage)

        print(f"Saving sampled fragments to {args.output}")
        save_subsequences(fragments, args.output)
        contigs = merge_subsequences(args.output)
        print_summary(contigs, subsequence)







