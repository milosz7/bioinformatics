from Bio import SeqIO
from Bio.Align import substitution_matrices, PairwiseAligner


def load_fasta(filename):
    sequences = SeqIO.parse(filename, "fasta")
    return [str(sequence.seq) for sequence in sequences]


def init_aligner(subst_mat, open_gap_score, extend_gap_score, mode = "global"):
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = subst_mat
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    return aligner


def get_best_subsequences(seq1, seq2, aligner):
    alignments = aligner.align(seq1, seq2)
    alignment = next(alignments)
    score = alignment.score
    aligned_seq1 = alignment.aligned[0]
    aligned_seq2 = alignment.aligned[1]

    if len(aligned_seq1) > 0:
        start1, end1 = aligned_seq1[0][0], aligned_seq1[-1][1]
        start2, end2 = aligned_seq2[0][0], aligned_seq2[-1][1]

        fragment1 = seq1[start1:end1]
        fragment2 = seq2[start2:end2]
        return score, fragment1, fragment2
    raise ValueError("No aligned subsequences found.")

if __name__ == "__main__":
    [seq1, seq2] = load_fasta("ex3.txt")
    subst_mat = substitution_matrices.load("PAM250")
    open_gap_penalty = -5
    extend_gap_penalty = -5
    aligner_mode = "local"
    aligner = init_aligner(subst_mat, open_gap_penalty, extend_gap_penalty, aligner_mode)
    score, *sequences = get_best_subsequences(seq1, seq2, aligner)
    print(score)
    print("\n".join(sequences))
