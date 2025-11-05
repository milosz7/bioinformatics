from Bio import SeqIO
from Bio.Align import substitution_matrices, PairwiseAligner


def load_fasta(filename):
    sequences = SeqIO.parse(filename, "fasta")
    return [str(sequence.seq) for sequence in sequences]


def init_aligner(subst_mat, open_gap_score, extend_gap_score, mode = "global"):
    aligner = PairwiseAligner(scoring="blastn")
    aligner.mode = mode
    aligner.substitution_matrix = subst_mat
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    return aligner


def global_alignment_score(seq1, seq2, aligner):
    alignments = aligner.align(seq1, seq2)
    return alignments[0].score


if __name__ == "__main__":
    [seq1, seq2] = load_fasta("ex2.txt")
    subst_mat = substitution_matrices.load("BLOSUM62")
    open_gap_penalty = -5
    extend_gap_penalty = 0
    aligner = init_aligner(subst_mat, open_gap_penalty, extend_gap_penalty)
    score = global_alignment_score(seq1, seq2, aligner)
    print(f"Match score: {score}")
