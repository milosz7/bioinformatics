from Bio import SeqIO


def load_fasta(filename):
    sequences = SeqIO.parse(filename, 'fasta')
    return [str(sequence.seq) for sequence in sequences]


def is_mutation_transition(nucleotide1, nucleotide2):
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}

    if nucleotide1 in purines and nucleotide2 in purines:
        return True
    if nucleotide1 in pyrimidines and nucleotide2 in pyrimidines:
        return True
    return False



def transition_transversion_ratio(seq1, seq2):
    mutations, transitions = 0, 0
    for nuc1, nuc2 in zip(seq1, seq2):
        if nuc1 != nuc2:
            mutations += 1
            transitions += int(is_mutation_transition(nuc1, nuc2))
    transversions = mutations - transitions
    return transitions / transversions


if __name__ == "__main__":
    [seq1, seq2] = load_fasta("ex3.txt")
    print(transition_transversion_ratio(seq1, seq2))
