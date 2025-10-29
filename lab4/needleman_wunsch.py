from collections import defaultdict
from Bio import Entrez, SeqIO
from Bio.Align import substitution_matrices
import numpy as np


class NeedlemanWunsch:
    def __init__(self, matrix, gap_penalty):
        self.matrix = matrix
        self.gap_penalty = gap_penalty
        self.paths = []

    def cost(self, elem_0, elem_1):
        return self.matrix[elem_0][elem_1]

    def align(self, v, w):
        n, m = len(v), len(w)
        S = np.zeros((n + 1, m + 1), dtype=int)
        backtrack = defaultdict(list)
        for j in range(1, m + 1):
            backtrack[(0, j)].append((0, j - 1))
        for i in range(1, n + 1):
            backtrack[(i, 0)].append((i - 1, 0))

        S[0, :] = np.arange(m + 1) * self.gap_penalty
        S[:, 0] = np.arange(n + 1) * self.gap_penalty
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diag_cost = S[i - 1, j - 1] + self.cost(v[i - 1], w[j - 1])
                col_cost = S[i - 1, j] + self.gap_penalty
                row_cost = S[i, j - 1] + self.gap_penalty
                max_value = max(diag_cost, row_cost, col_cost)
                S[i, j] = max_value
                if max_value == diag_cost:
                    backtrack[i, j].append((i - 1, j - 1))
                if max_value == row_cost:
                    backtrack[i, j].append((i, j - 1))
                if max_value == col_cost:
                    backtrack[i, j].append((i - 1, j))

        self._backtrack(backtrack, (n, m), (0, 0), [(n, m)])
        aligned_sequences = self._pad_sequences(v, w)
        return aligned_sequences

    def _pad_sequences(self, v, w):
        alignments = []
        for path in self.paths:
            v_copy, w_copy = v, w
            path = path[::-1]
            for n, (curr, prev) in enumerate(zip(path, path[1:])):
                if curr[0] == prev[0] and curr[1] < prev[1]:
                    v_copy = v_copy[:n] + "-" + v_copy[n:]
                if curr[0] < prev[0] and curr[1] == prev[1]:
                    w_copy = w_copy[:n] + "-" + w_copy[n:]
            alignments.append((v_copy, w_copy))
        self.paths = []
        return alignments

    def _backtrack(self, graph, current, target, path):
        if current == target:
            self.paths.append(path.copy())
            return

        for neighbor in graph[current]:
            if neighbor not in path:
                path.append(neighbor)
                self._backtrack(graph, neighbor, target, path)
                path.pop()


def load_gene(protein_id):
    with Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")
        return str(record.seq)


if __name__ == "__main__":
    Entrez.email = input("Please enter your email address: ")
    human_hemoglobin = load_gene(40886941)
    rat_hemoglobin = load_gene(34849618)
    blosum62 = substitution_matrices.load("BLOSUM62")

    nw = NeedlemanWunsch(blosum62, -7)
    alignments = nw.align(human_hemoglobin, rat_hemoglobin)
    for n, aligned_sequences in enumerate(alignments):
        print(f"==== Alignment {n} ====")
        print("\n".join(aligned_sequences))
