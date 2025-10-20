from Bio import SeqIO


def load_fasta(filename):
    sequences = SeqIO.parse(filename, 'fasta')
    return [str(sequence.seq) for sequence in sequences]


def dna_to_protein(dna, offset):
    gencode = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGC": "C",
        "TGT": "C",
        "TGA": "*",
        "TGG": "W",
    }
    output_sequence = ""
    for i in range(offset, len(dna), 3):
        code = dna[i : i + 3]
        next_aminoacid = gencode.get(code, "-")
        if next_aminoacid == "*":
            break
        output_sequence += next_aminoacid

    return output_sequence

def reverse_complement(dna):
    complement = str.maketrans("ATCG", "TAGC")
    return dna.translate(complement)[::-1]


def find_start_codons(dna):
    start_codon = "ATG"
    indices = []
    start = 0
    while True:
        start = dna.find(start_codon, start)
        if start == -1:
            break
        indices.append(start)
        start += len(start_codon)
    return indices


def extract_all_proteins(dna_sequence):
    start_codons = find_start_codons(dna_sequence)
    proteins = [dna_to_protein(dna_sequence, codon_idx) for codon_idx in start_codons]
    valid_proteins = {protein for protein in proteins if "-" not in protein}
    return valid_proteins

if __name__ == "__main__":
    [dna_sequence] = load_fasta("ex6.txt")
    reverse_complement_protein = extract_all_proteins(reverse_complement(dna_sequence))
    dna_sequence_protein = extract_all_proteins(dna_sequence)
    all_proteins = reverse_complement_protein | dna_sequence_protein
    print("\n".join(all_proteins))

