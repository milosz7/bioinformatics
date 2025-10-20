def load_input(filename):
    with open(filename, "r") as f:
        data = f.read().strip()
    return data


def transcribe_dna_to_rna(dna):
    return dna.replace("T", "U")


if __name__ == "__main__":
    dna = load_input("ex3.txt")
    print(transcribe_dna_to_rna(dna))
