import argparse

from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt


def read_clustal(path):
    alignment = AlignIO.read(path, "clustal")
    return alignment


def compute_distance(alignment):
    distance_constructor = DistanceCalculator("blosum62")
    return distance_constructor.get_distance(alignment)


def build_tree(distances):
    tree = DistanceTreeConstructor()
    return tree.upgma(distances)


def save_tree(tree, figsize=(20, 10), path="tree.png"):
    fig = plt.figure(figsize=figsize)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, branch_labels=lambda c: round(c.branch_length, 3), axes=axes, do_show=False)
    print(f"Saved tree to {path}")
    fig.savefig(path)


def calculate_distance_from(tree, species):
    terminals = tree.get_terminals()
    computed_distances = []
    for terminal in terminals:
        distance = tree.distance(terminal, species)
        computed_distances.append(((terminal, species), distance))
    computed_distances = sorted(computed_distances, key=lambda x: x[1])
    return computed_distances


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="path CLUSTAL to alignment")
    args = parser.parse_args()

    alignment = read_clustal(args.path)
    distances = compute_distance(alignment)
    tree = build_tree(distances)
    Phylo.draw_ascii(tree)
    save_tree(tree)
    distances_from_homo = calculate_distance_from(tree, "Homo-sapiens")
    for (other, homo), distance in distances_from_homo:
        print(f"{homo} - {other} | Distance: {distance:.4f}")

    print("""
    Wnioski:
    Na podstawie drzewa filogenetycznego opartego o gen hbb (hemoglobin subunit beta) możemy stwierdzić, że:
    - Największe podobieństwo (najmniejszy dystans) spośród badanych par człowiek - zwierze występuje człowiekiem
      a orangutanem / gibonem (taki sam dystans)
    - Najmniejsze podobieństwo występuje pomiędzy człowiekiem a szynszylą małą
    
    Większe podobieństwo tego białka u małp w porównaniu z gryzoniami jest wynikiem jakiego się spodziewałem.
    
    Dlaczego użyliśmy sekwencji aminokwasów a nie nukleotydów? 
    Kod genetyczny jest zdegenerowany - moglibyśmy otrzymać zaburzone wyniki pomimo kodowania 
    identycznych aminokwasów w sekwencji. W przypadku odległych od siebie organizmów takich mutacji mogło zajść
    wiele pogarszając jakość porównania.
    """)
