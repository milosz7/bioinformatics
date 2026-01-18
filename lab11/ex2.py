import argparse
from Bio.Cluster import kcluster
import pandas as pd


def load_and_parse_data(path, sep='\t'):
    df = pd.read_csv(path, sep=sep)
    new_columns = df.pop("GENE")
    df = df.T
    df.replace(to_replace='x', value=0.0, inplace=True)
    df.columns = new_columns
    return df


def cluster_timestamps(df, k=2):
    cluster_ids = kcluster(df.values, k)
    return cluster_ids


def print_results(cluster_ids, df, k=2):
    df['cluster_id'] = cluster_ids
    for n in range(k):
        print(f"Cluster {n}:")
        indices = df['cluster_id'] == n
        print(", ".join(df.loc[indices].index))


def main():
    parser = argparse.ArgumentParser(description='Cluster timestamps')
    parser.add_argument("--path", help="Path to gene expression file", default="./yeast_expression.txt")
    parser.add_argument("--sep", help="Separator between columns", default="\t")
    parser.add_argument("--k", help="Number of clusters", default=2)
    args = parser.parse_args()
    df = load_and_parse_data(args.path)
    cluster_ids, _error, _n_found = cluster_timestamps(df, args.k)
    print_results(cluster_ids, df, args.k)


if __name__ == "__main__":
    main()
