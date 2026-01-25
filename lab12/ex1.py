import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_data(filename):
    df = pd.read_csv(filename, sep=';', index_col=0, header=None)
    df.columns = ["V", "B", "P", "pI", "H1", "H2", "SAS", "FA"]
    df = df.replace(",", ".", regex=True)
    df = df.astype(float)
    return df


def calculate_pca(data, k=2):
    X = data.values
    Z = (X - X.mean(axis=0)) / X.std(axis=0)
    cov = np.cov(Z, rowvar=False)
    eig_vals, eig_vecs = np.linalg.eigh(cov)
    explained_variance = eig_vals / eig_vals.sum()

    idx = np.argsort(eig_vals)[::-1]
    eig_vecs = eig_vecs[:, idx][:, :k]
    feature_names = data.columns[idx][:k]

    explained_variance = explained_variance[idx][:k]
    mapped_data = Z.dot(eig_vecs)
    out = pd.DataFrame(mapped_data, index=data.index)
    return feature_names, explained_variance, out


def plot_data(data, explained_variance, feature_names, data_orig):
    if data.shape[1] > 2:
        raise ValueError("Can only plot 2D data")

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.suptitle(f'PCA explained variance ratio (PC1 - {feature_names[0]}, PC2 - {feature_names[1]})')

    groups = {
        "positive": set("RHK"),
        "negative": set("DE"),
        "polar": set("STNQ"),
        "hydrophobic": set("AILMFWYV"),
        "special": set("CUGP"),
    }
    colors = {
        "positive": "red",
        "negative": "blue",
        "polar": "green",
        "hydrophobic": "orange",
        "special": "purple",
    }
    default_color = "gray"

    def label_group(lbl):
        key = str(lbl)
        for gname, letters in groups.items():
            if key in letters:
                return gname
        return None

    group_names = [label_group(lbl) for lbl in data.index]
    color_list = [colors.get(g, default_color) for g in group_names]

    ax.scatter(data[0], data[1], c=color_list, edgecolor="black", s=60, alpha=0.9)

    for i, label in enumerate(data.index):
        pc1 = data_orig.iloc[i][feature_names[0]]
        pc2 = data_orig.iloc[i][feature_names[1]]
        pc_vals = f"({pc1}, {pc2})"
        x, y = data.iloc[i, 0], data.iloc[i, 1]

        ax.annotate(str(label) + pc_vals,
                    xy=(x, y),
                    xytext=(4, 4),
                    textcoords="offset points",
                    ha="left", va="bottom",
                    fontsize=8)


    for gname, col in colors.items():
        ax.scatter([], [], color=col, label=f"{gname} ({''.join(sorted(groups[gname]))})")
    ax.legend(title="Group", loc="best", fontsize=8)

    ax.set_xlabel(f"PC1 ({explained_variance[0] * 100:.2f}%)")
    ax.set_ylabel(f"PC2 ({explained_variance[1] * 100:.2f}%)")

    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()


def main():
    df = load_data("./aaProperties.txt")
    feature_names, explained_variance, transformed_data = calculate_pca(df)
    plot_data(transformed_data, explained_variance, feature_names, df)
    print("""
    WNIOSKI:
    Na wykresie możemy zauważyć, że badane aminokwasy zostały za pomocą
    PCA pogrupowane w zbiory zgodnie ze swoją grupą. Cechami decydującymi o 
    podziale były:
    - SAS (Pole powierzchni dostępnej dla wody)
    - FA (Ułamek pola powierzchni niedostępny dla wody po procesie fałdowania)
    Aminokwasy z ładunkiem elektrycznym mają niskie FA, negatywne mają średnią wartość SAS a pozytywne wysoką,
    widzimy je odseparowane od reszty na wykresie.
    Aminokwasy hydrofobowe wyróżnia wysoka wartość FA, są dobrze odseparowane na wykresie.
    Aminokwasy polarne również są zgrupowane razem, ich wartości FA oraz SAS są mniej więcej na środku przedziału.
    Aminokwasy specjalne nie zostały widocznie odseparowane jako osobna grupa na wykresie.
    Separacja nie jest idealna, gdyż ze wszystkich dostępnych cech wykorzystałem jedynie 2 w celu wizualizacji. 
    """)

if __name__ == '__main__':
    main()

