### Import necessary libraries
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

from pointpats.distance_statistics import g_test, f_test, l_test, k_test

sc.logging.print_header()

### Centralized path manager
def get_paths(dataset: str):
    """
    Return dataset-specific paths.

    Args:
        dataset (str): One of ["SN", "D", "AIH", "comparison"]
    """
    root = "/path/to/your/data/"
    return {
        "counts": root + "counts/",
        "meta": root + "meta/",
        "outdir": root + f"{dataset}_output/"
    }

### Define the functions needed to run the analysis ###

def canonical_pair(row):
    types = sorted([row['cell_type_1'], row['cell_type_2']])
    return f"{types[0]} ~ {types[1]}"

def split_by_fov(adata, fov_col='fov'):
    return [adata[adata.obs[fov_col] == fov].copy() for fov in adata.obs[fov_col].unique()]

def run_neighbourhood_enrichment_per_fov(adata, sample_name, radius=50, cluster_key="cell_type", fov_col="fov"):
    results = []
    fovs = adata.obs[fov_col].unique()

    for fov in fovs:
        adata_fov = adata[adata.obs[fov_col] == fov].copy()
        adata_fov.obs[cluster_key] = adata_fov.obs["final_cell_types"].astype("category")

        sq.gr.spatial_neighbors(adata_fov, radius=radius, coord_type="generic")
        sq.gr.nhood_enrichment(adata_fov, cluster_key=cluster_key)

        z = adata_fov.uns[f"{cluster_key}_nhood_enrichment"]['zscore']
        labels = adata_fov.obs[cluster_key].cat.categories.tolist()
        z_df = pd.DataFrame(z, index=labels, columns=labels)

        z_melt = z_df.reset_index().melt(id_vars='index')
        z_melt.columns = ['cell_type_1', 'cell_type_2', 'zscore']
        z_melt['fov'] = fov
        z_melt['sample'] = sample_name

        results.append(z_melt)

    return pd.concat(results, ignore_index=True)

def check_problematic_fovs_from_df(zscore_df, adata, fov_col="fov", cluster_key="cell_type", min_cells=50):
    problematic_fovs = []
    fov_counts = adata.obs[fov_col].value_counts()
    grouped = zscore_df.groupby('fov')

    for fov, group in grouped:
        issue = None
        n_cells = fov_counts.get(fov, 0)

        if n_cells < min_cells:
            issue = "too few cells"
            bad_pairs = []
        elif group['zscore'].isna().any():
            issue = "NaNs in z-score matrix"
            bad_pairs = group[group['zscore'].isna()][['cell_type_1', 'cell_type_2']].drop_duplicates().values.tolist()
        else:
            continue

        problematic_fovs.append({
            "FOV": fov,
            "n_cells": n_cells,
            "issue": issue,
            "affected_pairs": bad_pairs
        })

    return pd.DataFrame(problematic_fovs)

def neighbourhood_enrichment(adata, sample_name, radius=50, cluster_key="cell_type", save_plots=True, outdir=None):
    adata.obs[cluster_key] = adata.obs["final_cell_types"].astype("category")
    sq.gr.spatial_neighbors(adata, radius=radius, coord_type="generic")
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)

    if save_plots and outdir is not None:
        sq.pl.nhood_enrichment(
            adata,
            cluster_key=cluster_key,
            figsize=(5, 5),
            title=f"Neighborhood enrichment - {sample_name} (radius {radius})",
            save=f"{outdir}{sample_name}_neighborhood_enrichment_radius{radius}.png"
        )

    z = adata.uns[f"{cluster_key}_nhood_enrichment"]['zscore']
    labels = adata.obs[cluster_key].cat.categories.tolist()
    return pd.DataFrame(z, index=labels, columns=labels)

def average_neighbourhood_plot(zscore_dfs, output_path_png, output_path_csv, 
                               title="Average Neighborhood Enrichment", cmap="bwr"):
    z_avg = sum(zscore_dfs) / len(zscore_dfs)
    z_avg.to_csv(output_path_csv)

    plt.figure(figsize=(15, 15))
    ax = sns.heatmap(
        z_avg, cmap=cmap, center=0, annot=True, fmt=".2f", square=True,
        cbar_kws={"label": "Z-score"}
    )
    ax.set_title(title, fontsize=20, weight="bold")
    ax.set_xlabel("Neighbor cell type", fontsize=16, weight="bold")
    ax.set_ylabel("Index cell type", fontsize=16, weight="bold")
    ax.tick_params(axis="x", labelrotation=90, labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    plt.tight_layout()
    plt.savefig(output_path_png, dpi=300)
    plt.close()
    return z_avg

def plot_neighbourhood_difference(diff_df, title, output_path, cmap="bwr"):
    plt.figure(figsize=(10, 10))
    sns.heatmap(
        diff_df, cmap='PRGn', center=0, annot=True, fmt=".2f", square=True,
        cbar_kws={"label": "Δ Z-score"}
    )
    plt.title(title)
    plt.xlabel("Neighbor cell type")
    plt.ylabel("Index cell type")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)

def plot_centrality_scores(
    adata, sample_name, cluster_key="cell_type", celltype_palette=None,
    save_path=None, radius=50, figsize=(25, 8)
):
    sq.gr.centrality_scores(adata, cluster_key=cluster_key)
    sq.pl.centrality_scores(adata, cluster_key=cluster_key, figsize=(25, 6))
    plt.close()

    centrality_key = f"{cluster_key}_centrality_scores"
    centrality_df = pd.DataFrame(adata.uns[centrality_key])

    cols = centrality_df.columns.tolist()
    if "average_clustering" in cols:
        cols.insert(2, cols.pop(cols.index("average_clustering")))
    centrality_df = centrality_df[cols]

    cell_types = centrality_df.index.tolist()
    if celltype_palette is None:
        palette = plt.get_cmap("tab20").colors
        celltype_palette = dict(zip(cell_types, palette[:len(cell_types)]))

    fig, axes = plt.subplots(1, len(centrality_df.columns), figsize=figsize, sharey=True)
    for i, score in enumerate(centrality_df.columns):
        ax = axes[i]
        sorted_df = centrality_df.sort_values(by=score, ascending=True)
        bar_colors_sorted = [celltype_palette[ct] for ct in sorted_df.index]
        ax.barh(sorted_df.index, sorted_df[score], color=bar_colors_sorted)
        ax.set_title(score.replace('_', ' ').title())
        ax.set_xlabel("Score")
        ax.set_ylabel("Cell Type" if i == 0 else "")

    plt.tight_layout()
    if save_path is None:
        save_path = f"{sample_name}_centrality_scores_radius{radius}.png"
    plt.savefig(save_path, dpi=300)
    return centrality_df

def plot_average_centrality(
    centrality_dfs: list, sample_name: str, celltype_palette: dict,
    save_path: str = None, radius=50, figsize=(25, 8)
):
    centrality_dfs = [df.sort_index() for df in centrality_dfs]
    centrality_avg_df = sum(centrality_dfs) / len(centrality_dfs)

    cols = centrality_avg_df.columns.tolist()
    if "average_clustering" in cols:
        cols.insert(2, cols.pop(cols.index("average_clustering")))
    centrality_avg_df = centrality_avg_df[cols]

    fig, axes = plt.subplots(1, len(centrality_avg_df.columns), figsize=figsize, sharey=True)
    for i, score in enumerate(centrality_avg_df.columns):
        ax = axes[i]
        sorted_df = centrality_avg_df.sort_values(by=score, ascending=True)
        bar_colors_sorted = [celltype_palette[ct] for ct in sorted_df.index]
        ax.barh(sorted_df.index, sorted_df[score], color=bar_colors_sorted)
        ax.set_title(f"{sample_name} Avg {score.replace('_', ' ').title()}")
        ax.set_xlabel("Score")
        ax.set_ylabel("Cell Type" if i == 0 else "")

    plt.tight_layout()
    if save_path is None:
        save_path = f"{sample_name}_average_centrality_scores_radius{radius}.png"
    plt.savefig(save_path, dpi=300)
    return centrality_avg_df

def plot_ripley_all_clusters(coords, clusters, sample_name, celltype_color_map,
                             mode='G', interval=np.linspace(0, 80, 50), n_simulations=99,
                             save_plots=True, outdir=None):
    unique_clusters = np.unique(clusters)
    plt.figure(figsize=(9, 6))

    for clust in unique_clusters:
        mask = clusters == clust
        cluster_coords = coords[mask]

        if mode == 'G':
            res = g_test(cluster_coords, support=interval, n_simulations=n_simulations, keep_simulations=True)
        elif mode == 'F':
            res = f_test(cluster_coords, support=interval, n_simulations=n_simulations, keep_simulations=True)
        elif mode == 'K':
            res = k_test(cluster_coords, support=interval, n_simulations=n_simulations, keep_simulations=True)
        elif mode == 'L':
            res = l_test(cluster_coords, support=interval, n_simulations=n_simulations, keep_simulations=True)
        else:
            raise ValueError("mode must be one of ['G', 'F', 'K', 'L']")

        color = celltype_color_map.get(clust, None)
        obs_line, = plt.plot(res.support, res.statistic, label=f'{clust} observed', color=color)

        if hasattr(res, 'simulations') and res.simulations is not None:
            lower = np.percentile(res.simulations, 2.5, axis=0)
            upper = np.percentile(res.simulations, 97.5, axis=0)
            median = np.median(res.simulations, axis=0)
            plt.fill_between(res.support, lower, upper, alpha=0.1, color=obs_line.get_color())
            plt.plot(res.support, median, linestyle='--', alpha=0.7, color=obs_line.get_color(), label=f'{clust} median sim')

    plt.xlabel('Distance')
    ylabels = {
        'G': 'G function\n(CDF of NN distance)',
        'F': 'F function\n(CDF of empty space distance)',
        'K': 'K function',
        'L': 'L function\n(square-root transformed K)'
    }
    plt.ylabel(ylabels.get(mode, 'Value'))
    plt.title(f"Ripley's {mode} function for all cell types - {sample_name}")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    if save_plots and outdir is not None:
        plt.savefig(f"{outdir}{sample_name}_Ripley_{mode}_all_clusters.png", dpi=300)

def process_sims_stat(ripley_dict, source_label):
    df = ripley_dict['sims_stat'].copy()
    df['source'] = source_label
    return df

