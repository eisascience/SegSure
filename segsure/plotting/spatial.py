"""Spatial visualization for SegSure."""

import os
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_uncertainty_heatmap(
    disagreement_metrics: pd.DataFrame,
    cell_coords: pd.DataFrame,
    coord_cols: Dict[str, str],
    color_metric: str = "discord_severity",
    output_path: Optional[str] = None,
    figsize: tuple = (12, 10),
    logger=None
) -> Optional[plt.Figure]:
    """Plot spatial uncertainty heatmap.
    
    Parameters
    ----------
    disagreement_metrics : pd.DataFrame
        Disagreement metrics with cell IDs.
    cell_coords : pd.DataFrame
        Cell coordinate data.
    coord_cols : dict
        Coordinate column names.
    color_metric : str, optional
        Metric to use for coloring. Default is discord_severity.
    output_path : str, optional
        Path to save figure.
    figsize : tuple, optional
        Figure size.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    matplotlib.Figure or None
        Figure object if output_path not provided.
    """
    if logger:
        logger.info(f"Creating uncertainty heatmap (colored by {color_metric})")
    
    # Merge coordinates with metrics
    plot_data = disagreement_metrics.copy()
    
    # Add coordinates
    for idx in plot_data.index:
        atomx_id = plot_data.loc[idx, "atomx_cell"]
        if atomx_id in cell_coords.index:
            coord_row = cell_coords.loc[atomx_id]
            plot_data.loc[idx, "x"] = coord_row[coord_cols["x"]]
            plot_data.loc[idx, "y"] = coord_row[coord_cols["y"]]
    
    plot_data = plot_data.dropna(subset=["x", "y"])
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Determine colormap
    if color_metric == "discord_severity":
        # Categorical
        severity_order = ["none", "low", "medium", "high"]
        colors = ["green", "yellow", "orange", "red"]
        
        # Plot by category
        for severity, color in zip(severity_order, colors):
            mask = plot_data[color_metric] == severity
            if mask.any():
                ax.scatter(
                    plot_data.loc[mask, "x"],
                    plot_data.loc[mask, "y"],
                    c=color,
                    label=severity,
                    alpha=0.6,
                    s=50,
                )
    else:
        # Continuous
        scatter = ax.scatter(
            plot_data["x"],
            plot_data["y"],
            c=plot_data[color_metric],
            cmap="RdYlGn_r",
            s=50,
            alpha=0.6,
        )
        plt.colorbar(scatter, ax=ax, label=color_metric)
    
    ax.set_xlabel(coord_cols["x"])
    ax.set_ylabel(coord_cols["y"])
    ax.set_title(f"Segmentation Uncertainty Map ({color_metric})")
    ax.legend()
    ax.set_aspect("equal")
    
    plt.tight_layout()
    
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        if logger:
            logger.info(f"Saved uncertainty heatmap to {output_path}")
        plt.close()
        return None
    
    return fig


def plot_discord_distribution(
    disagreement_metrics: pd.DataFrame,
    output_path: Optional[str] = None,
    figsize: tuple = (10, 6),
    logger=None
) -> Optional[plt.Figure]:
    """Plot distribution of discord metrics.
    
    Parameters
    ----------
    disagreement_metrics : pd.DataFrame
        Disagreement metrics.
    output_path : str, optional
        Path to save figure.
    figsize : tuple, optional
        Figure size.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    matplotlib.Figure or None
        Figure object if output_path not provided.
    """
    if logger:
        logger.info("Creating discord distribution plot")
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    
    # Jaccard index
    if "jaccard_index" in disagreement_metrics.columns:
        axes[0, 0].hist(
            disagreement_metrics["jaccard_index"].dropna(), 
            bins=30, 
            color="blue", 
            alpha=0.7
        )
        axes[0, 0].set_xlabel("Jaccard Index")
        axes[0, 0].set_ylabel("Count")
        axes[0, 0].set_title("Transcript Overlap (Jaccard Index)")
    
    # Discord severity
    if "discord_severity" in disagreement_metrics.columns:
        severity_counts = disagreement_metrics["discord_severity"].value_counts()
        axes[0, 1].bar(severity_counts.index, severity_counts.values, color="orange", alpha=0.7)
        axes[0, 1].set_xlabel("Discord Severity")
        axes[0, 1].set_ylabel("Count")
        axes[0, 1].set_title("Discord Severity Distribution")
    
    # Centroid distance
    if "centroid_distance" in disagreement_metrics.columns:
        axes[1, 0].hist(
            disagreement_metrics["centroid_distance"].dropna(), 
            bins=30, 
            color="red", 
            alpha=0.7
        )
        axes[1, 0].set_xlabel("Centroid Distance")
        axes[1, 0].set_ylabel("Count")
        axes[1, 0].set_title("Centroid Distance Distribution")
    
    # Neighbor discord density
    if "neighbor_discord_density" in disagreement_metrics.columns:
        axes[1, 1].hist(
            disagreement_metrics["neighbor_discord_density"].dropna(), 
            bins=30, 
            color="purple", 
            alpha=0.7
        )
        axes[1, 1].set_xlabel("Neighbor Discord Density")
        axes[1, 1].set_ylabel("Count")
        axes[1, 1].set_title("Neighborhood Discord Density")
    
    plt.tight_layout()
    
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        if logger:
            logger.info(f"Saved discord distribution plot to {output_path}")
        plt.close()
        return None
    
    return fig
