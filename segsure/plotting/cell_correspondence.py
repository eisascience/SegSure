"""Cell correspondence visualization for SegSure."""

import os
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from shapely.geometry import Polygon


# Color mapping for relationship types
RELATIONSHIP_COLORS = {
    "one_to_one": "#2ecc71",     # Green
    "split": "#f39c12",          # Orange
    "merge": "#e74c3c",          # Red
    "lost_atomx": "#95a5a6",     # Gray
    "new_proseg": "#3498db",     # Blue
    "complex": "#9b59b6",        # Purple
}


def plot_polygon(ax, polygon: Polygon, facecolor='lightblue', edgecolor='blue', alpha=0.3, linewidth=1):
    """Plot a Polygon on matplotlib axes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    polygon : Polygon
        Shapely polygon to plot.
    facecolor : str
        Face color.
    edgecolor : str
        Edge color.
    alpha : float
        Transparency.
    linewidth : float
        Line width.
    """
    if polygon.geom_type == "Polygon":
        x, y = polygon.exterior.xy
        ax.fill(x, y, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, linewidth=linewidth)


def create_correspondence_figure(
    atomx_polygons: Dict,
    proseg_polygons: Dict,
    matches: pd.DataFrame,
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (16, 12),
    logger=None
) -> plt.Figure:
    """Create diagnostic figure for cell correspondence.
    
    Shows AtoMx boundaries, Proseg boundaries, overlays, and relationships colored by type.
    
    Parameters
    ----------
    atomx_polygons : dict
        Mapping of AtoMx cell IDs to Polygon objects.
    proseg_polygons : dict
        Mapping of Proseg cell IDs to Polygon objects.
    matches : pd.DataFrame
        Correspondence results with relationship_type column.
    output_path : str, optional
        Path to save figure.
    figsize : tuple
        Figure size.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure.
    """
    fig = plt.figure(figsize=figsize)
    
    # Get relationship type mapping
    relationship_type_col = "relationship_type"
    if relationship_type_col not in matches.columns:
        if logger:
            logger.warning("No relationship_type column found, using raw matches")
        relationships = {}
    else:
        # Create mapping of (atomx_cell, proseg_cell) -> relationship_type
        relationships = {}
        for _, row in matches.iterrows():
            atom_id = row.get("atomx_cell")
            proseg_id = row.get("proseg_cell")
            if atom_id and proseg_id:
                relationships[(atom_id, proseg_id)] = row.get(relationship_type_col, "unknown")
    
    # Create subplots
    # 1. AtoMx boundaries
    ax1 = plt.subplot(2, 3, 1)
    ax1.set_title("AtoMx Segmentation", fontsize=12, fontweight='bold')
    ax1.set_aspect('equal')
    
    for cell_id, polygon in atomx_polygons.items():
        plot_polygon(ax1, polygon, facecolor='lightblue', edgecolor='blue', alpha=0.5, linewidth=1)
    
    ax1.set_xlabel("X coordinate")
    ax1.set_ylabel("Y coordinate")
    
    # 2. Proseg boundaries
    ax2 = plt.subplot(2, 3, 2)
    ax2.set_title("Proseg Segmentation", fontsize=12, fontweight='bold')
    ax2.set_aspect('equal')
    
    for cell_id, polygon in proseg_polygons.items():
        plot_polygon(ax2, polygon, facecolor='lightcoral', edgecolor='red', alpha=0.5, linewidth=1)
    
    ax2.set_xlabel("X coordinate")
    ax2.set_ylabel("Y coordinate")
    
    # 3. Overlay
    ax3 = plt.subplot(2, 3, 3)
    ax3.set_title("Overlay (AtoMx blue, Proseg red)", fontsize=12, fontweight='bold')
    ax3.set_aspect('equal')
    
    for cell_id, polygon in atomx_polygons.items():
        plot_polygon(ax3, polygon, facecolor='lightblue', edgecolor='blue', alpha=0.3, linewidth=0.5)
    
    for cell_id, polygon in proseg_polygons.items():
        plot_polygon(ax3, polygon, facecolor='lightcoral', edgecolor='red', alpha=0.3, linewidth=0.5)
    
    ax3.set_xlabel("X coordinate")
    ax3.set_ylabel("Y coordinate")
    
    # 4-6. Relationship-specific overlays
    relationship_types = set(RELATIONSHIP_COLORS.keys())
    rel_idx = 0
    
    for rel_type in ["one_to_one", "split", "merge"]:
        ax = plt.subplot(2, 3, 4 + rel_idx)
        ax.set_title(f"{rel_type.replace('_', ' ').title()}", fontsize=12, fontweight='bold')
        ax.set_aspect('equal')
        
        # Plot all boundaries with low opacity
        for cell_id, polygon in atomx_polygons.items():
            plot_polygon(ax, polygon, facecolor='lightblue', edgecolor='blue', alpha=0.1, linewidth=0.5)
        
        for cell_id, polygon in proseg_polygons.items():
            plot_polygon(ax, polygon, facecolor='lightcoral', edgecolor='red', alpha=0.1, linewidth=0.5)
        
        # Highlight matching pairs for this relationship type
        for (atom_id, proseg_id), rel in relationships.items():
            if rel == rel_type:
                if atom_id in atomx_polygons:
                    plot_polygon(ax, atomx_polygons[atom_id], 
                               facecolor=RELATIONSHIP_COLORS[rel_type], 
                               edgecolor='darkgreen' if rel_type == 'one_to_one' else 'black',
                               alpha=0.6, linewidth=1)
                if proseg_id in proseg_polygons:
                    plot_polygon(ax, proseg_polygons[proseg_id], 
                               facecolor=RELATIONSHIP_COLORS[rel_type], 
                               edgecolor='black', alpha=0.3, linewidth=1)
        
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
        rel_idx += 1
    
    # Add lost and new cell overlay
    ax_lost = plt.subplot(2, 3, 6)
    ax_lost.set_title("Lost (gray) / New (blue)", fontsize=12, fontweight='bold')
    ax_lost.set_aspect('equal')
    
    # Plot all boundaries
    for cell_id, polygon in atomx_polygons.items():
        plot_polygon(ax_lost, polygon, facecolor='lightblue', edgecolor='blue', alpha=0.1, linewidth=0.5)
    
    for cell_id, polygon in proseg_polygons.items():
        plot_polygon(ax_lost, polygon, facecolor='lightcoral', edgecolor='red', alpha=0.1, linewidth=0.5)
    
    # Highlight lost and new cells
    for _, row in matches.iterrows():
        rel_type = row.get(relationship_type_col, "unknown")
        if rel_type == "lost_atomx":
            atom_id = row.get("atomx_cell")
            if atom_id in atomx_polygons:
                plot_polygon(ax_lost, atomx_polygons[atom_id], 
                           facecolor=RELATIONSHIP_COLORS["lost_atomx"], 
                           edgecolor='black', alpha=0.6, linewidth=1)
        elif rel_type == "new_proseg":
            proseg_id = row.get("proseg_cell")
            if proseg_id in proseg_polygons:
                plot_polygon(ax_lost, proseg_polygons[proseg_id], 
                           facecolor=RELATIONSHIP_COLORS["new_proseg"], 
                           edgecolor='black', alpha=0.6, linewidth=1)
    
    ax_lost.set_xlabel("X coordinate")
    ax_lost.set_ylabel("Y coordinate")
    
    plt.tight_layout()
    
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        if logger:
            logger.info(f"Saved correspondence figure to: {output_path}")
    
    return fig


def create_relationship_summary_figure(
    matches: pd.DataFrame,
    output_path: Optional[str] = None,
    logger=None
) -> plt.Figure:
    """Create summary figure for relationship statistics.
    
    Parameters
    ----------
    matches : pd.DataFrame
        Correspondence results with relationship_type column.
    output_path : str, optional
        Path to save figure.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    relationship_type_col = "relationship_type"
    
    # 1. Relationship type distribution
    ax = axes[0, 0]
    if relationship_type_col in matches.columns:
        counts = matches[relationship_type_col].value_counts()
        colors = [RELATIONSHIP_COLORS.get(rel, '#cccccc') for rel in counts.index]
        counts.plot(kind='bar', ax=ax, color=colors)
        ax.set_title("Relationship Type Distribution")
        ax.set_ylabel("Count")
        ax.set_xlabel("Relationship Type")
        ax.tick_params(axis='x', rotation=45)
    
    # 2. IoU distribution for one-to-one
    ax = axes[0, 1]
    if "iou" in matches.columns and relationship_type_col in matches.columns:
        one_to_one_iou = matches[matches[relationship_type_col] == "one_to_one"]["iou"]
        if len(one_to_one_iou) > 0:
            ax.hist(one_to_one_iou, bins=20, color=RELATIONSHIP_COLORS["one_to_one"], alpha=0.7, edgecolor='black')
            ax.set_title(f"IoU Distribution (One-to-One, n={len(one_to_one_iou)})")
            ax.set_xlabel("Intersection over Union (IoU)")
            ax.set_ylabel("Count")
            ax.axvline(one_to_one_iou.median(), color='red', linestyle='--', label=f'Median: {one_to_one_iou.median():.3f}')
            ax.legend()
    
    # 3. Intersection area vs union area
    ax = axes[1, 0]
    if "intersection_area" in matches.columns and "union_area" in matches.columns:
        ax.scatter(matches["union_area"], matches["intersection_area"], alpha=0.5)
        ax.set_xlabel("Union Area")
        ax.set_ylabel("Intersection Area")
        ax.set_title("Intersection vs Union Area")
        ax.set_yscale('log')
        ax.set_xscale('log')
    
    # 4. Coverage statistics
    ax = axes[1, 1]
    if "coverage_a" in matches.columns and "coverage_b" in matches.columns and relationship_type_col in matches.columns:
        one_to_one = matches[matches[relationship_type_col] == "one_to_one"]
        if len(one_to_one) > 0:
            coverage_data = [
                one_to_one["coverage_a"].values,
                one_to_one["coverage_b"].values,
            ]
            bp = ax.boxplot(coverage_data)
            ax.set_xticklabels(["AtoMx Coverage", "Proseg Coverage"])
            ax.set_title("Coverage Distribution (One-to-One)")
            ax.set_ylabel("Coverage Fraction")
            ax.set_ylim(0, 1)
    
    plt.tight_layout()
    
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        if logger:
            logger.info(f"Saved summary figure to: {output_path}")
    
    return fig
