"""Neighborhood visualization for SegSure."""

import os
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd


def create_neighborhood_panels(
    high_discord_cells: pd.DataFrame,
    cell_coords: pd.DataFrame,
    coord_cols: Dict[str, str],
    disagreement_metrics: pd.DataFrame,
    neighbors: Dict[str, List[str]],
    output_dir: Optional[str] = None,
    panel_size: float = 100.0,
    logger=None
) -> List[plt.Figure]:
    """Create diagnostic panels for high-discord neighborhoods.
    
    Parameters
    ----------
    high_discord_cells : pd.DataFrame
        Cells in high-discord regions.
    cell_coords : pd.DataFrame
        Cell coordinate data.
    coord_cols : dict
        Coordinate column names.
    disagreement_metrics : pd.DataFrame
        Full disagreement metrics.
    neighbors : dict
        Neighbor mappings.
    output_dir : str, optional
        Directory to save panels.
    panel_size : float, optional
        Size of neighborhood to show. Default is 100.0.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    list
        List of figure objects.
    """
    if logger:
        logger.info(f"Creating {len(high_discord_cells)} neighborhood panels")
    
    figures = []
    
    for idx, (_, row) in enumerate(high_discord_cells.iterrows()):
        cell_id = row["atomx_cell"]
        x = row["x"]
        y = row["y"]
        
        # Get neighborhood
        neighbor_ids = neighbors.get(cell_id, [])
        
        # Create figure
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Plot all cells in region
        x_min, x_max = x - panel_size / 2, x + panel_size / 2
        y_min, y_max = y - panel_size / 2, y + panel_size / 2
        
        region_cells = cell_coords[
            (cell_coords[coord_cols["x"]] >= x_min) &
            (cell_coords[coord_cols["x"]] <= x_max) &
            (cell_coords[coord_cols["y"]] >= y_min) &
            (cell_coords[coord_cols["y"]] <= y_max)
        ]
        
        # Color by discord
        for rid, rrow in region_cells.iterrows():
            metrics = disagreement_metrics[
                disagreement_metrics["atomx_cell"] == rid
            ]
            
            if len(metrics) > 0:
                severity = metrics.iloc[0].get("discord_severity", "unknown")
                color_map = {
                    "none": "green",
                    "low": "yellow",
                    "medium": "orange",
                    "high": "red",
                    "unknown": "gray",
                }
                color = color_map.get(severity, "gray")
            else:
                color = "gray"
            
            # Highlight central cell
            size = 200 if rid == cell_id else 50
            ax.scatter(
                rrow[coord_cols["x"]],
                rrow[coord_cols["y"]],
                c=color,
                s=size,
                alpha=0.7,
                edgecolors="black" if rid == cell_id else None,
            )
        
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")
        ax.set_title(f"Neighborhood of cell {cell_id}")
        ax.set_xlabel(coord_cols["x"])
        ax.set_ylabel(coord_cols["y"])
        
        figures.append(fig)
        
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"neighborhood_{cell_id}.png")
            plt.savefig(output_path, dpi=100, bbox_inches="tight")
            if logger and idx % 10 == 0:
                logger.info(f"Saved neighborhood panel {idx + 1}/{len(high_discord_cells)}")
            plt.close(fig)
    
    if logger:
        logger.info(f"Created {len(figures)} neighborhood panels")
    
    return figures
