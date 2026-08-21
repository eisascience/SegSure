"""Neighborhood disagreement metrics for SegSure."""

from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


def identify_neighbors(
    cell_coords: pd.DataFrame,
    coord_cols: Dict[str, str],
    neighbor_distance: float = 50.0,
    logger=None
) -> Dict[str, List[str]]:
    """Identify neighboring cells within distance threshold.
    
    Parameters
    ----------
    cell_coords : pd.DataFrame
        Cell coordinate data.
    coord_cols : dict
        Coordinate column names (x, y, z).
    neighbor_distance : float, optional
        Maximum distance for neighborhood. Default is 50.0.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    dict
        Dictionary mapping cell IDs to list of neighbor cell IDs.
    """
    if logger:
        logger.info(f"Identifying neighbors (distance <= {neighbor_distance})")
    
    coords = cell_coords[[coord_cols["x"], coord_cols["y"]]].values
    distances = cdist(coords, coords, metric="euclidean")
    
    neighbors = {}
    for i, cell_id in enumerate(cell_coords.index):
        neighbor_indices = np.where(
            (distances[i] > 0) & (distances[i] <= neighbor_distance)
        )[0]
        neighbors[cell_id] = [cell_coords.index[j] for j in neighbor_indices]
    
    if logger:
        avg_neighbors = np.mean([len(v) for v in neighbors.values()])
        logger.info(f"Average neighbors per cell: {avg_neighbors:.1f}")
    
    return neighbors


def compute_neighborhood_discord(
    assignment_comparisons: pd.DataFrame,
    neighbors_atomx: Dict[str, List[str]],
    neighbors_proseg: Dict[str, List[str]],
    cell_matches: pd.DataFrame,
    logger=None
) -> pd.DataFrame:
    """Compute disagreement density in neighborhoods.
    
    Parameters
    ----------
    assignment_comparisons : pd.DataFrame
        Per-cell discord metrics.
    neighbors_atomx : dict
        AtoMx neighbor mappings.
    neighbors_proseg : dict
        Proseg neighbor mappings.
    cell_matches : pd.DataFrame
        Cell matching results.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Neighborhood disagreement metrics.
    """
    if logger:
        logger.info("Computing neighborhood discord")
    
    # Create mapping from cell ID to discord
    discord_map = dict(
        zip(assignment_comparisons["atomx_cell"], 
            assignment_comparisons["discord_severity"])
    )
    
    results = []
    
    for _, match_row in cell_matches.iterrows():
        atomx_id = match_row["atomx_cell"]
        proseg_id = match_row["proseg_cell"]
        
        # Get neighbors
        atomx_neighbors = neighbors_atomx.get(atomx_id, [])
        proseg_neighbors = neighbors_proseg.get(proseg_id, [])
        
        # Count neighbor discord
        high_discord_neighbors = sum(
            1 for n in atomx_neighbors 
            if discord_map.get(n) == "high"
        )
        
        results.append({
            "atomx_cell": atomx_id,
            "proseg_cell": proseg_id,
            "n_neighbors_atomx": len(atomx_neighbors),
            "n_neighbors_proseg": len(proseg_neighbors),
            "high_discord_neighbors": high_discord_neighbors,
            "neighbor_discord_density": (
                high_discord_neighbors / max(len(atomx_neighbors), 1)
            ),
        })
    
    return pd.DataFrame(results) if results else pd.DataFrame()


def identify_high_discord_regions(
    neighborhood_discord: pd.DataFrame,
    cell_coords: pd.DataFrame,
    coord_cols: Dict[str, str],
    threshold: float = 0.5,
    logger=None
) -> pd.DataFrame:
    """Identify regions with high neighborhood discord.
    
    Parameters
    ----------
    neighborhood_discord : pd.DataFrame
        Neighborhood discord metrics.
    cell_coords : pd.DataFrame
        Cell coordinate data.
    coord_cols : dict
        Coordinate column names.
    threshold : float, optional
        Discord density threshold. Default is 0.5.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        High-discord cells with location and neighborhood info.
    """
    if logger:
        logger.info(f"Identifying high-discord regions (threshold={threshold})")
    
    high_discord = neighborhood_discord[
        neighborhood_discord["neighbor_discord_density"] >= threshold
    ].copy()
    
    # Add coordinates
    for idx in high_discord.index:
        atomx_id = high_discord.loc[idx, "atomx_cell"]
        if atomx_id in cell_coords.index:
            coord_row = cell_coords.loc[atomx_id]
            high_discord.loc[idx, "x"] = coord_row[coord_cols["x"]]
            high_discord.loc[idx, "y"] = coord_row[coord_cols["y"]]
    
    if logger:
        logger.info(f"Found {len(high_discord)} cells in high-discord regions")
    
    return high_discord
