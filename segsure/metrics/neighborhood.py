"""Neighborhood disagreement metrics for SegSure."""

from typing import Dict, List, Optional, Set

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


def compute_mapped_neighbor_jaccard(
    atomx_neighbors: Dict[str, List[str]],
    proseg_neighbors: Dict[str, List[str]],
    cell_matches: pd.DataFrame,
    logger=None
) -> float:
    """Compute Jaccard similarity on mapped neighbor sets.
    
    Parameters
    ----------
    atomx_neighbors : dict
        AtoMx neighbor mappings (cell_id -> list of neighbor_ids).
    proseg_neighbors : dict
        Proseg neighbor mappings.
    cell_matches : pd.DataFrame
        Cell matching results with atomx_cell and proseg_cell columns.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    float
        Jaccard similarity (0-1) for all matched neighborhoods.
    """
    if logger:
        logger.info("Computing mapped-neighbor Jaccard similarity")
    
    # Build correspondence map
    atomx_to_proseg = {}
    for _, row in cell_matches.iterrows():
        atomx_id = row["atomx_cell"]
        proseg_id = row["proseg_cell"]
        atomx_to_proseg[atomx_id] = proseg_id
    
    proseg_to_atomx = {v: k for k, v in atomx_to_proseg.items()}
    
    jaccard_scores = []
    
    for atomx_cell, atomx_neighbor_ids in atomx_neighbors.items():
        proseg_cell = atomx_to_proseg.get(atomx_cell)
        if proseg_cell is None:
            continue
        
        proseg_neighbor_ids = proseg_neighbors.get(proseg_cell, [])
        
        # Map AtoMx neighbors to Proseg via correspondence
        mapped_atomx_neighbors = set()
        for nb_id in atomx_neighbor_ids:
            mapped_id = atomx_to_proseg.get(nb_id)
            if mapped_id is not None:
                mapped_atomx_neighbors.add(mapped_id)
        
        proseg_neighbors_set = set(proseg_neighbor_ids)
        
        # Jaccard
        if len(mapped_atomx_neighbors | proseg_neighbors_set) > 0:
            jaccard = len(mapped_atomx_neighbors & proseg_neighbors_set) / len(
                mapped_atomx_neighbors | proseg_neighbors_set
            )
            jaccard_scores.append(jaccard)
    
    if jaccard_scores:
        result = np.mean(jaccard_scores)
        if logger:
            logger.info(
                f"Mapped-neighbor Jaccard: {result:.3f} "
                f"(across {len(jaccard_scores)} matched cells)"
            )
        return result
    
    return np.nan


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
        Neighborhood disagreement metrics with columns:
        - atomx_cell, proseg_cell
        - n_neighbors_atomx, n_neighbors_proseg
        - shared_neighbors, gained_neighbors, lost_neighbors
        - neighbor_jaccard
        - high_discord_neighbors
        - neighbor_discord_density
    """
    if logger:
        logger.info("Computing neighborhood discord metrics")
    
    # Build cell correspondence map
    atomx_to_proseg = {}
    proseg_to_atomx = {}
    for _, row in cell_matches.iterrows():
        atomx_id = row["atomx_cell"]
        proseg_id = row["proseg_cell"]
        atomx_to_proseg[atomx_id] = proseg_id
        proseg_to_atomx[proseg_id] = atomx_id
    
    # Create discord map
    discord_map = {}
    if "discord_severity" in assignment_comparisons.columns:
        for _, row in assignment_comparisons.iterrows():
            atomx_id = row.get("atomx_cell")
            severity = row.get("discord_severity", "unknown")
            discord_map[atomx_id] = severity
    
    results = []
    
    for _, match_row in cell_matches.iterrows():
        atomx_id = match_row["atomx_cell"]
        proseg_id = match_row["proseg_cell"]
        
        # Get neighbors
        atomx_neighbors = neighbors_atomx.get(atomx_id, [])
        proseg_neighbors = neighbors_proseg.get(proseg_id, [])
        
        # Map neighbors
        mapped_atomx_neighbors = set()
        for nb_id in atomx_neighbors:
            mapped_id = atomx_to_proseg.get(nb_id)
            if mapped_id is not None:
                mapped_atomx_neighbors.add(mapped_id)
        
        proseg_neighbors_set = set(proseg_neighbors)
        
        # Shared, gained, lost
        shared = mapped_atomx_neighbors & proseg_neighbors_set
        gained = proseg_neighbors_set - mapped_atomx_neighbors  # Proseg only
        lost = mapped_atomx_neighbors - proseg_neighbors_set    # AtoMx only
        
        # Neighbor Jaccard
        neighbor_union = mapped_atomx_neighbors | proseg_neighbors_set
        neighbor_jaccard = (
            len(shared) / len(neighbor_union)
            if len(neighbor_union) > 0
            else np.nan
        )
        
        # Count high-discord neighbors
        high_discord_neighbors = sum(
            1 for n in atomx_neighbors
            if discord_map.get(n) == "high"
        )
        
        results.append({
            "atomx_cell": atomx_id,
            "proseg_cell": proseg_id,
            "n_neighbors_atomx": len(atomx_neighbors),
            "n_neighbors_proseg": len(proseg_neighbors),
            "shared_neighbors": len(shared),
            "gained_neighbors": len(gained),
            "lost_neighbors": len(lost),
            "neighbor_jaccard": neighbor_jaccard,
            "high_discord_neighbors": high_discord_neighbors,
            "neighbor_discord_density": (
                high_discord_neighbors / max(len(atomx_neighbors), 1)
            ),
        })
    
    result_df = pd.DataFrame(results) if results else pd.DataFrame()
    
    if logger and not result_df.empty:
        logger.info(f"Computed neighborhood metrics for {len(result_df)} cell pairs")
    
    return result_df


def compute_local_split_merge_density(
    cell_matches: pd.DataFrame,
    neighbors_atomx: Dict[str, List[str]],
    neighbors_proseg: Dict[str, List[str]],
    logger=None
) -> pd.DataFrame:
    """Compute local split/merge relationship density in neighborhoods.
    
    Parameters
    ----------
    cell_matches : pd.DataFrame
        Cell matching results.
    neighbors_atomx : dict
        AtoMx neighbor mappings.
    neighbors_proseg : dict
        Proseg neighbor mappings.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Local split/merge density per cell pair with columns:
        - atomx_cell, proseg_cell
        - n_split_related_neighbors
        - n_merge_related_neighbors
        - local_split_density
        - local_merge_density
    """
    if logger:
        logger.info("Computing local split/merge density")
    
    # Classify cell relationships
    splits = set()
    merges = set()
    
    # Find split cells (1 AtoMx -> multiple Proseg)
    for _, row in cell_matches.iterrows():
        atomx_id = row["atomx_cell"]
        matches_for_cell = cell_matches[cell_matches["atomx_cell"] == atomx_id]
        if len(matches_for_cell) > 1:
            splits.add(atomx_id)
    
    # Find merge cells (multiple AtoMx -> 1 Proseg)
    for _, row in cell_matches.iterrows():
        proseg_id = row["proseg_cell"]
        matches_for_cell = cell_matches[cell_matches["proseg_cell"] == proseg_id]
        if len(matches_for_cell) > 1:
            merges.add(proseg_id)
    
    results = []
    
    for _, match_row in cell_matches.iterrows():
        atomx_id = match_row["atomx_cell"]
        proseg_id = match_row["proseg_cell"]
        
        # Get neighbors
        atomx_neighbors = neighbors_atomx.get(atomx_id, [])
        proseg_neighbors = neighbors_proseg.get(proseg_id, [])
        
        # Count split-related neighbors (neighbors involved in splits)
        split_neighbors = sum(1 for n in atomx_neighbors if n in splits)
        
        # Count merge-related neighbors (neighbors involved in merges)
        merge_neighbors = sum(1 for n in proseg_neighbors if n in merges)
        
        results.append({
            "atomx_cell": atomx_id,
            "proseg_cell": proseg_id,
            "n_split_related_neighbors": split_neighbors,
            "n_merge_related_neighbors": merge_neighbors,
            "local_split_density": (
                split_neighbors / max(len(atomx_neighbors), 1)
            ),
            "local_merge_density": (
                merge_neighbors / max(len(proseg_neighbors), 1)
            ),
        })
    
    result_df = pd.DataFrame(results) if results else pd.DataFrame()
    
    if logger and not result_df.empty:
        logger.info(
            f"Computed local split/merge density for {len(result_df)} cell pairs"
        )
    
    return result_df


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
