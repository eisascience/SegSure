"""Geometric disagreement metrics for SegSure."""

from typing import Dict, Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


def compute_centroid_distance(
    atomx_cells: pd.DataFrame,
    proseg_cells: pd.DataFrame,
    atomx_coords: Dict[str, str],
    proseg_coords: Dict[str, str],
    cell_matches: pd.DataFrame,
    logger=None
) -> pd.DataFrame:
    """Compute centroid distances for matched cell pairs.
    
    Parameters
    ----------
    atomx_cells : pd.DataFrame
        AtoMx cell data with centroids.
    proseg_cells : pd.DataFrame
        Proseg cell data with centroids.
    atomx_coords : dict
        AtoMx coordinate column names.
    proseg_coords : dict
        Proseg coordinate column names.
    cell_matches : pd.DataFrame
        Cell matching results.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Centroid distances for each match.
    """
    if logger:
        logger.info("Computing centroid distances")
    
    results = []
    
    for _, match_row in cell_matches.iterrows():
        atomx_id = match_row["atomx_cell"]
        proseg_id = match_row["proseg_cell"]
        
        atomx_row = atomx_cells[atomx_cells.index == atomx_id]
        proseg_row = proseg_cells[proseg_cells.index == proseg_id]
        
        if len(atomx_row) > 0 and len(proseg_row) > 0:
            ax = atomx_row[atomx_coords["x"]].values[0]
            ay = atomx_row[atomx_coords["y"]].values[0]
            px = proseg_row[proseg_coords["x"]].values[0]
            py = proseg_row[proseg_coords["y"]].values[0]
            
            distance = ((ax - px) ** 2 + (ay - py) ** 2) ** 0.5
            
            results.append({
                "atomx_cell": atomx_id,
                "proseg_cell": proseg_id,
                "centroid_distance": distance,
            })
    
    return pd.DataFrame(results) if results else pd.DataFrame()


def compute_polygon_overlap(
    atomx_polygons: Optional[pd.DataFrame],
    proseg_polygons: Optional[pd.DataFrame],
    cell_matches: pd.DataFrame,
    logger=None
) -> pd.DataFrame:
    """Compute polygon overlap (Intersection over Union) for matched cells.
    
    Parameters
    ----------
    atomx_polygons : pd.DataFrame, optional
        AtoMx polygon vertices.
    proseg_polygons : pd.DataFrame, optional
        Proseg polygon vertices.
    cell_matches : pd.DataFrame
        Cell matching results.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Polygon overlap metrics for each match.
    """
    if logger:
        logger.info("Computing polygon overlap")
    
    # Note: Full polygon overlap computation requires shapely and
    # properly structured polygon vertex data.
    # This is a placeholder that returns the structure.
    
    results = []
    for _, match_row in cell_matches.iterrows():
        results.append({
            "atomx_cell": match_row["atomx_cell"],
            "proseg_cell": match_row["proseg_cell"],
            "iou": None,  # Requires polygon geometry
            "overlap_area": None,
            "union_area": None,
        })
    
    if logger:
        logger.warning("Polygon overlap computation requires polygon geometry data")
    
    return pd.DataFrame(results) if results else pd.DataFrame()


def compute_boundary_disagreement(
    atomx_cells: pd.DataFrame,
    proseg_cells: pd.DataFrame,
    cell_matches: pd.DataFrame,
    logger=None
) -> pd.DataFrame:
    """Compute boundary disagreement metrics.
    
    Parameters
    ----------
    atomx_cells : pd.DataFrame
        AtoMx cell data.
    proseg_cells : pd.DataFrame
        Proseg cell data.
    cell_matches : pd.DataFrame
        Cell matching results.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Boundary disagreement metrics.
    """
    if logger:
        logger.info("Computing boundary disagreement")
    
    results = []
    for _, match_row in cell_matches.iterrows():
        results.append({
            "atomx_cell": match_row["atomx_cell"],
            "proseg_cell": match_row["proseg_cell"],
            "hausdorff_distance": None,  # Requires polygon boundaries
            "mean_boundary_distance": None,
        })
    
    return pd.DataFrame(results) if results else pd.DataFrame()
