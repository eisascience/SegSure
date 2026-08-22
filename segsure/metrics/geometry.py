"""Geometric disagreement metrics for SegSure."""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

try:
    from shapely.geometry import Polygon, MultiPolygon
    from shapely.ops import unary_union
    HAS_SHAPELY = True
except ImportError:
    HAS_SHAPELY = False


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
        AtoMx coordinate column names {'x': col, 'y': col}.
    proseg_coords : dict
        Proseg coordinate column names {'x': col, 'y': col}.
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


def polygon_from_vertices(
    vertices: List[Tuple[float, float]]
) -> Optional[Polygon]:
    """Create a Shapely polygon from vertices.
    
    Parameters
    ----------
    vertices : list of tuples
        Vertices as (x, y) coordinates.
    
    Returns
    -------
    Polygon or None
        Valid polygon if > 2 vertices, else None.
    """
    if not HAS_SHAPELY:
        return None
    
    if len(vertices) < 3:
        return None
    
    try:
        poly = Polygon(vertices)
        if poly.is_valid:
            return poly
        # Try to fix self-intersecting polygons
        return poly.buffer(0)
    except Exception:
        return None


def compute_polygon_metrics(
    atomx_polygon: Optional[Polygon],
    proseg_polygon: Optional[Polygon],
) -> Dict[str, Optional[float]]:
    """Compute geometry metrics between two polygons.
    
    Parameters
    ----------
    atomx_polygon : Polygon or None
        AtoMx cell polygon.
    proseg_polygon : Polygon or None
        Proseg cell polygon.
    
    Returns
    -------
    dict
        Dictionary with keys:
        - iou: Intersection over Union
        - coverage_atomx: fraction of AtoMx covered by Proseg
        - coverage_proseg: fraction of Proseg covered by AtoMx
        - area_difference: abs(atomx_area - proseg_area)
        - area_ratio: atomx_area / proseg_area
        - atomx_area: AtoMx polygon area
        - proseg_area: Proseg polygon area
    """
    metrics = {
        "iou": None,
        "coverage_atomx": None,
        "coverage_proseg": None,
        "area_difference": None,
        "area_ratio": None,
        "atomx_area": None,
        "proseg_area": None,
    }
    
    if not HAS_SHAPELY or atomx_polygon is None or proseg_polygon is None:
        return metrics
    
    try:
        # Compute intersection and union
        intersection = atomx_polygon.intersection(proseg_polygon)
        union = atomx_polygon.union(proseg_polygon)
        
        # IoU
        if union.area > 0:
            metrics["iou"] = intersection.area / union.area
        
        # Coverage fractions
        if atomx_polygon.area > 0:
            metrics["coverage_atomx"] = intersection.area / atomx_polygon.area
        if proseg_polygon.area > 0:
            metrics["coverage_proseg"] = intersection.area / proseg_polygon.area
        
        # Areas
        metrics["atomx_area"] = atomx_polygon.area
        metrics["proseg_area"] = proseg_polygon.area
        
        # Area metrics
        metrics["area_difference"] = abs(
            atomx_polygon.area - proseg_polygon.area
        )
        if proseg_polygon.area > 0:
            metrics["area_ratio"] = atomx_polygon.area / proseg_polygon.area
    
    except Exception:
        pass
    
    return metrics


def compute_polygon_overlap(
    atomx_polygons: Optional[pd.DataFrame],
    proseg_polygons: Optional[pd.DataFrame],
    cell_matches: pd.DataFrame,
    logger=None
) -> pd.DataFrame:
    """Compute polygon overlap metrics for matched cells.
    
    Parameters
    ----------
    atomx_polygons : pd.DataFrame, optional
        AtoMx polygon vertices with columns:
        - cell_id: cell identifier
        - x, y: vertex coordinates (one row per vertex)
    proseg_polygons : pd.DataFrame, optional
        Proseg polygon vertices (same schema).
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
        logger.info("Computing polygon overlap metrics")
    
    if not HAS_SHAPELY:
        if logger:
            logger.warning("Shapely not available; skipping polygon metrics")
        return pd.DataFrame()
    
    if atomx_polygons is None or proseg_polygons is None:
        if logger:
            logger.warning("Polygon data not provided; skipping polygon metrics")
        return pd.DataFrame()
    
    # Build polygon dict for each method
    atomx_poly_dict = _build_polygon_dict(atomx_polygons, logger)
    proseg_poly_dict = _build_polygon_dict(proseg_polygons, logger)
    
    results = []
    for _, match_row in cell_matches.iterrows():
        atomx_id = match_row["atomx_cell"]
        proseg_id = match_row["proseg_cell"]
        
        atomx_poly = atomx_poly_dict.get(atomx_id)
        proseg_poly = proseg_poly_dict.get(proseg_id)
        
        metrics = compute_polygon_metrics(atomx_poly, proseg_poly)
        metrics["atomx_cell"] = atomx_id
        metrics["proseg_cell"] = proseg_id
        
        results.append(metrics)
    
    if logger:
        logger.info(f"Computed polygon metrics for {len(results)} cell pairs")
    
    return pd.DataFrame(results) if results else pd.DataFrame()


def _build_polygon_dict(
    polygon_df: pd.DataFrame,
    logger=None
) -> Dict[str, Optional[Polygon]]:
    """Build dictionary mapping cell ID to Polygon.
    
    Parameters
    ----------
    polygon_df : pd.DataFrame
        Polygon vertex data.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    dict
        Mapping cell_id -> Polygon.
    """
    if not HAS_SHAPELY:
        return {}
    
    # Infer cell and coordinate columns
    cell_col = None
    for col in ["cell_id", "cellID", "cell", "nucleus_id"]:
        if col in polygon_df.columns:
            cell_col = col
            break
    
    x_col = None
    y_col = None
    for col in polygon_df.columns:
        if "x" in col.lower() and x_col is None:
            x_col = col
        elif "y" in col.lower() and y_col is None:
            y_col = col
    
    if cell_col is None or x_col is None or y_col is None:
        if logger:
            logger.warning("Could not infer cell/coordinate columns in polygon data")
        return {}
    
    poly_dict = {}
    
    for cell_id in polygon_df[cell_col].unique():
        cell_vertices = polygon_df[polygon_df[cell_col] == cell_id]
        vertices = list(zip(cell_vertices[x_col], cell_vertices[y_col]))
        
        poly = polygon_from_vertices(vertices)
        if poly is not None:
            poly_dict[cell_id] = poly
    
    if logger:
        logger.info(f"Built polygons for {len(poly_dict)} cells")
    
    return poly_dict


def compute_boundary_disagreement(
    atomx_cells: pd.DataFrame,
    proseg_cells: pd.DataFrame,
    atomx_polygons: Optional[pd.DataFrame],
    proseg_polygons: Optional[pd.DataFrame],
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
        Boundary disagreement metrics.
    """
    if logger:
        logger.info("Computing boundary disagreement metrics")
    
    if not HAS_SHAPELY or atomx_polygons is None or proseg_polygons is None:
        if logger:
            logger.warning(
                "Shapely or polygon data not available; "
                "skipping boundary metrics"
            )
        results = []
        for _, match_row in cell_matches.iterrows():
            results.append({
                "atomx_cell": match_row["atomx_cell"],
                "proseg_cell": match_row["proseg_cell"],
                "boundary_hausdorff": None,
                "mean_boundary_distance": None,
            })
        return pd.DataFrame(results) if results else pd.DataFrame()
    
    # Build polygon dict
    atomx_poly_dict = _build_polygon_dict(atomx_polygons, logger)
    proseg_poly_dict = _build_polygon_dict(proseg_polygons, logger)
    
    results = []
    for _, match_row in cell_matches.iterrows():
        atomx_id = match_row["atomx_cell"]
        proseg_id = match_row["proseg_cell"]
        
        atomx_poly = atomx_poly_dict.get(atomx_id)
        proseg_poly = proseg_poly_dict.get(proseg_id)
        
        boundary_metrics = _compute_boundary_metrics(atomx_poly, proseg_poly)
        boundary_metrics["atomx_cell"] = atomx_id
        boundary_metrics["proseg_cell"] = proseg_id
        
        results.append(boundary_metrics)
    
    if logger:
        logger.info(f"Computed boundary metrics for {len(results)} cell pairs")
    
    return pd.DataFrame(results) if results else pd.DataFrame()


def _compute_boundary_metrics(
    atomx_poly: Optional[Polygon],
    proseg_poly: Optional[Polygon]
) -> Dict[str, Optional[float]]:
    """Compute boundary distance metrics between two polygons.
    
    Parameters
    ----------
    atomx_poly : Polygon or None
        AtoMx boundary.
    proseg_poly : Polygon or None
        Proseg boundary.
    
    Returns
    -------
    dict
        Dictionary with keys:
        - boundary_hausdorff: Hausdorff distance between boundaries
        - mean_boundary_distance: Mean distance between boundaries
    """
    metrics = {
        "boundary_hausdorff": None,
        "mean_boundary_distance": None,
    }
    
    if not HAS_SHAPELY or atomx_poly is None or proseg_poly is None:
        return metrics
    
    try:
        # Hausdorff distance
        metrics["boundary_hausdorff"] = atomx_poly.boundary.hausdorff_distance(
            proseg_poly.boundary
        )
        
        # Mean distance (approximate)
        # Sample points on boundary and compute mean distance
        atomx_boundary_coords = np.array(atomx_poly.boundary.coords)
        proseg_boundary_coords = np.array(proseg_poly.boundary.coords)
        
        if len(atomx_boundary_coords) > 0 and len(proseg_boundary_coords) > 0:
            distances = cdist(
                atomx_boundary_coords[:, :2],
                proseg_boundary_coords[:, :2],
                metric="euclidean"
            )
            metrics["mean_boundary_distance"] = np.mean(np.min(distances, axis=1))
    
    except Exception:
        pass
    
    return metrics
