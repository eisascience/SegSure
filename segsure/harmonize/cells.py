"""Cell identification and matching for SegSure using polygon-based methods."""

from typing import Dict, List, Optional, Set, Tuple, Any
import warnings

import pandas as pd
import numpy as np
from scipy.spatial import distance
from shapely.geometry import shape, Polygon
from shapely.strtree import STRtree
import networkx as nx

# Configurable thresholds (documented)
DEFAULT_THRESHOLDS = {
    # Minimum IoU for one-to-one match
    "iou_threshold_one_to_one": 0.5,
    # Minimum coverage for considering an overlap
    "min_coverage": 0.01,
    # Maximum distance for centroid-based fallback
    "max_centroid_distance": 100.0,
    # Minimum overlap area
    "min_overlap_area": 1.0,
}


def get_cell_id_column(df: pd.DataFrame, logger=None) -> Optional[str]:
    """Infer cell ID column name from DataFrame.
    
    Looks for columns with common cell ID naming patterns.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to inspect.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    str or None
        Inferred cell ID column name.
    """
    cell_names = ["cell_id", "cellID", "cell_index", "cluster", "nucleus_id", "segmentation_id"]
    
    for col in df.columns:
        if col.lower() in [n.lower() for n in cell_names]:
            if logger:
                logger.info(f"Inferred cell ID column: {col}")
            return col
    
    # Check index
    if df.index.name:
        if logger:
            logger.info(f"Using index as cell ID: {df.index.name}")
        return df.index.name
    
    if logger:
        logger.warning("Could not infer cell ID column")
    
    return None


def normalize_cell_ids(
    df: pd.DataFrame,
    cell_id_col: str,
    method: str = "pass_through",
    prefix: Optional[str] = None,
    logger=None
) -> pd.DataFrame:
    """Normalize cell ID values for consistent matching.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with cell IDs.
    cell_id_col : str
        Cell ID column name.
    method : str, optional
        Normalization method ("pass_through", "string", "hash").
    prefix : str, optional
        Prefix to add to cell IDs.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with normalized cell IDs.
    """
    df_norm = df.copy()
    
    if method == "string":
        df_norm[cell_id_col] = df_norm[cell_id_col].astype(str)
    
    if prefix:
        df_norm[cell_id_col] = prefix + df_norm[cell_id_col].astype(str)
    
    if logger:
        logger.info(f"Normalized cell IDs using method={method}, prefix={prefix}")
    
    return df_norm


def parse_polygons(
    polygon_df: pd.DataFrame,
    cell_id_col: str = "cell_id",
    geometry_col: str = "geometry",
    logger=None
) -> Dict[str, Optional[Polygon]]:
    """Parse polygon data from DataFrame.
    
    Assumes geometry column contains WKT or similar format.
    
    Parameters
    ----------
    polygon_df : pd.DataFrame
        DataFrame with polygon data.
    cell_id_col : str
        Cell ID column name.
    geometry_col : str
        Geometry column name (WKT format).
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    dict
        Mapping of cell IDs to Polygon objects.
    """
    polygons = {}
    
    if geometry_col not in polygon_df.columns:
        if logger:
            logger.warning(f"Geometry column '{geometry_col}' not found in polygon DataFrame")
        return polygons
    
    for idx, row in polygon_df.iterrows():
        cell_id = row[cell_id_col]
        geom_data = row[geometry_col]
        
        try:
            # Try parsing as WKT
            if isinstance(geom_data, str):
                from shapely.wkt import loads as wkt_loads
                poly = wkt_loads(geom_data)
            else:
                # Try parsing as GeoJSON-like dict
                poly = shape(geom_data)
            
            if poly.is_valid and poly.geom_type == "Polygon":
                polygons[cell_id] = poly
            else:
                if logger:
                    logger.warning(f"Invalid polygon for cell {cell_id}")
        except Exception as e:
            if logger:
                logger.warning(f"Failed to parse polygon for cell {cell_id}: {e}")
    
    if logger:
        logger.info(f"Parsed {len(polygons)} polygons from {len(polygon_df)} rows")
    
    return polygons


def calculate_overlap_metrics(
    poly_a: Polygon,
    poly_b: Polygon
) -> Dict[str, float]:
    """Calculate overlap metrics between two polygons.
    
    Parameters
    ----------
    poly_a : Polygon
        First polygon.
    poly_b : Polygon
        Second polygon.
    
    Returns
    -------
    dict
        Dictionary with overlap metrics:
        - intersection_area: Area of intersection
        - union_area: Area of union
        - iou: Intersection over union
        - coverage_a: Fraction of polygon A covered by B
        - coverage_b: Fraction of polygon B covered by A
        - centroid_distance: Distance between centroids
        - area_ratio: Ratio of areas (A/B)
    """
    try:
        intersection = poly_a.intersection(poly_b)
        intersection_area = intersection.area
        
        union = poly_a.union(poly_b)
        union_area = union.area
        
        # Avoid division by zero
        iou = intersection_area / union_area if union_area > 0 else 0.0
        
        coverage_a = intersection_area / poly_a.area if poly_a.area > 0 else 0.0
        coverage_b = intersection_area / poly_b.area if poly_b.area > 0 else 0.0
        
        centroid_a = poly_a.centroid
        centroid_b = poly_b.centroid
        centroid_distance = centroid_a.distance(centroid_b)
        
        area_ratio = poly_a.area / poly_b.area if poly_b.area > 0 else 0.0
        
        return {
            "intersection_area": float(intersection_area),
            "union_area": float(union_area),
            "iou": float(iou),
            "coverage_a": float(coverage_a),
            "coverage_b": float(coverage_b),
            "centroid_distance": float(centroid_distance),
            "area_ratio": float(area_ratio),
        }
    except Exception as e:
        warnings.warn(f"Error calculating overlap metrics: {e}")
        return {
            "intersection_area": 0.0,
            "union_area": 0.0,
            "iou": 0.0,
            "coverage_a": 0.0,
            "coverage_b": 0.0,
            "centroid_distance": float("inf"),
            "area_ratio": 0.0,
        }


def match_cells_polygon_based(
    atomx_cells: pd.DataFrame,
    proseg_cells: pd.DataFrame,
    atomx_polygons: Dict[str, Polygon],
    proseg_polygons: Dict[str, Polygon],
    atomx_coords: Dict[str, str],
    proseg_coords: Dict[str, str],
    atomx_cell_id: str,
    proseg_cell_id: str,
    thresholds: Optional[Dict[str, float]] = None,
    logger=None
) -> pd.DataFrame:
    """Match cells using polygon overlap with spatial indexing.
    
    Parameters
    ----------
    atomx_cells : pd.DataFrame
        AtoMx cell data.
    proseg_cells : pd.DataFrame
        Proseg cell data.
    atomx_polygons : dict
        Mapping of AtoMx cell IDs to Polygon objects.
    proseg_polygons : dict
        Mapping of Proseg cell IDs to Polygon objects.
    atomx_coords : dict
        AtoMx coordinate column names.
    proseg_coords : dict
        Proseg coordinate column names.
    atomx_cell_id : str
        AtoMx cell ID column name.
    proseg_cell_id : str
        Proseg cell ID column name.
    thresholds : dict, optional
        Matching thresholds.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Overlap metrics for all candidate pairs.
    """
    if thresholds is None:
        thresholds = DEFAULT_THRESHOLDS
    
    if logger:
        logger.info(f"Matching {len(atomx_cells)} AtoMx to {len(proseg_cells)} Proseg cells")
        logger.info(f"AtoMx polygons available: {len(atomx_polygons)}")
        logger.info(f"Proseg polygons available: {len(proseg_polygons)}")
    
    # Fall back to centroid-based matching if polygons unavailable
    if not atomx_polygons or not proseg_polygons:
        if logger:
            logger.warning("Using centroid-based fallback matching")
        return match_cells_centroid_based(
            atomx_cells, proseg_cells, atomx_coords, proseg_coords,
            atomx_cell_id, proseg_cell_id, thresholds, logger
        )
    
    # Build spatial index for Proseg polygons
    proseg_ids = list(proseg_polygons.keys())
    proseg_geoms = [proseg_polygons[cid] for cid in proseg_ids]
    tree = STRtree(proseg_geoms)
    
    candidate_pairs = []
    
    # Find candidate pairs using spatial index
    for atomx_id in atomx_cells[atomx_cell_id].unique():
        if atomx_id not in atomx_polygons:
            continue
        
        atomx_poly = atomx_polygons[atomx_id]
        
        # Use envelope for broad query
        envelope = atomx_poly.envelope
        candidates_indices = tree.query(envelope)
        
        # For each candidate, calculate metrics
        for cand_idx in candidates_indices:
            proseg_id = proseg_ids[cand_idx]
            proseg_poly = proseg_polygons[proseg_id]
            
            metrics = calculate_overlap_metrics(atomx_poly, proseg_poly)
            
            # Only keep pairs with meaningful overlap
            if metrics["intersection_area"] >= thresholds["min_overlap_area"]:
                candidate_pairs.append({
                    "atomx_cell": atomx_id,
                    "proseg_cell": proseg_id,
                    **metrics,
                })
    
    df_candidates = pd.DataFrame(candidate_pairs) if candidate_pairs else pd.DataFrame()
    
    if logger:
        logger.info(f"Found {len(df_candidates)} candidate pairs")
    
    return df_candidates


def match_cells_centroid_based(
    atomx_cells: pd.DataFrame,
    proseg_cells: pd.DataFrame,
    atomx_coords: Dict[str, str],
    proseg_coords: Dict[str, str],
    atomx_cell_id: str,
    proseg_cell_id: str,
    thresholds: Optional[Dict[str, float]] = None,
    logger=None
) -> pd.DataFrame:
    """Fallback centroid-based matching when polygons unavailable.
    
    Parameters
    ----------
    atomx_cells : pd.DataFrame
        AtoMx cell data with coordinates.
    proseg_cells : pd.DataFrame
        Proseg cell data with coordinates.
    atomx_coords : dict
        AtoMx coordinate column names.
    proseg_coords : dict
        Proseg coordinate column names.
    atomx_cell_id : str
        AtoMx cell ID column name.
    proseg_cell_id : str
        Proseg cell ID column name.
    thresholds : dict, optional
        Matching thresholds.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Distance metrics for all candidate pairs.
    """
    if thresholds is None:
        thresholds = DEFAULT_THRESHOLDS
    
    if logger:
        logger.info("Using centroid-based fallback matching")
    
    distance_threshold = thresholds.get("max_centroid_distance", 100.0)
    
    pairs = []
    
    for _, atomx_row in atomx_cells.iterrows():
        atomx_id = atomx_row[atomx_cell_id]
        atomx_x = atomx_row[atomx_coords["x"]]
        atomx_y = atomx_row[atomx_coords["y"]]
        
        for _, proseg_row in proseg_cells.iterrows():
            proseg_id = proseg_row[proseg_cell_id]
            proseg_x = proseg_row[proseg_coords["x"]]
            proseg_y = proseg_row[proseg_coords["y"]]
            
            dist = np.sqrt((atomx_x - proseg_x) ** 2 + (atomx_y - proseg_y) ** 2)
            
            if dist <= distance_threshold:
                pairs.append({
                    "atomx_cell": atomx_id,
                    "proseg_cell": proseg_id,
                    "intersection_area": 0.0,
                    "union_area": 0.0,
                    "iou": 0.0,
                    "coverage_a": 0.0,
                    "coverage_b": 0.0,
                    "centroid_distance": float(dist),
                    "area_ratio": 0.0,
                })
    
    df_pairs = pd.DataFrame(pairs) if pairs else pd.DataFrame()
    
    if logger:
        logger.info(f"Found {len(df_pairs)} candidate pairs via centroid matching")
    
    return df_pairs


def classify_relationships(
    candidates: pd.DataFrame,
    atomx_all: Set[str],
    proseg_all: Set[str],
    atomx_cell_col: str = "atomx_cell",
    proseg_cell_col: str = "proseg_cell",
    iou_threshold: float = 0.5,
    logger=None
) -> pd.DataFrame:
    """Classify cell relationships from candidate pairs.
    
    Parameters
    ----------
    candidates : pd.DataFrame
        Candidate pairs with overlap metrics.
    atomx_all : set
        All AtoMx cell IDs in universe.
    proseg_all : set
        All Proseg cell IDs in universe.
    atomx_cell_col : str
        AtoMx cell ID column name.
    proseg_cell_col : str
        Proseg cell ID column name.
    iou_threshold : float
        Minimum IoU for considering an overlap.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Classification results with relationship types.
    """
    # Filter to meaningful overlaps
    df = candidates.copy()
    
    if len(df) == 0:
        df["relationship_type"] = []
        return df
    
    # Find the actual column names (they're always "atomx_cell" and "proseg_cell" from match_cells_polygon_based)
    if "atomx_cell" not in df.columns:
        # Try to find them
        atom_col = next((c for c in df.columns if "atomx" in c.lower()), "atomx_cell")
        proseg_col = next((c for c in df.columns if "proseg" in c.lower()), "proseg_cell")
    else:
        atom_col = "atomx_cell"
        proseg_col = "proseg_cell"
    
    # Override with passed parameters if they exist in the dataframe
    if atomx_cell_col in df.columns:
        atom_col = atomx_cell_col
    if proseg_cell_col in df.columns:
        proseg_col = proseg_cell_col
    
    # Build bipartite graph
    G = nx.Graph()
    
    # Add nodes
    for cell_id in atomx_all:
        G.add_node(f"atomx_{cell_id}", source="atomx")
    for cell_id in proseg_all:
        G.add_node(f"proseg_{cell_id}", source="proseg")
    
    # Add edges from candidates with meaningful overlap
    for _, row in df.iterrows():
        if row.get("iou", 0) >= iou_threshold or row.get("intersection_area", 0) > 0:
            G.add_edge(
                f"atomx_{row[atom_col]}",
                f"proseg_{row[proseg_col]}",
                iou=row.get("iou", 0),
                intersection_area=row.get("intersection_area", 0),
            )
    
    # Classify relationships
    relationships = []
    processed_atomx = set()
    processed_proseg = set()
    
    for _, row in df.iterrows():
        atomx_id = row[atom_col]
        proseg_id = row[proseg_col]
        
        atomx_node = f"atomx_{atomx_id}"
        proseg_node = f"proseg_{proseg_id}"
        
        # Get neighbors in graph
        atomx_neighbors = set(G.neighbors(atomx_node)) if atomx_node in G else set()
        proseg_neighbors = set(G.neighbors(proseg_node)) if proseg_node in G else set()
        
        # Count neighbors (excluding current pairing)
        atomx_neighbor_count = len([n for n in atomx_neighbors if n.startswith("proseg_")])
        proseg_neighbor_count = len([n for n in proseg_neighbors if n.startswith("atomx_")])
        
        # Classify
        if atomx_neighbor_count == 1 and proseg_neighbor_count == 1:
            relationship_type = "one_to_one"
        elif atomx_neighbor_count > 1 and proseg_neighbor_count == 1:
            relationship_type = "split"
        elif atomx_neighbor_count == 1 and proseg_neighbor_count > 1:
            relationship_type = "merge"
        else:
            relationship_type = "complex"
        
        row_dict = row.to_dict()
        row_dict["relationship_type"] = relationship_type
        relationships.append(row_dict)
        
        processed_atomx.add(atomx_id)
        processed_proseg.add(proseg_id)
    
    df_results = pd.DataFrame(relationships)
    
    # Add lost and new cells
    lost_atomx = atomx_all - processed_atomx
    new_proseg = proseg_all - processed_proseg
    
    # Create rows for lost cells
    for cell_id in lost_atomx:
        df_results = pd.concat([
            df_results,
            pd.DataFrame([{
                atom_col: cell_id,
                proseg_col: None,
                "relationship_type": "lost_atomx",
                "iou": 0.0,
                "intersection_area": 0.0,
                "union_area": 0.0,
                "coverage_a": 0.0,
                "coverage_b": 0.0,
                "centroid_distance": np.nan,
                "area_ratio": 0.0,
            }]),
        ], ignore_index=True)
    
    # Create rows for new cells
    for cell_id in new_proseg:
        df_results = pd.concat([
            df_results,
            pd.DataFrame([{
                atom_col: None,
                proseg_col: cell_id,
                "relationship_type": "new_proseg",
                "iou": 0.0,
                "intersection_area": 0.0,
                "union_area": 0.0,
                "coverage_a": 0.0,
                "coverage_b": 0.0,
                "centroid_distance": np.nan,
                "area_ratio": 0.0,
            }]),
        ], ignore_index=True)
    
    if logger:
        relationship_counts = df_results["relationship_type"].value_counts()
        logger.info("Relationship type counts:")
        for rel_type, count in relationship_counts.items():
            logger.info(f"  {rel_type}: {count}")
    
    return df_results


def match_cells(
    atomx_cells: pd.DataFrame,
    proseg_cells: pd.DataFrame,
    atomx_coords: Dict[str, str],
    proseg_coords: Dict[str, str],
    atomx_cell_id: str,
    proseg_cell_id: str,
    atomx_polygons: Optional[Dict[str, Polygon]] = None,
    proseg_polygons: Optional[Dict[str, Polygon]] = None,
    thresholds: Optional[Dict[str, float]] = None,
    distance_threshold: float = 10.0,
    logger=None
) -> pd.DataFrame:
    """Match cells between AtoMx and Proseg using polygon or centroid methods.
    
    Parameters
    ----------
    atomx_cells : pd.DataFrame
        AtoMx cell data with coordinates.
    proseg_cells : pd.DataFrame
        Proseg cell data with coordinates.
    atomx_coords : dict
        AtoMx coordinate column names.
    proseg_coords : dict
        Proseg coordinate column names.
    atomx_cell_id : str
        AtoMx cell ID column name.
    proseg_cell_id : str
        Proseg cell ID column name.
    atomx_polygons : dict, optional
        AtoMx cell polygons.
    proseg_polygons : dict, optional
        Proseg cell polygons.
    thresholds : dict, optional
        Matching thresholds.
    distance_threshold : float, optional
        Fallback distance threshold.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Matching results with classification.
    """
    if thresholds is None:
        thresholds = DEFAULT_THRESHOLDS.copy()
        thresholds["max_centroid_distance"] = distance_threshold
    
    # Get candidate pairs
    if atomx_polygons and proseg_polygons:
        candidates = match_cells_polygon_based(
            atomx_cells, proseg_cells, atomx_polygons, proseg_polygons,
            atomx_coords, proseg_coords, atomx_cell_id, proseg_cell_id,
            thresholds, logger
        )
    else:
        candidates = match_cells_centroid_based(
            atomx_cells, proseg_cells, atomx_coords, proseg_coords,
            atomx_cell_id, proseg_cell_id, thresholds, logger
        )
    
    # Classify relationships
    atomx_all = set(atomx_cells[atomx_cell_id].unique())
    proseg_all = set(proseg_cells[proseg_cell_id].unique())
    
    results = classify_relationships(
        candidates, atomx_all, proseg_all,
        atomx_cell_col=atomx_cell_id,
        proseg_cell_col=proseg_cell_id,
        iou_threshold=thresholds.get("iou_threshold_one_to_one", 0.5),
        logger=logger
    )
    
    return results


def identify_split_merges(
    matches: pd.DataFrame,
    atomx_cell_id: str = "atomx_cell",
    proseg_cell_id: str = "proseg_cell",
    logger=None
) -> Dict[str, List]:
    """Identify split and merge events from cell matches (legacy function).
    
    Parameters
    ----------
    matches : pd.DataFrame
        Cell matching results.
    atomx_cell_id : str, optional
        Column name for AtoMx cell IDs.
    proseg_cell_id : str, optional
        Column name for Proseg cell IDs.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    dict
        Dictionary with splits, merges, one_to_one, lost, and newly_inferred.
    """
    if "relationship_type" in matches.columns:
        # New format
        relationship_types = matches["relationship_type"].value_counts()
        
        one_to_one = matches[matches["relationship_type"] == "one_to_one"].to_dict("records")
        splits = matches[matches["relationship_type"] == "split"].to_dict("records")
        merges = matches[matches["relationship_type"] == "merge"].to_dict("records")
        lost = matches[matches["relationship_type"] == "lost_atomx"].to_dict("records")
        newly_inferred = matches[matches["relationship_type"] == "new_proseg"].to_dict("records")
    else:
        # Old format fallback
        atomx_cells = set(matches[atomx_cell_id].unique()) if len(matches) > 0 else set()
        proseg_cells = set(matches[proseg_cell_id].unique()) if len(matches) > 0 else set()
        
        one_to_one = []
        splits = []
        merges = []
        
        atomx_match_counts = matches[atomx_cell_id].value_counts()
        proseg_match_counts = matches[proseg_cell_id].value_counts()
        
        for atomx_id in atomx_match_counts.index:
            n_matches = atomx_match_counts[atomx_id]
            proseg_ids = matches[matches[atomx_cell_id] == atomx_id][proseg_cell_id].tolist()
            
            if n_matches == 1:
                proseg_id = proseg_ids[0]
                if proseg_match_counts[proseg_id] == 1:
                    one_to_one.append({"atomx_cell": atomx_id, "proseg_cell": proseg_id})
                else:
                    merges.append({"atomx_cell": atomx_id, "proseg_cells": proseg_ids})
            else:
                splits.append({"atomx_cell": atomx_id, "proseg_cells": proseg_ids})
        
        lost = list(atomx_cells - set(matches[atomx_cell_id]))
        newly_inferred = list(proseg_cells - set(matches[proseg_cell_id]))
    
    result = {
        "one_to_one": one_to_one,
        "splits": splits,
        "merges": merges,
        "lost_cells": lost,
        "newly_inferred_cells": newly_inferred,
    }
    
    if logger:
        logger.info("Cell relationship summary:")
        logger.info(f"  One-to-one: {len(one_to_one)}")
        logger.info(f"  Splits: {len(splits)}")
        logger.info(f"  Merges: {len(merges)}")
        logger.info(f"  Lost cells: {len(lost)}")
        logger.info(f"  Newly inferred cells: {len(newly_inferred)}")
    
    return result
