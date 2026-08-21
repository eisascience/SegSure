"""Cell identification and matching for SegSure."""

from typing import Dict, List, Optional, Set, Tuple

import pandas as pd


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


def match_cells(
    atomx_cells: pd.DataFrame,
    proseg_cells: pd.DataFrame,
    atomx_coords: Dict[str, str],
    proseg_coords: Dict[str, str],
    atomx_cell_id: str,
    proseg_cell_id: str,
    distance_threshold: float = 10.0,
    logger=None
) -> pd.DataFrame:
    """Match cells between AtoMx and Proseg based on spatial proximity.
    
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
    distance_threshold : float, optional
        Maximum distance for cell matching. Default is 10.0.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Matching results with atomx_cell, proseg_cell, distance.
    """
    if logger:
        logger.info(
            f"Matching {len(atomx_cells)} AtoMx cells to {len(proseg_cells)} "
            f"Proseg cells (threshold={distance_threshold})"
        )
    
    matches = []
    
    for _, atomx_row in atomx_cells.iterrows():
        atomx_x = atomx_row[atomx_coords["x"]]
        atomx_y = atomx_row[atomx_coords["y"]]
        atomx_id = atomx_row[atomx_cell_id]
        
        # Calculate distances to all Proseg cells
        distances = []
        for _, proseg_row in proseg_cells.iterrows():
            proseg_x = proseg_row[proseg_coords["x"]]
            proseg_y = proseg_row[proseg_coords["y"]]
            proseg_id = proseg_row[proseg_cell_id]
            
            dist = ((atomx_x - proseg_x) ** 2 + (atomx_y - proseg_y) ** 2) ** 0.5
            
            if dist <= distance_threshold:
                distances.append({
                    "atomx_cell": atomx_id,
                    "proseg_cell": proseg_id,
                    "distance": dist,
                })
        
        # Keep best match if exists
        if distances:
            best_match = min(distances, key=lambda x: x["distance"])
            matches.append(best_match)
    
    if logger:
        logger.info(f"Found {len(matches)} cell matches")
    
    return pd.DataFrame(matches) if matches else pd.DataFrame()


def identify_split_merges(
    matches: pd.DataFrame,
    atomx_cell_id: str = "atomx_cell",
    proseg_cell_id: str = "proseg_cell",
    logger=None
) -> Dict[str, List]:
    """Identify split and merge events from cell matches.
    
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
    atomx_cells = set(matches[atomx_cell_id].unique())
    proseg_cells = set(matches[proseg_cell_id].unique())
    
    # Identify relationship types
    one_to_one = []
    splits = []
    merges = []
    
    # Count matches per cell
    atomx_match_counts = matches[atomx_cell_id].value_counts()
    proseg_match_counts = matches[proseg_cell_id].value_counts()
    
    for atomx_id in atomx_match_counts.index:
        n_matches = atomx_match_counts[atomx_id]
        proseg_ids = matches[matches[atomx_cell_id] == atomx_id][proseg_cell_id].tolist()
        
        if n_matches == 1:
            # Check if the Proseg cell also has only one match
            proseg_id = proseg_ids[0]
            if proseg_match_counts[proseg_id] == 1:
                one_to_one.append({"atomx_cell": atomx_id, "proseg_cell": proseg_id})
            else:
                # Merge event (one AtoMx mapped to Proseg with multiple matches)
                merges.append({"atomx_cell": atomx_id, "proseg_cells": proseg_ids})
        else:
            # Split event
            splits.append({"atomx_cell": atomx_id, "proseg_cells": proseg_ids})
    
    # Identify lost and newly inferred cells
    all_atomx = set(matches[atomx_cell_id]) if len(matches) > 0 else set()
    all_proseg = set(matches[proseg_cell_id]) if len(matches) > 0 else set()
    
    lost = list(atomx_cells - all_atomx) if atomx_cells else []
    newly_inferred = list(proseg_cells - all_proseg) if proseg_cells else []
    
    result = {
        "one_to_one": one_to_one,
        "splits": splits,
        "merges": merges,
        "lost_cells": lost,
        "newly_inferred_cells": newly_inferred,
    }
    
    if logger:
        logger.info(f"Cell relationship summary:")
        logger.info(f"  One-to-one: {len(one_to_one)}")
        logger.info(f"  Splits: {len(splits)}")
        logger.info(f"  Merges: {len(merges)}")
        logger.info(f"  Lost cells: {len(lost)}")
        logger.info(f"  Newly inferred cells: {len(newly_inferred)}")
    
    return result
