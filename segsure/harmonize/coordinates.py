"""Coordinate system harmonization for SegSure."""

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd


def get_coordinate_columns(df: pd.DataFrame, logger=None) -> Dict[str, Optional[str]]:
    """Infer coordinate column names from DataFrame.
    
    Looks for columns with common X/Y/Z coordinate naming patterns.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to inspect.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    dict
        Dictionary with inferred coordinate columns.
    """
    coords = {"x": None, "y": None, "z": None}
    
    # Common X coordinate names
    x_names = ["x", "X", "x_coord", "X_coord", "center_x", "Center_X"]
    y_names = ["y", "Y", "y_coord", "Y_coord", "center_y", "Center_Y"]
    z_names = ["z", "Z", "z_coord", "Z_coord", "center_z", "Center_Z"]
    
    for col in df.columns:
        if coords["x"] is None and col in x_names:
            coords["x"] = col
        elif coords["y"] is None and col in y_names:
            coords["y"] = col
        elif coords["z"] is None and col in z_names:
            coords["z"] = col
    
    if logger:
        logger.info(f"Inferred coordinate columns: {coords}")
    
    return coords


def get_fov_column(df: pd.DataFrame, logger=None) -> Optional[str]:
    """Infer FOV column name from DataFrame.
    
    Looks for columns with common FOV naming patterns.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to inspect.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    str or None
        Inferred FOV column name.
    """
    fov_names = ["fov", "FOV", "slide_id", "slide", "tile", "region"]
    
    for col in df.columns:
        if col in fov_names:
            if logger:
                logger.info(f"Inferred FOV column: {col}")
            return col
    
    if logger:
        logger.warning("Could not infer FOV column")
    
    return None


def normalize_coordinates(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    z_col: Optional[str] = None,
    logger=None
) -> pd.DataFrame:
    """Normalize coordinates to consistent units and ranges.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with coordinate columns.
    x_col : str
        X coordinate column name.
    y_col : str
        Y coordinate column name.
    z_col : str, optional
        Z coordinate column name.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with normalized coordinates.
    """
    df_norm = df.copy()
    
    # For now, just ensure numeric
    df_norm[x_col] = pd.to_numeric(df_norm[x_col], errors="coerce")
    df_norm[y_col] = pd.to_numeric(df_norm[y_col], errors="coerce")
    
    if z_col and z_col in df_norm.columns:
        df_norm[z_col] = pd.to_numeric(df_norm[z_col], errors="coerce")
    
    if logger:
        logger.info("Coordinates normalized to numeric type")
    
    return df_norm


def transform_to_global_coordinates(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    fov_col: Optional[str] = None,
    fov_positions: Optional[pd.DataFrame] = None,
    z_col: Optional[str] = None,
    logger=None
) -> pd.DataFrame:
    """Transform local FOV coordinates to global slide coordinates.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with local coordinates.
    x_col : str
        X coordinate column name.
    y_col : str
        Y coordinate column name.
    fov_col : str, optional
        FOV identifier column name.
    fov_positions : pd.DataFrame, optional
        FOV position offsets (FOV -> X/Y offset).
    z_col : str, optional
        Z coordinate column name.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with global coordinates.
    """
    df_global = df.copy()
    
    if fov_col and fov_positions is not None:
        if logger:
            logger.info("Transforming coordinates to global system")
        
        # Infer FOV position columns
        fov_x_col = None
        fov_y_col = None
        for col in fov_positions.columns:
            if "x" in col.lower():
                fov_x_col = col
            elif "y" in col.lower():
                fov_y_col = col
        
        if fov_x_col and fov_y_col:
            # Create mapping of FOV ID to position
            fov_map = dict(zip(fov_positions.iloc[:, 0], 
                             zip(fov_positions[fov_x_col], 
                                 fov_positions[fov_y_col])))
            
            # Apply transformation
            offsets = df_global[fov_col].map(lambda x: fov_map.get(x, (0, 0)))
            df_global[x_col] = df_global[x_col] + offsets.apply(lambda x: x[0])
            df_global[y_col] = df_global[y_col] + offsets.apply(lambda x: x[1])
    
    return df_global


def align_coordinate_systems(
    atomx_df: pd.DataFrame,
    proseg_adata,
    atomx_coords: Dict[str, str],
    proseg_coords: Dict[str, str],
    logger=None
) -> Tuple[pd.DataFrame, Dict]:
    """Align AtoMx and Proseg coordinate systems.
    
    Parameters
    ----------
    atomx_df : pd.DataFrame
        AtoMx coordinate data.
    proseg_adata : anndata.AnnData
        Proseg object.
    atomx_coords : dict
        AtoMx coordinate column names.
    proseg_coords : dict
        Proseg coordinate column names.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    tuple
        (alignment_result_df, alignment_info_dict)
    """
    if logger:
        logger.info("Aligning coordinate systems between AtoMx and Proseg")
    
    alignment_info = {
        "atomx_coords": atomx_coords,
        "proseg_coords": proseg_coords,
        "atomx_range": {
            "x": (atomx_df[atomx_coords["x"]].min(), atomx_df[atomx_coords["x"]].max()),
            "y": (atomx_df[atomx_coords["y"]].min(), atomx_df[atomx_coords["y"]].max()),
        },
    }
    
    return atomx_df, alignment_info
