"""Validation utilities for SegSure."""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def validate_file_exists(filepath: str, logger=None) -> bool:
    """Check if file exists.
    
    Parameters
    ----------
    filepath : str
        Path to file.
    logger : logging.Logger, optional
        Logger instance for reporting.
    
    Returns
    -------
    bool
        True if file exists, False otherwise.
    """
    exists = os.path.isfile(filepath)
    if not exists and logger:
        logger.warning(f"File not found: {filepath}")
    return exists


def validate_directory_exists(dirpath: str, logger=None) -> bool:
    """Check if directory exists.
    
    Parameters
    ----------
    dirpath : str
        Path to directory.
    logger : logging.Logger, optional
        Logger instance for reporting.
    
    Returns
    -------
    bool
        True if directory exists, False otherwise.
    """
    exists = os.path.isdir(dirpath)
    if not exists and logger:
        logger.warning(f"Directory not found: {dirpath}")
    return exists


def get_file_info(filepath: str) -> Dict:
    """Get basic file information.
    
    Parameters
    ----------
    filepath : str
        Path to file.
    
    Returns
    -------
    dict
        Dictionary with file info (exists, size, path).
    """
    info = {
        "path": filepath,
        "exists": os.path.isfile(filepath),
        "size_bytes": os.path.getsize(filepath) if os.path.isfile(filepath) else None,
        "size_mb": (
            os.path.getsize(filepath) / (1024 * 1024)
            if os.path.isfile(filepath)
            else None
        ),
    }
    return info


def validate_dataframe_columns(
    df: pd.DataFrame,
    required_columns: List[str],
    logger=None
) -> Tuple[bool, List[str]]:
    """Validate that DataFrame has required columns.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate.
    required_columns : list
        List of required column names.
    logger : logging.Logger, optional
        Logger instance for reporting.
    
    Returns
    -------
    bool
        True if all required columns present.
    list
        List of missing columns.
    """
    missing = [col for col in required_columns if col not in df.columns]
    
    if missing and logger:
        logger.warning(f"Missing columns: {missing}")
    
    return len(missing) == 0, missing


def get_dataframe_summary(df: pd.DataFrame) -> Dict:
    """Get summary statistics for a DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to summarize.
    
    Returns
    -------
    dict
        Dictionary with summary info.
    """
    summary = {
        "shape": df.shape,
        "columns": list(df.columns),
        "dtypes": df.dtypes.to_dict(),
        "missing_values": df.isnull().sum().to_dict(),
        "missing_fraction": (df.isnull().sum() / len(df)).to_dict(),
        "memory_usage_mb": df.memory_usage(deep=True).sum() / (1024 * 1024),
    }
    
    # Add numeric column statistics
    numeric_cols = df.select_dtypes(include=np.number).columns
    if len(numeric_cols) > 0:
        summary["numeric_summary"] = df[numeric_cols].describe().to_dict()
    
    return summary


def get_coordinate_ranges(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    z_col: Optional[str] = None
) -> Dict:
    """Get coordinate ranges from DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with coordinate columns.
    x_col : str
        Name of X coordinate column.
    y_col : str
        Name of Y coordinate column.
    z_col : str, optional
        Name of Z coordinate column.
    
    Returns
    -------
    dict
        Dictionary with coordinate ranges.
    """
    ranges = {
        "x": (df[x_col].min(), df[x_col].max()),
        "y": (df[y_col].min(), df[y_col].max()),
    }
    
    if z_col and z_col in df.columns:
        ranges["z"] = (df[z_col].min(), df[z_col].max())
    
    return ranges


def validate_coordinates(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    z_col: Optional[str] = None,
    logger=None
) -> bool:
    """Validate coordinate columns are numeric and non-null.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with coordinate columns.
    x_col : str
        Name of X coordinate column.
    y_col : str
        Name of Y coordinate column.
    z_col : str, optional
        Name of Z coordinate column.
    logger : logging.Logger, optional
        Logger instance for reporting.
    
    Returns
    -------
    bool
        True if coordinates are valid.
    """
    coord_cols = [x_col, y_col]
    if z_col:
        coord_cols.append(z_col)
    
    valid = True
    for col in coord_cols:
        if col not in df.columns:
            if logger:
                logger.error(f"Coordinate column not found: {col}")
            valid = False
        elif not pd.api.types.is_numeric_dtype(df[col]):
            if logger:
                logger.error(f"Coordinate column not numeric: {col}")
            valid = False
        elif df[col].isnull().any():
            null_count = df[col].isnull().sum()
            if logger:
                logger.warning(f"Coordinate column has nulls: {col} ({null_count})")
            valid = False
    
    return valid
