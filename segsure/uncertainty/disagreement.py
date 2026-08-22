"""Disagreement aggregation for SegSure."""

from typing import Dict, Optional

import pandas as pd


def aggregate_disagreement_metrics(
    centroid_metrics: Optional[pd.DataFrame] = None,
    transcript_metrics: Optional[pd.DataFrame] = None,
    geometry_metrics: Optional[pd.DataFrame] = None,
    neighborhood_metrics: Optional[pd.DataFrame] = None,
    cell_relationships: Optional[Dict] = None,
    logger=None
) -> pd.DataFrame:
    """Aggregate all disagreement metrics into comprehensive table.
    
    Parameters
    ----------
    centroid_metrics : pd.DataFrame, optional
        Centroid distance metrics.
    transcript_metrics : pd.DataFrame, optional
        Transcript assignment discord metrics.
    geometry_metrics : pd.DataFrame, optional
        Geometric disagreement metrics.
    neighborhood_metrics : pd.DataFrame, optional
        Neighborhood discord metrics.
    cell_relationships : dict, optional
        Cell relationship summary (splits, merges, etc.).
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Aggregated disagreement metrics per cell pair.
    """
    if logger:
        logger.info("Aggregating disagreement metrics")
    
    # Start with transcript metrics as base
    result = transcript_metrics.copy() if transcript_metrics is not None else pd.DataFrame()
    
    # Merge other metrics
    if centroid_metrics is not None and len(result) > 0:
        result = result.merge(
            centroid_metrics[["atomx_cell", "centroid_distance"]],
            on="atomx_cell",
            how="left"
        )
    
    if geometry_metrics is not None and len(result) > 0:
        result = result.merge(
            geometry_metrics[["atomx_cell", "iou", "overlap_area"]],
            on="atomx_cell",
            how="left"
        )
    
    if neighborhood_metrics is not None and len(result) > 0:
        result = result.merge(
            neighborhood_metrics[[
                "atomx_cell", "neighbor_discord_density", "high_discord_neighbors"
            ]],
            on="atomx_cell",
            how="left"
        )
    
    # Add cell relationship info
    if cell_relationships is not None and len(result) > 0:
        def get_relationship_type(atomx_id):
            for match in cell_relationships.get("one_to_one", []):
                if match.get("atomx_cell") == atomx_id:
                    return "one_to_one"
            for match in cell_relationships.get("splits", []):
                if match.get("atomx_cell") == atomx_id:
                    return "split"
            for match in cell_relationships.get("merges", []):
                if match.get("atomx_cell") == atomx_id:
                    return "merge"
            return "unknown"
        
        result["relationship_type"] = result["atomx_cell"].apply(
            get_relationship_type
        )
    
    if logger:
        logger.info(f"Aggregated metrics for {len(result)} cell pairs")
    
    return result


def compute_overall_uncertainty_score(
    disagreement_metrics: pd.DataFrame,
    weights: Optional[Dict[str, float]] = None,
    logger=None
) -> pd.DataFrame:
    """Compute overall uncertainty score (composite metric).
    
    Note: This is optional and kept as a single composite metric
    alongside the decomposed disagreement components.
    
    Parameters
    ----------
    disagreement_metrics : pd.DataFrame
        Aggregated disagreement metrics.
    weights : dict, optional
        Weights for different metric components. If None, uses equal weighting.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Metrics with added uncertainty_score column.
    """
    if logger:
        logger.info("Computing overall uncertainty score")
    
    result = disagreement_metrics.copy()
    
    # Define default weights
    if weights is None:
        weights = {
            "discord_rate": 0.3,
            "centroid_distance": 0.2,
            "neighbor_discord_density": 0.2,
            "iou": -0.3,  # Negative because higher IoU is better
        }
    
    # Normalize and combine metrics
    # This is a simple linear combination; more sophisticated methods can be added
    score = 0.0
    weight_sum = 0.0
    
    for metric, weight in weights.items():
        if metric in result.columns:
            normalized = result[metric].fillna(0).values
            if metric == "iou":
                # IoU is inverted (lower = more discord)
                normalized = 1.0 - normalized
            else:
                # Normalize to [0, 1]
                max_val = normalized.max()
                if max_val > 0:
                    normalized = normalized / max_val
            
            score = score + weight * normalized
            weight_sum += abs(weight)
    
    # Scale to [0, 1]
    if weight_sum > 0:
        result["uncertainty_score"] = score / weight_sum
    else:
        result["uncertainty_score"] = 0.0
    
    if logger:
        logger.info(
            f"Uncertainty scores range: [{result['uncertainty_score'].min():.3f}, "
            f"{result['uncertainty_score'].max():.3f}]"
        )
    
    return result


def classify_cell_uncertainty(
    disagreement_metrics: pd.DataFrame,
    uncertainty_col: str = "discord_severity",
    logger=None
) -> pd.DataFrame:
    """Classify cells by uncertainty category.
    
    Parameters
    ----------
    disagreement_metrics : pd.DataFrame
        Disagreement metrics.
    uncertainty_col : str, optional
        Column to use for classification. Default is discord_severity.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Metrics with uncertainty classification.
    """
    if logger:
        logger.info(f"Classifying cell uncertainty using '{uncertainty_col}'")
    
    result = disagreement_metrics.copy()
    
    if uncertainty_col in result.columns:
        uncertainty_counts = result[uncertainty_col].value_counts()
        if logger:
            logger.info(f"Uncertainty distribution:\n{uncertainty_counts}")
    
    return result
