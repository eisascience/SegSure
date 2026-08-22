#!/usr/bin/env python
"""
Match cells between AtoMx and Proseg segmentations.

This script identifies corresponding cells across segmentation methods
and classifies relationship types (one-to-one, splits, merges, lost, new, complex).
"""

import argparse
import json
import os

import pandas as pd
import yaml

from segsure.harmonize import cells, coordinates
from segsure.io import atomx, proseg
from segsure.utils import logging as logging_util


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def filter_by_fov(
    df: pd.DataFrame,
    fov: str,
    fov_col: str = "fov",
    logger=None
) -> pd.DataFrame:
    """Filter data by FOV.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data to filter.
    fov : str
        FOV identifier to keep.
    fov_col : str
        FOV column name.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Filtered data.
    """
    if fov_col in df.columns:
        df_fov = df[df[fov_col] == fov].copy()
        if logger:
            logger.info(f"Filtered to FOV {fov}: {len(df_fov)} rows")
        return df_fov
    else:
        if logger:
            logger.warning(f"FOV column '{fov_col}' not found, using all data")
        return df


def match_cells_main(
    atomx_root: str,
    atomx_files: dict,
    proseg_root: str,
    proseg_files: dict,
    output_dir: str,
    fov: str = None,
    distance_threshold: float = 10.0,
    logger=None,
) -> dict:
    """Match cells between AtoMx and Proseg."""
    logger.info("=" * 80)
    logger.info("MATCHING CELLS BETWEEN ATOMX AND PROSEG")
    logger.info("=" * 80)
    
    result = {"matches": None, "relationships": None}
    
    # Load AtoMx data
    logger.info("\n--- Loading AtoMx Data ---")
    atomx_loader = atomx.AtoMxLoader(atomx_root, logger=logger)
    atomx_data = atomx_loader.load_all(
        expression=atomx_files["expression"],
        fov_positions=atomx_files["fov_positions"],
        metadata=atomx_files["metadata"],
        transcripts=atomx_files["transcripts"],
        polygons=atomx_files["polygons"],
    )
    
    atomx_meta = atomx_data["metadata"]
    if atomx_meta is None:
        logger.error("Could not load AtoMx metadata")
        return result
    
    # Get coordinate and cell ID columns for AtoMx
    atomx_coords = coordinates.get_coordinate_columns(atomx_meta, logger=logger)
    atomx_cell_id_col = cells.get_cell_id_column(atomx_meta, logger=logger)
    if atomx_cell_id_col is None:
        atomx_cell_id_col = atomx_meta.index.name or "atomx_cell_id"
        atomx_meta[atomx_cell_id_col] = atomx_meta.index
    
    # Load Proseg data
    logger.info("\n--- Loading Proseg Data ---")
    proseg_loader = proseg.ProsegLoader(proseg_root, logger=logger)
    
    proseg_h5ad = proseg_loader.load_h5ad(proseg_files["h5ad"])
    if proseg_h5ad is None:
        proseg_h5ad = proseg_loader.load_zarr(proseg_files["zarr"])
    
    if proseg_h5ad is None:
        logger.error("Could not load Proseg data")
        return result
    
    # Extract Proseg cell data
    proseg_meta = proseg_h5ad.obs.copy()
    proseg_meta["proseg_cell_id"] = proseg_h5ad.obs.index
    
    # Get spatial coordinates from Proseg
    if "spatial" in proseg_h5ad.obsm:
        spatial = proseg_h5ad.obsm["spatial"]
        proseg_meta["spatial_x"] = spatial[:, 0]
        proseg_meta["spatial_y"] = spatial[:, 1]
        proseg_coords = {"x": "spatial_x", "y": "spatial_y"}
    else:
        proseg_coords = {"x": None, "y": None}
        for col in proseg_meta.columns:
            if "x" in col.lower() and proseg_coords["x"] is None:
                proseg_coords["x"] = col
            elif "y" in col.lower() and proseg_coords["y"] is None:
                proseg_coords["y"] = col
    
    # Filter by FOV if specified
    if fov:
        logger.info(f"\n--- Filtering to FOV {fov} ---")
        atomx_meta = filter_by_fov(atomx_meta, fov, logger=logger)
        proseg_meta = filter_by_fov(proseg_meta, fov, logger=logger)
        if len(atomx_meta) == 0 or len(proseg_meta) == 0:
            logger.error(f"No data found for FOV {fov}")
            return result
    
    # Parse polygons if available
    logger.info("\n--- Parsing Polygons ---")
    atomx_polygons = {}
    proseg_polygons = {}
    
    if atomx_data["polygons"] is not None:
        atomx_polygons = cells.parse_polygons(
            atomx_data["polygons"],
            cell_id_col=atomx_cell_id_col,
            geometry_col="geometry",
            logger=logger
        )
    
    # For Proseg, check if polygon data available in uns or obsm
    if hasattr(proseg_h5ad, "uns") and "polygons" in proseg_h5ad.uns:
        proseg_poly_data = proseg_h5ad.uns["polygons"]
        if isinstance(proseg_poly_data, pd.DataFrame):
            proseg_polygons = cells.parse_polygons(
                proseg_poly_data,
                cell_id_col="proseg_cell_id",
                geometry_col="geometry",
                logger=logger
            )
    
    # Match cells
    logger.info("\n--- Matching Cells ---")
    
    if (
        atomx_coords["x"]
        and atomx_coords["y"]
        and proseg_coords["x"]
        and proseg_coords["y"]
    ):
        matches_df = cells.match_cells(
            atomx_cells=atomx_meta,
            proseg_cells=proseg_meta,
            atomx_coords=atomx_coords,
            proseg_coords=proseg_coords,
            atomx_cell_id=atomx_cell_id_col,
            proseg_cell_id="proseg_cell_id",
            atomx_polygons=atomx_polygons if atomx_polygons else None,
            proseg_polygons=proseg_polygons if proseg_polygons else None,
            distance_threshold=distance_threshold,
            logger=logger,
        )
        
        result["matches"] = matches_df
        
        # Identify split/merge relationships (legacy)
        logger.info("\n--- Identifying Split/Merge Events ---")
        relationships = cells.identify_split_merges(
            matches_df,
            atomx_cell_id=atomx_cell_id_col,
            proseg_cell_id="proseg_cell_id",
            logger=logger,
        )
        
        result["relationships"] = relationships
        
        # Save results
        logger.info("\n--- Saving Results ---")
        os.makedirs(output_dir, exist_ok=True)
        
        # Add FOV to results if available
        if atomx_coords.get("fov"):
            if len(atomx_meta) > 0:
                atomx_fov = atomx_meta[atomx_coords["fov"]].iloc[0]
                matches_df["fov"] = atomx_fov
        
        # CSV format
        matches_file_csv = os.path.join(output_dir, "cell_matches.csv.gz")
        matches_df.to_csv(matches_file_csv, index=False, compression="gzip")
        logger.info(f"Saved cell matches (CSV) to: {matches_file_csv}")
        
        # Parquet format
        matches_file_parquet = os.path.join(output_dir, "cell_matches.parquet")
        matches_df.to_parquet(matches_file_parquet, index=False)
        logger.info(f"Saved cell matches (Parquet) to: {matches_file_parquet}")
        
        # Relationships JSON
        relationships_file = os.path.join(output_dir, "cell_relationships.json")
        with open(relationships_file, "w") as f:
            json.dump(relationships, f, indent=2, default=str)
        logger.info(f"Saved cell relationships to: {relationships_file}")
        
        # Summary stats
        logger.info("\n--- Summary Statistics ---")
        logger.info(f"Total cells in correspondence: {len(matches_df)}")
        if "relationship_type" in matches_df.columns:
            rel_counts = matches_df["relationship_type"].value_counts()
            logger.info("Relationship type counts:")
            for rel_type, count in rel_counts.items():
                logger.info(f"  {rel_type}: {count}")
        
        if "iou" in matches_df.columns:
            one_to_one = matches_df[matches_df.get("relationship_type") == "one_to_one"]
            if len(one_to_one) > 0:
                median_iou = one_to_one["iou"].median()
                logger.info(f"Median IoU (one-to-one): {median_iou:.4f}")
    
    logger.info("=" * 80)
    logger.info("CELL MATCHING COMPLETE")
    logger.info("=" * 80)
    
    return result


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Match cells between AtoMx and Proseg segmentations"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config/datasets.yaml",
        help="Path to configuration file",
    )
    parser.add_argument(
        "--dataset-id",
        type=str,
        default="33710_32578_37826_37374_36135",
        help="Dataset ID",
    )
    parser.add_argument(
        "--fov",
        type=str,
        default=None,
        help="FOV identifier (optional)",
    )
    parser.add_argument(
        "--output-root",
        type=str,
        default="results",
        help="Root directory for results",
    )
    parser.add_argument(
        "--distance-threshold",
        type=float,
        default=10.0,
        help="Maximum distance for cell matching",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    
    args = parser.parse_args()
    
    # Setup logging
    output_dir = os.path.join(args.output_root, "tables", "logs")
    os.makedirs(output_dir, exist_ok=True)
    logger = logging_util.setup_logger(
        "match_cells",
        output_dir=output_dir,
        verbose=args.verbose,
    )
    
    logger.info("Starting cell matching")
    logger.info(f"Config: {args.config}")
    logger.info(f"Dataset ID: {args.dataset_id}")
    logger.info(f"FOV: {args.fov}")
    logger.info(f"Distance threshold: {args.distance_threshold}")
    
    # Load configuration
    config = load_config(args.config)
    dataset_config = config["datasets"][args.dataset_id]
    
    # Match cells
    match_cells_main(
        atomx_root=dataset_config["atomx"]["root"],
        atomx_files=dataset_config["atomx"],
        proseg_root=dataset_config["proseg"]["root"],
        proseg_files=dataset_config["proseg"],
        output_dir=os.path.join(args.output_root, "tables"),
        fov=args.fov,
        distance_threshold=args.distance_threshold,
        logger=logger,
    )


if __name__ == "__main__":
    main()
