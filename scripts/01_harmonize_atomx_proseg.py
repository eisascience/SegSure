#!/usr/bin/env python
"""
Harmonize coordinate and identifier systems between AtoMx and Proseg.

This script loads AtoMx and Proseg data, normalizes coordinates and IDs,
and produces harmonized datasets ready for matching and comparison.
"""

import argparse
import json
import os
import sys

import pandas as pd
import yaml

from segsure.harmonize import cells, coordinates, transcripts
from segsure.io import atomx, proseg
from segsure.utils import logging as logging_util


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def harmonize_atomx_proseg(
    atomx_root: str,
    atomx_files: dict,
    proseg_root: str,
    proseg_files: dict,
    output_dir: str,
    logger,
) -> dict:
    """Harmonize AtoMx and Proseg data."""
    logger.info("=" * 80)
    logger.info("HARMONIZING ATOMX AND PROSEG DATA")
    logger.info("=" * 80)
    
    result = {
        "atomx": {},
        "proseg": {},
        "harmonization": {},
    }
    
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
    
    # Extract coordinates
    logger.info("\n--- Analyzing AtoMx Coordinates ---")
    atomx_meta = atomx_data["metadata"]
    atomx_coords = coordinates.get_coordinate_columns(atomx_meta, logger=logger)
    atomx_fov_col = coordinates.get_fov_column(atomx_meta, logger=logger)
    
    result["atomx"]["coordinate_columns"] = atomx_coords
    result["atomx"]["fov_column"] = atomx_fov_col
    
    # Load Proseg data
    logger.info("\n--- Loading Proseg Data ---")
    proseg_loader = proseg.ProsegLoader(proseg_root, logger=logger)
    
    proseg_h5ad = proseg_loader.load_h5ad(proseg_files["h5ad"])
    if proseg_h5ad is None:
        proseg_h5ad = proseg_loader.load_zarr(proseg_files["zarr"])
    
    if proseg_h5ad is None:
        logger.error("Could not load Proseg data")
        return result
    
    # Extract Proseg coordinates
    logger.info("\n--- Analyzing Proseg Coordinates ---")
    proseg_coords = {"x": None, "y": None, "z": None}
    
    if "spatial" in proseg_h5ad.obsm:
        spatial = proseg_h5ad.obsm["spatial"]
        proseg_coords["x"] = "spatial_x"
        proseg_coords["y"] = "spatial_y"
        if spatial.shape[1] >= 3:
            proseg_coords["z"] = "spatial_z"
        logger.info(f"Found spatial coordinates in obsm: {proseg_coords}")
    
    # Check obs columns for coordinate info
    for col in proseg_h5ad.obs.columns:
        if "x" in col.lower() and proseg_coords["x"] is None:
            proseg_coords["x"] = col
        elif "y" in col.lower() and proseg_coords["y"] is None:
            proseg_coords["y"] = col
    
    result["proseg"]["coordinate_columns"] = proseg_coords
    
    # Harmonize coordinates
    logger.info("\n--- Harmonizing Coordinates ---")
    
    if atomx_coords["x"] and atomx_coords["y"]:
        atomx_norm = coordinates.normalize_coordinates(
            atomx_meta,
            atomx_coords["x"],
            atomx_coords["y"],
            atomx_coords.get("z"),
            logger=logger,
        )
    
    if proseg_coords["x"] and proseg_coords["y"]:
        # Create DataFrame from Proseg coordinates
        proseg_coords_df = pd.DataFrame({
            proseg_coords["x"]: proseg_h5ad.obsm.get("spatial", pd.DataFrame()).T[0]
            if "spatial" in proseg_h5ad.obsm else [],
            proseg_coords["y"]: proseg_h5ad.obsm.get("spatial", pd.DataFrame()).T[1]
            if "spatial" in proseg_h5ad.obsm else [],
        })
    
    result["harmonization"]["atomx_normalized_coordinates"] = (
        atomx_norm[[atomx_coords["x"], atomx_coords["y"]]].describe().to_dict()
        if atomx_coords["x"] and atomx_coords["y"]
        else None
    )
    
    # Extract and harmonize cell IDs
    logger.info("\n--- Harmonizing Cell IDs ---")
    
    atomx_cell_id_col = cells.get_cell_id_column(atomx_meta, logger=logger)
    if atomx_cell_id_col is None:
        atomx_cell_id_col = atomx_meta.index.name or "cell_id"
    
    proseg_cell_id_col = "cell_id"
    if proseg_h5ad.obs.index.name:
        proseg_cell_id_col = proseg_h5ad.obs.index.name
    
    result["atomx"]["cell_id_column"] = atomx_cell_id_col
    result["proseg"]["cell_id_column"] = proseg_cell_id_col
    
    logger.info(f"AtoMx cell ID column: {atomx_cell_id_col}")
    logger.info(f"Proseg cell ID column: {proseg_cell_id_col}")
    
    # Extract and harmonize transcript assignments
    logger.info("\n--- Harmonizing Transcript Assignments ---")
    
    atomx_tx = atomx_data["transcripts"]
    if atomx_tx is not None:
        atomx_gene_col = transcripts.get_transcript_id_column(atomx_tx, logger=logger)
        atomx_cell_assign_col = transcripts.get_cell_assignment_column(
            atomx_tx, logger=logger
        )
        
        result["atomx"]["gene_column"] = atomx_gene_col
        result["atomx"]["cell_assignment_column"] = atomx_cell_assign_col
    
    # Save harmonized data summaries
    logger.info("\n--- Saving Harmonization Results ---")
    
    os.makedirs(output_dir, exist_ok=True)
    
    harmonization_file = os.path.join(output_dir, "harmonization_info.json")
    with open(harmonization_file, "w") as f:
        json.dump(result, f, indent=2, default=str)
    
    logger.info(f"Saved harmonization info to: {harmonization_file}")
    
    # Save harmonized metadata tables
    if atomx_meta is not None and atomx_coords["x"]:
        atomx_meta_out = os.path.join(output_dir, "atomx_metadata_harmonized.csv")
        atomx_meta.to_csv(atomx_meta_out)
        logger.info(f"Saved harmonized AtoMx metadata to: {atomx_meta_out}")
    
    logger.info("=" * 80)
    logger.info("HARMONIZATION COMPLETE")
    logger.info("=" * 80)
    
    return result


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Harmonize coordinate and identifier systems between AtoMx and Proseg"
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
        help="Dataset ID to harmonize",
    )
    parser.add_argument(
        "--output-root",
        type=str,
        default="results",
        help="Root directory for results",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    
    args = parser.parse_args()
    
    # Setup logging
    output_dir = os.path.join(args.output_root, "harmonized", "logs")
    os.makedirs(output_dir, exist_ok=True)
    logger = logging_util.setup_logger(
        "harmonize_atomx_proseg",
        output_dir=output_dir,
        verbose=args.verbose,
    )
    
    logger.info("Starting SegSure data harmonization")
    logger.info(f"Config: {args.config}")
    logger.info(f"Dataset ID: {args.dataset_id}")
    logger.info(f"Output root: {args.output_root}")
    
    # Load configuration
    config = load_config(args.config)
    dataset_config = config["datasets"][args.dataset_id]
    
    # Harmonize data
    harmonize_atomx_proseg(
        atomx_root=dataset_config["atomx"]["root"],
        atomx_files=dataset_config["atomx"],
        proseg_root=dataset_config["proseg"]["root"],
        proseg_files=dataset_config["proseg"],
        output_dir=os.path.join(args.output_root, "harmonized"),
        logger=logger,
    )


if __name__ == "__main__":
    main()
