#!/usr/bin/env python
"""
Compute disagreement metrics for segmentation comparison.

This script computes geometric, transcript assignment, and neighborhood
disagreement metrics across all cell pairs.
"""

import argparse
import json
import os

import pandas as pd
import yaml

from segsure.io import atomx, proseg
from segsure.metrics import geometry, neighborhood, transcript_assignment
from segsure.utils import logging as logging_util


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def compute_metrics_main(
    atomx_root: str,
    atomx_files: dict,
    proseg_root: str,
    proseg_files: dict,
    cell_matches_file: str,
    transcript_comparison_file: str,
    output_dir: str,
    logger=None,
) -> dict:
    """Compute disagreement metrics."""
    logger.info("=" * 80)
    logger.info("COMPUTING DISAGREEMENT METRICS")
    logger.info("=" * 80)
    
    result = {}
    
    # Load cell matches
    logger.info("\n--- Loading Data ---")
    if not os.path.isfile(cell_matches_file):
        logger.error(f"Cell matches file not found: {cell_matches_file}")
        return result
    
    matches = pd.read_csv(cell_matches_file)
    logger.info(f"Loaded {len(matches)} cell matches")
    
    # Load transcript comparisons if available
    transcript_comparisons = None
    if os.path.isfile(transcript_comparison_file):
        transcript_comparisons = pd.read_csv(transcript_comparison_file)
        logger.info(f"Loaded {len(transcript_comparisons)} transcript comparisons")
    
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
    
    # Load Proseg data
    logger.info("\n--- Loading Proseg Data ---")
    proseg_loader = proseg.ProsegLoader(proseg_root, logger=logger)
    proseg_h5ad = proseg_loader.load_h5ad(proseg_files["h5ad"])
    if proseg_h5ad is None:
        proseg_h5ad = proseg_loader.load_zarr(proseg_files["zarr"])
    
    if proseg_h5ad is None:
        logger.error("Could not load Proseg data")
        return result
    
    proseg_meta = proseg_h5ad.obs.copy()
    
    # Extract coordinates
    atomx_coords = {"x": None, "y": None}
    for col in atomx_meta.columns:
        if "x" in col.lower() and atomx_coords["x"] is None:
            atomx_coords["x"] = col
        elif "y" in col.lower() and atomx_coords["y"] is None:
            atomx_coords["y"] = col
    
    proseg_coords = {"x": None, "y": None}
    if "spatial" in proseg_h5ad.obsm:
        spatial = proseg_h5ad.obsm["spatial"]
        proseg_meta["spatial_x"] = spatial[:, 0]
        proseg_meta["spatial_y"] = spatial[:, 1]
        proseg_coords["x"] = "spatial_x"
        proseg_coords["y"] = "spatial_y"
    
    # Compute centroid distances
    logger.info("\n--- Computing Centroid Distances ---")
    centroid_metrics = geometry.compute_centroid_distance(
        atomx_cells=atomx_meta.set_index(atomx_meta.index),
        proseg_cells=proseg_meta.set_index(proseg_meta.index),
        atomx_coords=atomx_coords,
        proseg_coords=proseg_coords,
        cell_matches=matches,
        logger=logger,
    )
    
    result["centroid_metrics"] = {
        "computed": len(centroid_metrics) > 0,
        "count": len(centroid_metrics),
    }
    
    # Compute transcript overlap metrics
    logger.info("\n--- Computing Transcript Overlap Metrics ---")
    if transcript_comparisons is not None:
        tx_metrics = transcript_assignment.compute_transcript_overlap_metrics(
            assignment_comparisons=transcript_comparisons,
            logger=logger,
        )
        
        # Classify discord severity
        tx_metrics = transcript_assignment.compute_discord_severity(
            assignment_comparisons=tx_metrics,
            logger=logger,
        )
        
        result["transcript_metrics"] = {
            "computed": len(tx_metrics) > 0,
            "count": len(tx_metrics),
        }
    else:
        tx_metrics = matches.copy()
        tx_metrics["jaccard_index"] = 0.0
        result["transcript_metrics"] = {"computed": False, "count": 0}
    
    # Identify neighbors
    logger.info("\n--- Identifying Neighborhoods ---")
    if atomx_coords["x"] and atomx_coords["y"]:
        neighbors_atomx = neighborhood.identify_neighbors(
            cell_coords=atomx_meta,
            coord_cols=atomx_coords,
            neighbor_distance=50.0,
            logger=logger,
        )
        
        result["neighbors_atomx"] = {
            "identified": True,
            "count": len(neighbors_atomx),
        }
    else:
        neighbors_atomx = {}
        result["neighbors_atomx"] = {"identified": False}
    
    if proseg_coords["x"] and proseg_coords["y"]:
        neighbors_proseg = neighborhood.identify_neighbors(
            cell_coords=proseg_meta,
            coord_cols=proseg_coords,
            neighbor_distance=50.0,
            logger=logger,
        )
        
        result["neighbors_proseg"] = {
            "identified": True,
            "count": len(neighbors_proseg),
        }
    else:
        neighbors_proseg = {}
        result["neighbors_proseg"] = {"identified": False}
    
    # Compute neighborhood discord
    logger.info("\n--- Computing Neighborhood Discord ---")
    neighborhood_metrics = neighborhood.compute_neighborhood_discord(
        assignment_comparisons=tx_metrics,
        neighbors_atomx=neighbors_atomx,
        neighbors_proseg=neighbors_proseg,
        cell_matches=matches,
        logger=logger,
    )
    
    result["neighborhood_metrics"] = {
        "computed": len(neighborhood_metrics) > 0,
        "count": len(neighborhood_metrics),
    }
    
    # Save results
    logger.info("\n--- Saving Results ---")
    os.makedirs(output_dir, exist_ok=True)
    
    if len(centroid_metrics) > 0:
        centroid_file = os.path.join(output_dir, "centroid_distances.csv")
        centroid_metrics.to_csv(centroid_file, index=False)
        logger.info(f"Saved centroid distances to: {centroid_file}")
    
    if len(tx_metrics) > 0:
        tx_file = os.path.join(output_dir, "transcript_metrics.csv")
        tx_metrics.to_csv(tx_file, index=False)
        logger.info(f"Saved transcript metrics to: {tx_file}")
    
    if len(neighborhood_metrics) > 0:
        neigh_file = os.path.join(output_dir, "neighborhood_metrics.csv")
        neighborhood_metrics.to_csv(neigh_file, index=False)
        logger.info(f"Saved neighborhood metrics to: {neigh_file}")
    
    summary_file = os.path.join(output_dir, "metrics_summary.json")
    with open(summary_file, "w") as f:
        json.dump(result, f, indent=2, default=str)
    logger.info(f"Saved summary to: {summary_file}")
    
    logger.info("=" * 80)
    logger.info("METRICS COMPUTATION COMPLETE")
    logger.info("=" * 80)
    
    return result


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compute disagreement metrics for segmentation comparison"
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
    output_dir = os.path.join(args.output_root, "tables", "logs")
    os.makedirs(output_dir, exist_ok=True)
    logger = logging_util.setup_logger(
        "compute_metrics",
        output_dir=output_dir,
        verbose=args.verbose,
    )
    
    logger.info("Starting metrics computation")
    logger.info(f"Dataset ID: {args.dataset_id}")
    
    # Load configuration
    config = load_config(args.config)
    dataset_config = config["datasets"][args.dataset_id]
    
    # Compute metrics
    compute_metrics_main(
        atomx_root=dataset_config["atomx"]["root"],
        atomx_files=dataset_config["atomx"],
        proseg_root=dataset_config["proseg"]["root"],
        proseg_files=dataset_config["proseg"],
        cell_matches_file=os.path.join(args.output_root, "tables", "cell_matches.csv"),
        transcript_comparison_file=os.path.join(
            args.output_root, "tables", "transcript_comparisons.csv"
        ),
        output_dir=os.path.join(args.output_root, "tables"),
        logger=logger,
    )


if __name__ == "__main__":
    main()
