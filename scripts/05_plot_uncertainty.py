#!/usr/bin/env python
"""
Plot segmentation uncertainty and disagreement maps.

This script generates spatial visualizations of disagreement metrics
and diagnostic neighborhoods for high-discord regions.
"""

import argparse
import json
import os

import pandas as pd
import yaml

from segsure.io import atomx, proseg
from segsure.plotting import neighborhood, spatial
from segsure.uncertainty import disagreement
from segsure.utils import logging as logging_util


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def plot_uncertainty_main(
    atomx_root: str,
    atomx_files: dict,
    proseg_root: str,
    proseg_files: dict,
    metrics_dir: str,
    output_dir: str,
    logger=None,
) -> dict:
    """Generate uncertainty visualizations."""
    logger.info("=" * 80)
    logger.info("PLOTTING UNCERTAINTY AND DISAGREEMENT MAPS")
    logger.info("=" * 80)
    
    result = {"plots": []}
    
    # Load metrics
    logger.info("\n--- Loading Metrics ---")
    
    transcript_file = os.path.join(metrics_dir, "transcript_metrics.csv")
    if not os.path.isfile(transcript_file):
        logger.warning(f"Transcript metrics not found: {transcript_file}")
        transcript_metrics = pd.DataFrame()
    else:
        transcript_metrics = pd.read_csv(transcript_file)
        logger.info(f"Loaded {len(transcript_metrics)} transcript metrics")
    
    neighborhood_file = os.path.join(metrics_dir, "neighborhood_metrics.csv")
    if not os.path.isfile(neighborhood_file):
        logger.warning(f"Neighborhood metrics not found: {neighborhood_file}")
        neighborhood_metrics = pd.DataFrame()
    else:
        neighborhood_metrics = pd.read_csv(neighborhood_file)
        logger.info(f"Loaded {len(neighborhood_metrics)} neighborhood metrics")
    
    centroid_file = os.path.join(metrics_dir, "centroid_distances.csv")
    if not os.path.isfile(centroid_file):
        logger.warning(f"Centroid metrics not found: {centroid_file}")
        centroid_metrics = pd.DataFrame()
    else:
        centroid_metrics = pd.read_csv(centroid_file)
        logger.info(f"Loaded {len(centroid_metrics)} centroid metrics")
    
    # Load AtoMx coordinates
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
    
    # Infer coordinate columns
    atomx_coords = {"x": None, "y": None}
    for col in atomx_meta.columns:
        if "x" in col.lower() and atomx_coords["x"] is None:
            atomx_coords["x"] = col
        elif "y" in col.lower() and atomx_coords["y"] is None:
            atomx_coords["y"] = col
    
    if not atomx_coords["x"] or not atomx_coords["y"]:
        logger.warning("Could not infer AtoMx coordinate columns")
        return result
    
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "neighborhoods"), exist_ok=True)
    
    # Plot discord distribution
    logger.info("\n--- Creating Discord Distribution Plot ---")
    if len(transcript_metrics) > 0:
        dist_plot_path = os.path.join(output_dir, "discord_distribution.png")
        spatial.plot_discord_distribution(
            disagreement_metrics=transcript_metrics,
            output_path=dist_plot_path,
            logger=logger,
        )
        result["plots"].append("discord_distribution.png")
    
    # Plot spatial uncertainty heatmap
    logger.info("\n--- Creating Uncertainty Heatmap ---")
    if len(transcript_metrics) > 0 and atomx_coords["x"]:
        heatmap_path = os.path.join(output_dir, "uncertainty_heatmap.png")
        spatial.plot_uncertainty_heatmap(
            disagreement_metrics=transcript_metrics,
            cell_coords=atomx_meta,
            coord_cols=atomx_coords,
            color_metric="discord_severity",
            output_path=heatmap_path,
            logger=logger,
        )
        result["plots"].append("uncertainty_heatmap.png")
    
    # Identify high-discord regions
    logger.info("\n--- Identifying High-Discord Regions ---")
    if len(neighborhood_metrics) > 0:
        high_discord_cells = neighborhood.identify_high_discord_regions(
            neighborhood_discord=neighborhood_metrics,
            cell_coords=atomx_meta,
            coord_cols=atomx_coords,
            threshold=0.5,
            logger=logger,
        )
        
        result["high_discord_regions"] = len(high_discord_cells)
        
        # Create neighborhood panels
        if len(high_discord_cells) > 0:
            logger.info("\n--- Creating Neighborhood Panels ---")
            
            # Get neighbors
            neighbors_atomx = {}
            # Would need to compute neighborhoods from coordinates
            
            neighborhood.create_neighborhood_panels(
                high_discord_cells=high_discord_cells,
                cell_coords=atomx_meta,
                coord_cols=atomx_coords,
                disagreement_metrics=transcript_metrics,
                neighbors=neighbors_atomx,
                output_dir=os.path.join(output_dir, "neighborhoods"),
                logger=logger,
            )
            result["plots"].append("neighborhoods/")
    
    # Save plot summary
    logger.info("\n--- Saving Results ---")
    summary_file = os.path.join(output_dir, "plots_summary.json")
    with open(summary_file, "w") as f:
        json.dump(result, f, indent=2, default=str)
    logger.info(f"Saved plot summary to: {summary_file}")
    
    logger.info("=" * 80)
    logger.info("UNCERTAINTY PLOTTING COMPLETE")
    logger.info("=" * 80)
    
    return result


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Plot segmentation uncertainty and disagreement maps"
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
    output_dir = os.path.join(args.output_root, "spatial_maps", "logs")
    os.makedirs(output_dir, exist_ok=True)
    logger = logging_util.setup_logger(
        "plot_uncertainty",
        output_dir=output_dir,
        verbose=args.verbose,
    )
    
    logger.info("Starting uncertainty visualization")
    logger.info(f"Dataset ID: {args.dataset_id}")
    
    # Load configuration
    config = load_config(args.config)
    dataset_config = config["datasets"][args.dataset_id]
    
    # Generate plots
    plot_uncertainty_main(
        atomx_root=dataset_config["atomx"]["root"],
        atomx_files=dataset_config["atomx"],
        proseg_root=dataset_config["proseg"]["root"],
        proseg_files=dataset_config["proseg"],
        metrics_dir=os.path.join(args.output_root, "tables"),
        output_dir=os.path.join(args.output_root, "spatial_maps"),
        logger=logger,
    )


if __name__ == "__main__":
    main()
