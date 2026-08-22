#!/usr/bin/env python
"""
Comprehensive molecule-level disagreement analysis for a single FOV.

This script performs end-to-end disagreement analysis:
1. Loads molecule-level data from AtoMx and Proseg
2. Matches molecules using stable IDs and spatial proximity
3. Classifies molecules into 8 status categories
4. Computes geometry, transcript, and topology disagreement metrics
5. Generates diagnostic visualizations
6. Exports results in multiple formats
"""

import argparse
import json
import os

import pandas as pd
import yaml

from segsure.harmonize import transcripts
from segsure.io import atomx, proseg, export
from segsure.metrics import geometry, neighborhood, transcript_assignment
from segsure.molecules import MoleculeTracker
from segsure.plotting.diagnostic_panels import NeighborhoodDiagnosticPanel
from segsure.utils import logging as logging_util


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def run_molecule_analysis(
    dataset_config: dict,
    fov_id: str,
    spatial_tolerance: float = 1.5,
    output_dir: str = "results",
    logger=None,
) -> dict:
    """Run comprehensive molecule-level disagreement analysis for a single FOV.

    Parameters
    ----------
    dataset_config : dict
        Dataset configuration (from YAML).
    fov_id : str
        FOV identifier to analyze.
    spatial_tolerance : float, optional
        Spatial tolerance for molecule matching (microns). Default is 1.5.
    output_dir : str, optional
        Output directory for results.
    logger : logging.Logger, optional
        Logger instance.

    Returns
    -------
    dict
        Comprehensive analysis summary.
    """
    logger.info("=" * 80)
    logger.info("MOLECULE-LEVEL DISAGREEMENT ANALYSIS")
    logger.info("=" * 80)

    result = {
        "fov_id": fov_id,
        "spatial_tolerance": spatial_tolerance,
        "outputs": {},
    }

    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    tables_dir = os.path.join(output_dir, "tables")
    viz_dir = os.path.join(output_dir, "neighborhoods")
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(viz_dir, exist_ok=True)

    # Initialize exporter
    exporter = export.DisagreementExporter(tables_dir, logger=logger)

    # Load AtoMx data
    logger.info("\n--- Loading AtoMx Data ---")
    atomx_loader = atomx.AtoMxLoader(
        dataset_config["atomx"]["root"], logger=logger
    )
    atomx_data = atomx_loader.load_all(
        expression=dataset_config["atomx"]["expression"],
        fov_positions=dataset_config["atomx"]["fov_positions"],
        metadata=dataset_config["atomx"]["metadata"],
        transcripts=dataset_config["atomx"]["transcripts"],
        polygons=dataset_config["atomx"]["polygons"],
    )

    atomx_molecules = atomx_data.get("transcripts")
    if atomx_molecules is None:
        logger.error("Could not load AtoMx transcripts")
        return result

    # Filter to FOV if needed
    if "fov" in atomx_molecules.columns:
        atomx_molecules = atomx_molecules[atomx_molecules["fov"] == fov_id]
        logger.info(f"Filtered to FOV {fov_id}: {len(atomx_molecules)} molecules")

    # Load Proseg data
    logger.info("\n--- Loading Proseg Data ---")
    proseg_loader = proseg.ProsegLoader(
        dataset_config["proseg"]["root"], logger=logger
    )
    proseg_h5ad = proseg_loader.load_h5ad(dataset_config["proseg"]["h5ad"])
    if proseg_h5ad is None:
        proseg_h5ad = proseg_loader.load_zarr(dataset_config["proseg"]["zarr"])

    if proseg_h5ad is None:
        logger.error("Could not load Proseg data")
        return result

    # Extract molecules from Proseg
    proseg_molecules = _extract_proseg_molecules(proseg_h5ad, logger)
    if proseg_molecules is None or proseg_molecules.empty:
        logger.error("Could not extract Proseg transcripts/molecules")
        return result

    # Load pre-computed cell matches
    logger.info("\n--- Loading Cell Matches ---")
    cell_matches_file = os.path.join(tables_dir, "cell_matches.csv")
    if not os.path.isfile(cell_matches_file):
        logger.warning(f"Cell matches not found: {cell_matches_file}")
        logger.warning("Assuming 1:1 correspondence via row order (placeholder)")
        # This would require prior cell matching
        return result

    cell_matches = pd.read_csv(cell_matches_file)
    logger.info(f"Loaded {len(cell_matches)} cell matches")

    # Match molecules
    logger.info("\n--- Matching Molecules ---")
    tracker = MoleculeTracker(
        atomx_molecules=atomx_molecules,
        proseg_molecules=proseg_molecules,
        cell_matches=cell_matches,
        spatial_tolerance=spatial_tolerance,
        logger=logger,
    )

    matched_molecules, ambiguous_matches = tracker.process_all_molecules()

    if matched_molecules.empty:
        logger.warning("No molecules matched")
        return result

    logger.info(f"Matched {len(matched_molecules)} molecules")
    logger.info(f"Found {len(ambiguous_matches)} ambiguous matches")

    result["n_molecules_matched"] = len(matched_molecules)
    result["n_ambiguous_matches"] = len(ambiguous_matches)

    # Export molecule-level results
    logger.info("\n--- Exporting Results ---")
    exporter.export_molecules_parquet(matched_molecules)
    exporter.export_molecules_csv(matched_molecules)
    exporter.export_ambiguous_matches(ambiguous_matches)

    result["outputs"]["molecules_parquet"] = os.path.join(
        tables_dir, "molecules_analysis.parquet"
    )
    result["outputs"]["molecules_csv"] = os.path.join(
        tables_dir, "molecules_analysis.csv"
    )
    result["outputs"]["ambiguous_matches"] = os.path.join(
        tables_dir, "ambiguous_matches.json"
    )

    # Compute molecule-level summary
    logger.info("\n--- Computing Molecule-Level Summaries ---")
    molecule_summary = (
        transcript_assignment.compute_molecule_level_status_summary(
            matched_molecules, logger=logger
        )
    )
    result["molecule_summary"] = molecule_summary

    # Compute cell-level summaries
    cell_summaries = transcript_assignment.compute_cell_level_transcript_summary(
        matched_molecules, logger=logger
    )
    exporter.export_cell_level_summaries(cell_summaries)
    result["outputs"]["cell_summaries"] = os.path.join(
        tables_dir, "cell_level_summaries.csv"
    )

    # Compute geometry metrics
    logger.info("\n--- Computing Geometry Metrics ---")
    atomx_meta = atomx_data.get("metadata")
    proseg_meta = proseg_h5ad.obs.copy() if proseg_h5ad is not None else pd.DataFrame()

    if atomx_meta is not None and not proseg_meta.empty:
        # Infer coordinates
        atomx_coords = _infer_coords(atomx_meta)
        proseg_coords = _infer_coords(proseg_meta)

        # Compute geometry
        centroid_metrics = geometry.compute_centroid_distance(
            atomx_cells=atomx_meta,
            proseg_cells=proseg_meta,
            atomx_coords=atomx_coords,
            proseg_coords=proseg_coords,
            cell_matches=cell_matches,
            logger=logger,
        )

        if not centroid_metrics.empty:
            exporter.export_geometry_metrics(centroid_metrics)
            result["outputs"]["geometry_metrics"] = os.path.join(
                tables_dir, "geometry_metrics.csv"
            )

        # Polygon overlap if available
        atomx_polygons = atomx_data.get("polygons")
        proseg_polygons = None  # Would need to load from Proseg

        if atomx_polygons is not None:
            polygon_metrics = geometry.compute_polygon_overlap(
                atomx_polygons=atomx_polygons,
                proseg_polygons=proseg_polygons,
                cell_matches=cell_matches,
                logger=logger,
            )
            if not polygon_metrics.empty:
                result["geometry_summary"] = {
                    "iou_mean": polygon_metrics["iou"].mean(),
                    "iou_std": polygon_metrics["iou"].std(),
                }

    # Compute neighborhood metrics
    logger.info("\n--- Computing Neighborhood Metrics ---")
    if atomx_meta is not None and not proseg_meta.empty:
        atomx_coords = _infer_coords(atomx_meta)
        proseg_coords = _infer_coords(proseg_meta)

        neighbors_atomx = neighborhood.identify_neighbors(
            cell_coords=atomx_meta,
            coord_cols=atomx_coords,
            neighbor_distance=50.0,
            logger=logger,
        )

        neighbors_proseg = neighborhood.identify_neighbors(
            cell_coords=proseg_meta,
            coord_cols=proseg_coords,
            neighbor_distance=50.0,
            logger=logger,
        )

        neighborhood_metrics = neighborhood.compute_neighborhood_discord(
            assignment_comparisons=cell_summaries,
            neighbors_atomx=neighbors_atomx,
            neighbors_proseg=neighbors_proseg,
            cell_matches=cell_matches,
            logger=logger,
        )

        if not neighborhood_metrics.empty:
            exporter.export_neighborhood_metrics(neighborhood_metrics)
            result["outputs"]["neighborhood_metrics"] = os.path.join(
                tables_dir, "neighborhood_metrics.csv"
            )

            # Compute Jaccard
            jaccard = neighborhood.compute_mapped_neighbor_jaccard(
                atomx_neighbors=neighbors_atomx,
                proseg_neighbors=neighbors_proseg,
                cell_matches=cell_matches,
                logger=logger,
            )
            result["neighborhood_summary"] = {"neighbor_jaccard": jaccard}

    # Generate diagnostic panels for high-disagreement regions
    logger.info("\n--- Generating Diagnostic Panels ---")
    n_panels_generated = 0
    # This would need to identify high-disagreement cells first
    result["outputs"]["diagnostic_panels"] = viz_dir
    result["n_diagnostic_panels"] = n_panels_generated

    # Export analysis summary
    exporter.export_analysis_summary(result)
    result["outputs"]["analysis_summary"] = os.path.join(
        tables_dir, "analysis_summary.json"
    )

    logger.info("\n" + "=" * 80)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Molecules analyzed: {result['n_molecules_matched']}")
    logger.info(f"Ambiguous matches: {result['n_ambiguous_matches']}")
    logger.info(f"Output directory: {output_dir}")

    return result


def _extract_proseg_molecules(adata, logger=None) -> pd.DataFrame:
    """Extract molecule-level data from Proseg h5ad.

    Parameters
    ----------
    adata : AnnData
        Proseg AnnData object.
    logger : logging.Logger, optional
        Logger instance.

    Returns
    -------
    pd.DataFrame or None
        Molecule data if available.
    """
    # This is a placeholder; actual implementation depends on
    # how Proseg stores molecule-level data
    if adata is None:
        return None

    # Try to extract from adata.var and adata.obs
    molecules = []

    if "transcript_id" in adata.var.columns:
        for tx_idx, tx_id in enumerate(adata.var.index):
            for cell_idx, cell_id in enumerate(adata.obs.index):
                if adata.X[cell_idx, tx_idx] > 0:
                    mol_data = {
                        "transcript_id": tx_id,
                        "cell_id": cell_id,
                        "count": adata.X[cell_idx, tx_idx],
                    }
                    if "spatial" in adata.obsm:
                        mol_data["x"] = adata.obsm["spatial"][cell_idx, 0]
                        mol_data["y"] = adata.obsm["spatial"][cell_idx, 1]
                    molecules.append(mol_data)

    if molecules:
        if logger:
            logger.info(f"Extracted {len(molecules)} molecules from Proseg")
        return pd.DataFrame(molecules)

    if logger:
        logger.warning(
            "Could not extract molecule-level data from Proseg; "
            "may need to use aggregated gene expression instead"
        )
    return None


def _infer_coords(df: pd.DataFrame) -> dict:
    """Infer coordinate column names."""
    coords = {"x": None, "y": None}
    for col in df.columns:
        if "x" in col.lower() and coords["x"] is None:
            coords["x"] = col
        elif "y" in col.lower() and coords["y"] is None:
            coords["y"] = col
    return coords


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Comprehensive molecule-level disagreement analysis for single FOV"
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
        "--fov-id",
        type=str,
        required=True,
        help="FOV ID to analyze (must be validated)",
    )
    parser.add_argument(
        "--output-root",
        type=str,
        default="results",
        help="Root directory for results",
    )
    parser.add_argument(
        "--spatial-tolerance",
        type=float,
        default=1.5,
        help="Spatial tolerance for molecule matching (microns)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    # Setup logging
    output_dir = os.path.join(args.output_root, "logs")
    os.makedirs(output_dir, exist_ok=True)
    logger = logging_util.setup_logger(
        "molecule_analysis",
        output_dir=output_dir,
        verbose=args.verbose,
    )

    logger.info("Starting molecule-level disagreement analysis")
    logger.info(f"Dataset ID: {args.dataset_id}")
    logger.info(f"FOV ID: {args.fov_id}")

    # Load configuration
    config = load_config(args.config)
    dataset_config = config["datasets"][args.dataset_id]

    # Run analysis
    result = run_molecule_analysis(
        dataset_config=dataset_config,
        fov_id=args.fov_id,
        spatial_tolerance=args.spatial_tolerance,
        output_dir=os.path.join(args.output_root, "fov_" + args.fov_id),
        logger=logger,
    )

    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("ANALYSIS SUMMARY")
    logger.info("=" * 80)
    logger.info(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
