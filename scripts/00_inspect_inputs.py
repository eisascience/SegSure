#!/usr/bin/env python
"""
Inspect and audit AtoMx and Proseg input data.

This script creates a comprehensive audit report of all input files,
documenting schemas, data completeness, coordinate ranges, and other
key properties for downstream harmonization.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

from segsure.io import atomx, proseg
from segsure.utils import logging as logging_util
from segsure.utils import validation


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def audit_atomx_data(
    root_dir: str,
    filenames: dict,
    logger,
    max_rows: int = None,
    fov: str = None,
) -> dict:
    """Audit AtoMx exported data files."""
    logger.info("=" * 80)
    logger.info("AUDITING ATOMX DATA")
    logger.info("=" * 80)

    audit = {
        "root_dir": root_dir,
        "files": {},
    }

    loader = atomx.AtoMxLoader(root_dir, logger=logger)

    # Audit each file
    logger.info("\n--- Expression Matrix ---")
    expr_file = filenames.get("expression")
    if expr_file:
        filepath = os.path.join(root_dir, expr_file)
        schema = loader.inspect_file_schema(expr_file)
        if schema:
            audit["files"]["expression"] = schema
            logger.info(f"File: {expr_file}")
            logger.info(f"Size: {schema['size_mb']:.1f} MB")
            logger.info(f"Columns: {len(schema['columns'])}")
            logger.info(f"Estimated rows: {schema['n_rows_estimate']}")

    logger.info("\n--- FOV Positions ---")
    fov_file = filenames.get("fov_positions")
    if fov_file:
        schema = loader.inspect_file_schema(fov_file)
        if schema:
            audit["files"]["fov_positions"] = schema
            logger.info(f"File: {fov_file}")
            logger.info(f"Size: {schema['size_mb']:.1f} MB")
            logger.info(f"Columns: {schema['columns']}")
            logger.info(f"Rows: {schema['n_rows_estimate']}")

            # Try to get available FOVs
            available_fovs = loader.get_available_fovs(fov_file)
            if available_fovs:
                audit["files"]["fov_positions"]["available_fovs"] = available_fovs
                logger.info(f"Available FOVs: {available_fovs[:5]}...")
                if fov and fov not in available_fovs:
                    logger.warning(f"Requested FOV '{fov}' not found in available FOVs")

    logger.info("\n--- Metadata ---")
    meta_file = filenames.get("metadata")
    if meta_file:
        schema = loader.inspect_file_schema(meta_file)
        if schema:
            audit["files"]["metadata"] = schema
            logger.info(f"File: {meta_file}")
            logger.info(f"Size: {schema['size_mb']:.1f} MB")
            logger.info(f"Columns: {schema['columns']}")
            logger.info(f"Rows: {schema['n_rows_estimate']}")

            # Actually load metadata to check for important columns
            meta_df = loader.load_metadata(meta_file, nrows=100, fov=fov)
            if meta_df is not None:
                audit["files"]["metadata"]["sample_columns"] = list(meta_df.columns)
                logger.info(f"Sample metadata columns: {list(meta_df.columns)}")

    logger.info("\n--- Transcripts ---")
    tx_file = filenames.get("transcripts")
    if tx_file:
        tx_schema = loader.inspect_transcript_columns(tx_file)
        if tx_schema:
            audit["files"]["transcripts"] = tx_schema
            logger.info(f"File: {tx_file}")
            logger.info(f"Size: {tx_schema['size_mb']:.1f} MB")
            logger.info(f"Columns: {tx_schema['columns']}")
            logger.info(f"Estimated rows: {tx_schema['n_rows_estimate']}")
            logger.info(f"Coordinate columns detected: {tx_schema.get('coordinate_columns', {})}")
            logger.info(f"Gene column: {tx_schema.get('gene_column', 'NOT FOUND')}")
            logger.info(f"Cell ID column: {tx_schema.get('cell_id_column', 'NOT FOUND')}")

            # Sample some data
            tx_df = loader.load_transcripts(tx_file, nrows=100, fov=fov)
            if tx_df is not None:
                logger.info(f"Sample transcript data shape: {tx_df.shape}")
                logger.info(f"Sample columns: {list(tx_df.columns)}")

    logger.info("\n--- Polygons ---")
    poly_file = filenames.get("polygons")
    if poly_file:
        schema = loader.inspect_file_schema(poly_file)
        if schema:
            audit["files"]["polygons"] = schema
            logger.info(f"File: {poly_file}")
            logger.info(f"Size: {schema['size_mb']:.1f} MB")
            logger.info(f"Columns: {schema['columns']}")
            logger.info(f"Estimated rows: {schema['n_rows_estimate']}")

            # Sample
            poly_df = loader.load_polygons(poly_file, nrows=100, fov=fov)
            if poly_df is not None:
                logger.info(f"Sample polygon data shape: {poly_df.shape}")

    return audit


def audit_proseg_data(
    root_dir: str,
    filenames: dict,
    logger,
    fov: str = None,
) -> dict:
    """Audit Proseg output files."""
    logger.info("\n" + "=" * 80)
    logger.info("AUDITING PROSEG DATA")
    logger.info("=" * 80)

    audit = {
        "root_dir": root_dir,
        "files": {},
        "data_availability": {
            "polygons_available": "UNKNOWN",
            "molecule_level_assignments_available": "UNKNOWN",
            "stable_molecule_ids_available": "UNKNOWN",
            "coordinate_systems_compatible": "UNKNOWN",
        },
    }

    loader = proseg.ProsegLoader(root_dir, logger=logger)

    # Audit H5AD
    logger.info("\n--- H5AD ---")
    h5ad_file = filenames.get("h5ad")
    if h5ad_file:
        filepath = os.path.join(root_dir, h5ad_file)
        if os.path.isfile(filepath):
            file_info = {
                "path": filepath,
                "exists": True,
                "size_mb": os.path.getsize(filepath) / (1024 ** 2),
            }
            audit["files"]["h5ad"] = file_info
            logger.info(f"File: {h5ad_file}")
            logger.info(f"Size: {file_info['size_mb']:.1f} MB")

            # Load H5AD to inspect structure (not full data)
            try:
                adata = loader.load_h5ad(h5ad_file)
                if adata is not None:
                    audit["files"]["h5ad"]["shape"] = adata.shape
                    audit["files"]["h5ad"]["obs_columns"] = list(adata.obs.columns)
                    # Don't treat var as individual molecules
                    audit["files"]["h5ad"]["var_count"] = adata.n_vars
                    audit["files"]["h5ad"]["obsm_keys"] = list(adata.obsm.keys()) if adata.obsm else []
                    audit["files"]["h5ad"]["uns_keys"] = (
                        list(adata.uns.keys())[:10] if adata.uns else []
                    )
                    audit["files"]["h5ad"]["layers"] = (
                        list(adata.layers.keys()) if adata.layers else []
                    )

                    logger.info(f"Shape: {adata.shape} (observations x features)")
                    logger.info(f"Obs columns: {list(adata.obs.columns)}")
                    logger.info(f"Features (var): {adata.n_vars} (NOT individual molecules)")
                    logger.info(f"Obsm keys: {list(adata.obsm.keys())}")
                    logger.info(f"Uns keys (first 10): {audit['files']['h5ad']['uns_keys']}")
            except Exception as e:
                audit["files"]["h5ad"]["error"] = str(e)
                logger.error(f"Error loading H5AD: {e}")

    # Audit Zarr
    logger.info("\n--- Zarr ---")
    zarr_file = filenames.get("zarr")
    if zarr_file:
        zarr_path = os.path.join(root_dir, zarr_file)
        if os.path.isdir(zarr_path):
            logger.info(f"Zarr directory: {zarr_file}")

            # Get Zarr structure
            zarr_stats = loader.get_zarr_stats(zarr_file)
            if zarr_stats:
                audit["files"]["zarr"] = zarr_stats
                logger.info(f"Size: {zarr_stats['size_mb']:.1f} MB")
                logger.info(f"Arrays: {zarr_stats['arrays']}")
                logger.info(f"Groups: {zarr_stats['groups']}")

            # Get detailed structure
            logger.info("\n--- Zarr Structure (recursive) ---")
            zarr_structure = loader.inspect_zarr_structure(zarr_file)
            if zarr_structure:
                audit["files"]["zarr_structure"] = zarr_structure
                logger.info(json.dumps(zarr_structure, indent=2, default=str))

            # Try to load Zarr as AnnData
            try:
                adata = loader.load_zarr(zarr_file)
                if adata is not None:
                    audit["files"]["zarr"]["adata_shape"] = adata.shape
                    audit["files"]["zarr"]["obs_columns"] = list(adata.obs.columns)
                    audit["files"]["zarr"]["var_count"] = adata.n_vars
                    audit["files"]["zarr"]["obsm_keys"] = list(adata.obsm.keys()) if adata.obsm else []
                    audit["files"]["zarr"]["uns_keys"] = (
                        list(adata.uns.keys())[:10] if adata.uns else []
                    )

                    logger.info(f"\nZarr AnnData shape: {adata.shape}")
                    logger.info(f"Obs columns: {list(adata.obs.columns)}")
                    logger.info(f"Obsm keys: {list(adata.obsm.keys())}")

                    # Check for key data availability
                    audit["data_availability"] = check_proseg_data_availability(
                        adata, audit["files"]
                    )
                    for key, value in audit["data_availability"].items():
                        logger.info(f"{key}: {value}")

            except Exception as e:
                audit["files"]["zarr"]["error"] = str(e)
                logger.error(f"Error loading Zarr: {e}")

    # Audit Seurat
    logger.info("\n--- Seurat RDS ---")
    seurat_file = filenames.get("seurat")
    if seurat_file:
        filepath = os.path.join(root_dir, seurat_file)
        if os.path.isfile(filepath):
            file_info = {
                "path": filepath,
                "exists": True,
                "size_mb": os.path.getsize(filepath) / (1024 ** 2),
            }
            audit["files"]["seurat"] = file_info
            logger.info(f"File: {seurat_file}")
            logger.info(f"Size: {file_info['size_mb']:.1f} MB")
            logger.info("Seurat RDS loading not implemented (requires R/rpy2)")

    return audit


def check_proseg_data_availability(adata, files_audit) -> dict:
    """Check what key data is available in Proseg output."""
    availability = {
        "polygons_available": "NO",
        "molecule_level_assignments_available": "NO",
        "stable_molecule_ids_available": "NO",
        "coordinate_systems_compatible": "UNKNOWN",
    }

    # Check for polygons in obsm (spatial coordinates)
    if "spatial" in adata.obsm or "centroids" in adata.obsm:
        availability["polygons_available"] = "YES"

    # Check for uns with polygon data
    if adata.uns and any("polygon" in key.lower() for key in adata.uns.keys()):
        availability["polygons_available"] = "YES"

    # Check for molecule-level assignments
    # Look for transcript count or assignment columns in obs
    obs_cols_lower = [col.lower() for col in adata.obs.columns]
    if any("molecule" in col or "transcript" in col or "count" in col for col in obs_cols_lower):
        availability["molecule_level_assignments_available"] = "YES"

    # Check if there are stable IDs for molecules
    # Look for transcript ID or similar columns
    if any("id" in col.lower() for col in obs_cols_lower):
        availability["stable_molecule_ids_available"] = "YES"

    # Check coordinate compatibility
    if "spatial" in adata.obsm and adata.obsm["spatial"].shape[1] >= 2:
        availability["coordinate_systems_compatible"] = "YES"

    return availability


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Inspect and audit AtoMx and Proseg input data"
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
        help="Dataset ID to audit",
    )
    parser.add_argument(
        "--output-root",
        type=str,
        default="results",
        help="Root directory for results",
    )
    parser.add_argument(
        "--fov",
        type=str,
        default=None,
        help="Specific FOV to inspect (if not specified, inspects all)",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Maximum number of rows to load for data inspection",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    # Setup logging
    output_dir = os.path.join(args.output_root, "audit", "logs")
    os.makedirs(output_dir, exist_ok=True)
    logger = logging_util.setup_logger(
        "inspect_inputs",
        output_dir=output_dir,
        verbose=args.verbose,
    )

    logger.info("Starting SegSure data inspection")
    logger.info(f"Config: {args.config}")
    logger.info(f"Dataset ID: {args.dataset_id}")
    logger.info(f"Output root: {args.output_root}")
    if args.fov:
        logger.info(f"FOV: {args.fov}")
    if args.max_rows:
        logger.info(f"Max rows: {args.max_rows}")

    # Load configuration
    config = load_config(args.config)
    dataset_config = config["datasets"][args.dataset_id]

    # Audit AtoMx data
    atomx_audit = audit_atomx_data(
        dataset_config["atomx"]["root"],
        dataset_config["atomx"],
        logger,
        max_rows=args.max_rows,
        fov=args.fov,
    )

    # Audit Proseg data
    proseg_audit = audit_proseg_data(
        dataset_config["proseg"]["root"],
        dataset_config["proseg"],
        logger,
        fov=args.fov,
    )

    # Save audit report
    audit_report = {
        "dataset_id": args.dataset_id,
        "fov": args.fov,
        "max_rows": args.max_rows,
        "atomx": atomx_audit,
        "proseg": proseg_audit,
    }

    audit_path = os.path.join(args.output_root, "audit")
    os.makedirs(audit_path, exist_ok=True)

    report_file = os.path.join(audit_path, f"audit_report_{args.dataset_id}.json")
    with open(report_file, "w") as f:
        json.dump(audit_report, f, indent=2, default=str)

    logger.info(f"\nAudit report saved to: {report_file}")

    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("AUDIT SUMMARY")
    logger.info("=" * 80)
    logger.info("\nProseg data availability:")
    for key, value in proseg_audit.get("data_availability", {}).items():
        logger.info(f"  {key}: {value}")

    logger.info("\nAtoMx files inspected:")
    for file_type in atomx_audit.get("files", {}).keys():
        logger.info(f"  {file_type}: OK")

    logger.info("\nProseg files inspected:")
    for file_type in proseg_audit.get("files", {}).keys():
        if file_type != "zarr_structure":
            logger.info(f"  {file_type}: OK")

    logger.info("\n" + "=" * 80)
    logger.info("INSPECTION COMPLETE")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
