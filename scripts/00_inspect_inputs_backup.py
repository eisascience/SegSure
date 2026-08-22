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


def audit_atomx_data(root_dir: str, filenames: dict, logger) -> dict:
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
        file_info = validation.get_file_info(filepath)
        audit["files"]["expression"] = file_info
        
        expr_df = loader.load_expression_matrix(expr_file)
        if expr_df is not None:
            audit["files"]["expression"]["shape"] = expr_df.shape
            audit["files"]["expression"]["columns"] = list(expr_df.columns)[:10]  # First 10
            audit["files"]["expression"]["dtypes"] = str(expr_df.dtypes.to_dict())
            audit["files"]["expression"]["memory_mb"] = expr_df.memory_usage(deep=True).sum() / (1024**2)
            logger.info(f"Shape: {expr_df.shape}")
            logger.info(f"First 10 columns: {list(expr_df.columns[:10])}")
            logger.info(f"Data types: {expr_df.dtypes.unique()}")
            logger.info(f"Memory: {file_info['size_mb']:.1f} MB")
    
    logger.info("\n--- FOV Positions ---")
    fov_file = filenames.get("fov_positions")
    if fov_file:
        filepath = os.path.join(root_dir, fov_file)
        file_info = validation.get_file_info(filepath)
        audit["files"]["fov_positions"] = file_info
        
        fov_df = loader.load_fov_positions(fov_file)
        if fov_df is not None:
            audit["files"]["fov_positions"]["shape"] = fov_df.shape
            audit["files"]["fov_positions"]["columns"] = list(fov_df.columns)
            summary = validation.get_dataframe_summary(fov_df)
            audit["files"]["fov_positions"]["n_rows"] = fov_df.shape[0]
            audit["files"]["fov_positions"]["n_cols"] = fov_df.shape[1]
            logger.info(f"Shape: {fov_df.shape}")
            logger.info(f"Columns: {list(fov_df.columns)}")
            logger.info(f"First rows:\n{fov_df.head()}")
    
    logger.info("\n--- Metadata ---")
    meta_file = filenames.get("metadata")
    if meta_file:
        filepath = os.path.join(root_dir, meta_file)
        file_info = validation.get_file_info(filepath)
        audit["files"]["metadata"] = file_info
        
        meta_df = loader.load_metadata(meta_file)
        if meta_df is not None:
            audit["files"]["metadata"]["shape"] = meta_df.shape
            audit["files"]["metadata"]["columns"] = list(meta_df.columns)
            audit["files"]["metadata"]["n_rows"] = meta_df.shape[0]
            audit["files"]["metadata"]["n_cols"] = meta_df.shape[1]
            logger.info(f"Shape: {meta_df.shape}")
            logger.info(f"Columns: {list(meta_df.columns)}")
            logger.info(f"First rows:\n{meta_df.head()}")
    
    logger.info("\n--- Transcripts ---")
    tx_file = filenames.get("transcripts")
    if tx_file:
        filepath = os.path.join(root_dir, tx_file)
        file_info = validation.get_file_info(filepath)
        audit["files"]["transcripts"] = file_info
        
        tx_df = loader.load_transcripts(tx_file)
        if tx_df is not None:
            audit["files"]["transcripts"]["shape"] = tx_df.shape
            audit["files"]["transcripts"]["columns"] = list(tx_df.columns)
            audit["files"]["transcripts"]["n_rows"] = tx_df.shape[0]
            audit["files"]["transcripts"]["n_cols"] = tx_df.shape[1]
            logger.info(f"Shape: {tx_df.shape}")
            logger.info(f"Columns: {list(tx_df.columns)}")
            logger.info(f"First rows:\n{tx_df.head()}")
            
            # Identify likely coordinate and assignment columns
            coord_cols = {}
            for col in tx_df.columns:
                if "x" in col.lower():
                    coord_cols["x"] = col
                elif "y" in col.lower():
                    coord_cols["y"] = col
            
            if coord_cols:
                logger.info(f"Detected coordinate columns: {coord_cols}")
                ranges = validation.get_coordinate_ranges(
                    tx_df, coord_cols["x"], coord_cols.get("y")
                )
                logger.info(f"Coordinate ranges: {ranges}")
    
    logger.info("\n--- Polygons ---")
    poly_file = filenames.get("polygons")
    if poly_file:
        filepath = os.path.join(root_dir, poly_file)
        file_info = validation.get_file_info(filepath)
        audit["files"]["polygons"] = file_info
        
        poly_df = loader.load_polygons(poly_file)
        if poly_df is not None:
            audit["files"]["polygons"]["shape"] = poly_df.shape
            audit["files"]["polygons"]["columns"] = list(poly_df.columns)
            audit["files"]["polygons"]["n_rows"] = poly_df.shape[0]
            audit["files"]["polygons"]["n_cols"] = poly_df.shape[1]
            logger.info(f"Shape: {poly_df.shape}")
            logger.info(f"Columns: {list(poly_df.columns)}")
            logger.info(f"First rows:\n{poly_df.head()}")
    
    return audit


def audit_proseg_data(root_dir: str, filenames: dict, logger) -> dict:
    """Audit Proseg output files."""
    logger.info("\n" + "=" * 80)
    logger.info("AUDITING PROSEG DATA")
    logger.info("=" * 80)
    
    audit = {
        "root_dir": root_dir,
        "files": {},
    }
    
    loader = proseg.ProsegLoader(root_dir, logger=logger)
    
    # Audit H5AD
    logger.info("\n--- H5AD ---")
    h5ad_file = filenames.get("h5ad")
    if h5ad_file:
        filepath = os.path.join(root_dir, h5ad_file)
        file_info = validation.get_file_info(filepath)
        audit["files"]["h5ad"] = file_info
        
        adata = loader.load_h5ad(h5ad_file)
        if adata is not None:
            audit["files"]["h5ad"]["shape"] = adata.shape
            audit["files"]["h5ad"]["obs_columns"] = list(adata.obs.columns)
            audit["files"]["h5ad"]["var_columns"] = list(adata.var.columns)[:10]
            audit["files"]["h5ad"]["obsm_keys"] = list(adata.obsm.keys()) if adata.obsm else []
            audit["files"]["h5ad"]["uns_keys"] = list(adata.uns.keys()) if adata.uns else []
            audit["files"]["h5ad"]["layers"] = list(adata.layers.keys()) if adata.layers else []
            
            logger.info(f"Shape (observations x features): {adata.shape}")
            logger.info(f"Obs columns: {list(adata.obs.columns)}")
            logger.info(f"First var: {list(adata.var.columns[:5])}")
            logger.info(f"Obsm keys: {list(adata.obsm.keys())}")
            logger.info(f"Uns keys: {list(adata.uns.keys())[:5]}")
            if adata.layers:
                logger.info(f"Layers: {list(adata.layers.keys())}")
            
            # Check for coordinates
            if "spatial" in adata.obsm:
                logger.info(f"Spatial coordinates shape: {adata.obsm['spatial'].shape}")
    
    # Audit Zarr
    logger.info("\n--- Zarr ---")
    zarr_file = filenames.get("zarr")
    if zarr_file:
        filepath = os.path.join(root_dir, zarr_file)
        if os.path.isdir(filepath):
            file_info = validation.get_file_info(filepath + "/adata.h5ad")
            audit["files"]["zarr"] = {"path": filepath, "exists": True}
            
            adata = loader.load_zarr(zarr_file)
            if adata is not None:
                audit["files"]["zarr"]["shape"] = adata.shape
                audit["files"]["zarr"]["obs_columns"] = list(adata.obs.columns)
                audit["files"]["zarr"]["obsm_keys"] = list(adata.obsm.keys()) if adata.obsm else []
                
                logger.info(f"Shape: {adata.shape}")
                logger.info(f"Obs columns: {list(adata.obs.columns)}")
                logger.info(f"Obsm keys: {list(adata.obsm.keys())}")
    
    # Note on Seurat
    logger.info("\n--- Seurat RDS ---")
    seurat_file = filenames.get("seurat")
    if seurat_file:
        filepath = os.path.join(root_dir, seurat_file)
        file_info = validation.get_file_info(filepath)
        audit["files"]["seurat"] = file_info
        logger.info("Seurat RDS loading not implemented (requires R/rpy2)")
    
    return audit


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
    
    # Load configuration
    config = load_config(args.config)
    dataset_config = config["datasets"][args.dataset_id]
    
    # Audit AtoMx data
    atomx_audit = audit_atomx_data(
        dataset_config["atomx"]["root"],
        dataset_config["atomx"],
        logger,
    )
    
    # Audit Proseg data
    proseg_audit = audit_proseg_data(
        dataset_config["proseg"]["root"],
        dataset_config["proseg"],
        logger,
    )
    
    # Save audit report
    audit_report = {
        "dataset_id": args.dataset_id,
        "atomx": atomx_audit,
        "proseg": proseg_audit,
    }
    
    audit_path = os.path.join(args.output_root, "audit")
    os.makedirs(audit_path, exist_ok=True)
    
    report_file = os.path.join(audit_path, f"audit_report_{args.dataset_id}.json")
    with open(report_file, "w") as f:
        json.dump(audit_report, f, indent=2, default=str)
    
    logger.info(f"\nAudit report saved to: {report_file}")
    logger.info("=" * 80)
    logger.info("INSPECTION COMPLETE")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
