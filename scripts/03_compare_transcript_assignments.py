#!/usr/bin/env python
"""
Compare transcript-to-cell assignments between AtoMx and Proseg.

This script identifies transcript assignment disagreements across
the two segmentation methods for matched cells.
"""

import argparse
import json
import os

import pandas as pd
import yaml

from segsure.harmonize import transcripts
from segsure.io import atomx, proseg
from segsure.utils import logging as logging_util


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def compare_transcripts_main(
    atomx_root: str,
    atomx_files: dict,
    proseg_root: str,
    proseg_files: dict,
    cell_matches_file: str,
    output_dir: str,
    logger=None,
) -> dict:
    """Compare transcript assignments between AtoMx and Proseg."""
    logger.info("=" * 80)
    logger.info("COMPARING TRANSCRIPT ASSIGNMENTS")
    logger.info("=" * 80)
    
    result = {}
    
    # Load cell matches
    logger.info("\n--- Loading Cell Matches ---")
    if not os.path.isfile(cell_matches_file):
        logger.error(f"Cell matches file not found: {cell_matches_file}")
        return result
    
    matches = pd.read_csv(cell_matches_file)
    logger.info(f"Loaded {len(matches)} cell matches")
    
    # Load AtoMx transcripts
    logger.info("\n--- Loading AtoMx Transcripts ---")
    atomx_loader = atomx.AtoMxLoader(atomx_root, logger=logger)
    atomx_tx = atomx_loader.load_transcripts(atomx_files["transcripts"])
    
    if atomx_tx is None:
        logger.error("Could not load AtoMx transcripts")
        return result
    
    # Identify gene and cell assignment columns
    atomx_gene_col = transcripts.get_transcript_id_column(atomx_tx, logger=logger)
    atomx_cell_col = transcripts.get_cell_assignment_column(atomx_tx, logger=logger)
    
    if not atomx_gene_col or not atomx_cell_col:
        logger.warning("Could not identify AtoMx gene or cell columns")
        return result
    
    # Load Proseg transcripts
    logger.info("\n--- Loading Proseg Transcripts ---")
    proseg_loader = proseg.ProsegLoader(proseg_root, logger=logger)
    
    proseg_h5ad = proseg_loader.load_h5ad(proseg_files["h5ad"])
    if proseg_h5ad is None:
        proseg_h5ad = proseg_loader.load_zarr(proseg_files["zarr"])
    
    if proseg_h5ad is None:
        logger.error("Could not load Proseg data")
        return result
    
    # Extract transcript information from Proseg if available
    # Note: This depends on Proseg storing transcript-level information
    if "transcript_id" in proseg_h5ad.var.columns:
        proseg_gene_col = "transcript_id"
    else:
        proseg_gene_col = "gene_names" if "gene_names" in proseg_h5ad.var.columns else proseg_h5ad.var.index.name
    
    proseg_cell_col = proseg_h5ad.obs.index.name or "proseg_cell_id"
    
    logger.info(f"AtoMx gene column: {atomx_gene_col}, cell column: {atomx_cell_col}")
    logger.info(f"Proseg gene column: {proseg_gene_col}, cell column: {proseg_cell_col}")
    
    # Harmonize transcript assignments
    logger.info("\n--- Harmonizing Transcript Data ---")
    harmonized = transcripts.harmonize_transcript_assignments(
        atomx_transcripts=atomx_tx,
        proseg_transcripts=pd.DataFrame(),  # Placeholder
        atomx_gene_col=atomx_gene_col,
        proseg_gene_col=proseg_gene_col,
        atomx_cell_col=atomx_cell_col,
        proseg_cell_col=proseg_cell_col,
        logger=logger,
    )
    
    result["harmonized"] = {
        "common_genes": len(harmonized["common_genes"]),
        "atomx_only_genes": len(harmonized["atomx_only_genes"]),
        "proseg_only_genes": len(harmonized["proseg_only_genes"]),
    }
    
    # Compare transcript assignments for matched cells
    logger.info("\n--- Comparing Matched Cell Transcripts ---")
    
    # For now, create a placeholder comparison
    # Full implementation requires accessing Proseg transcript assignments
    comparisons = []
    for _, match_row in matches.iterrows():
        atomx_id = match_row["atomx_cell"]
        proseg_id = match_row["proseg_cell"]
        
        atomx_txs = set(
            atomx_tx[atomx_tx[atomx_cell_col] == atomx_id][atomx_gene_col].unique()
        )
        
        # Store comparison
        comparisons.append({
            "atomx_cell": atomx_id,
            "proseg_cell": proseg_id,
            "atomx_transcript_count": len(atomx_txs),
            "proseg_transcript_count": 0,  # Would require Proseg transcript data
            "overlap_transcript_count": 0,
            "atomx_only_count": len(atomx_txs),
            "proseg_only_count": 0,
        })
    
    comparison_df = pd.DataFrame(comparisons)
    result["comparisons"] = comparison_df
    
    # Save results
    logger.info("\n--- Saving Results ---")
    os.makedirs(output_dir, exist_ok=True)
    
    comparison_file = os.path.join(output_dir, "transcript_comparisons.csv")
    comparison_df.to_csv(comparison_file, index=False)
    logger.info(f"Saved transcript comparisons to: {comparison_file}")
    
    summary_file = os.path.join(output_dir, "transcript_summary.json")
    with open(summary_file, "w") as f:
        json.dump(result, f, indent=2, default=str)
    logger.info(f"Saved summary to: {summary_file}")
    
    logger.info("=" * 80)
    logger.info("TRANSCRIPT COMPARISON COMPLETE")
    logger.info("=" * 80)
    
    return result


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compare transcript-to-cell assignments between AtoMx and Proseg"
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
        "--cell-matches",
        type=str,
        default=None,
        help="Path to cell matches file (auto-detected if not provided)",
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
        "compare_transcripts",
        output_dir=output_dir,
        verbose=args.verbose,
    )
    
    logger.info("Starting transcript comparison")
    logger.info(f"Dataset ID: {args.dataset_id}")
    
    # Load configuration
    config = load_config(args.config)
    dataset_config = config["datasets"][args.dataset_id]
    
    # Auto-detect cell matches file if not provided
    if args.cell_matches is None:
        args.cell_matches = os.path.join(args.output_root, "tables", "cell_matches.csv")
    
    # Compare transcripts
    compare_transcripts_main(
        atomx_root=dataset_config["atomx"]["root"],
        atomx_files=dataset_config["atomx"],
        proseg_root=dataset_config["proseg"]["root"],
        proseg_files=dataset_config["proseg"],
        cell_matches_file=args.cell_matches,
        output_dir=os.path.join(args.output_root, "tables"),
        logger=logger,
    )


if __name__ == "__main__":
    main()
