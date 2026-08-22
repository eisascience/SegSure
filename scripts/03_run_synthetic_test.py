#!/usr/bin/env python
"""
Generate synthetic cell correspondence test data and run the matching pipeline.

This script creates synthetic polygon data to demonstrate the cell matching algorithm
and generates diagnostic figures.
"""

import os
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path

from shapely.geometry import box, Polygon
from shapely.wkt import dumps as wkt_dumps

from segsure.harmonize import cells
from segsure.plotting import cell_correspondence
from segsure.utils import logging as logging_util


def create_synthetic_test_data(output_dir: str, n_one_to_one: int = 10, 
                               n_splits: int = 3, n_merges: int = 2, logger=None):
    """Create synthetic cell data for testing.
    
    Parameters
    ----------
    output_dir : str
        Output directory for test data.
    n_one_to_one : int
        Number of one-to-one matching cells.
    n_splits : int
        Number of splitting events.
    n_merges : int
        Number of merging events.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    dict
        Test data with atomx_df, proseg_df, atomx_polygons, proseg_polygons
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if logger:
        logger.info(f"Creating synthetic test data in {output_dir}")
        logger.info(f"  One-to-one: {n_one_to_one}")
        logger.info(f"  Splits: {n_splits}")
        logger.info(f"  Merges: {n_merges}")
    
    atomx_rows = []
    atomx_polys = {}
    proseg_rows = []
    proseg_polys = {}
    
    cell_id = 0
    x_offset = 0
    
    # Create one-to-one matching cells
    for i in range(n_one_to_one):
        atomx_id = f"a{cell_id}"
        proseg_id = f"p{cell_id}"
        
        # Create overlapping polygons
        poly_a = box(x_offset, 0, x_offset + 10, 10)
        poly_p = box(x_offset + 0.5, 0.5, x_offset + 10.5, 10.5)
        
        atomx_rows.append({
            "cell_id": atomx_id,
            "x": x_offset + 5,
            "y": 5,
        })
        atomx_polys[atomx_id] = poly_a
        
        proseg_rows.append({
            "cell_id": proseg_id,
            "x": x_offset + 5.5,
            "y": 5.5,
        })
        proseg_polys[proseg_id] = poly_p
        
        cell_id += 1
        x_offset += 15
    
    # Create split events (1 AtoMx -> 2 Proseg)
    for i in range(n_splits):
        atomx_id = f"a{cell_id}"
        
        # Large AtoMx cell
        poly_a = box(x_offset, 0, x_offset + 20, 10)
        atomx_rows.append({
            "cell_id": atomx_id,
            "x": x_offset + 10,
            "y": 5,
        })
        atomx_polys[atomx_id] = poly_a
        
        # Two smaller Proseg cells
        for j in range(2):
            proseg_id = f"p{cell_id}"
            poly_p = box(x_offset + j * 10, 0, x_offset + (j + 1) * 10, 10)
            proseg_rows.append({
                "cell_id": proseg_id,
                "x": x_offset + j * 10 + 5,
                "y": 5,
            })
            proseg_polys[proseg_id] = poly_p
            cell_id += 1
        
        x_offset += 25
    
    # Create merge events (2 AtoMx -> 1 Proseg)
    for i in range(n_merges):
        # Two AtoMx cells
        for j in range(2):
            atomx_id = f"a{cell_id}"
            poly_a = box(x_offset + j * 10, 0, x_offset + (j + 1) * 10, 10)
            atomx_rows.append({
                "cell_id": atomx_id,
                "x": x_offset + j * 10 + 5,
                "y": 5,
            })
            atomx_polys[atomx_id] = poly_a
            cell_id += 1
        
        # One large Proseg cell
        proseg_id = f"p{cell_id}"
        poly_p = box(x_offset, 0, x_offset + 20, 10)
        proseg_rows.append({
            "cell_id": proseg_id,
            "x": x_offset + 10,
            "y": 5,
        })
        proseg_polys[proseg_id] = poly_p
        cell_id += 1
        x_offset += 25
    
    # Add a lost cell (AtoMx with no match)
    atomx_id = f"a{cell_id}"
    poly_a = box(x_offset + 100, 0, x_offset + 110, 10)
    atomx_rows.append({
        "cell_id": atomx_id,
        "x": x_offset + 105,
        "y": 5,
    })
    atomx_polys[atomx_id] = poly_a
    cell_id += 1
    
    # Add a new cell (Proseg with no match)
    proseg_id = f"p{cell_id}"
    poly_p = box(x_offset + 100, 0, x_offset + 110, 10)
    proseg_rows.append({
        "cell_id": proseg_id,
        "x": x_offset + 105,
        "y": 5,
    })
    proseg_polys[proseg_id] = poly_p
    
    atomx_df = pd.DataFrame(atomx_rows)
    proseg_df = pd.DataFrame(proseg_rows)
    
    if logger:
        logger.info(f"Created {len(atomx_df)} AtoMx cells and {len(proseg_df)} Proseg cells")
    
    return {
        "atomx_df": atomx_df,
        "proseg_df": proseg_df,
        "atomx_polygons": atomx_polys,
        "proseg_polygons": proseg_polys,
    }


def run_synthetic_test(output_dir: str = None, logger=None):
    """Run synthetic test of cell matching pipeline.
    
    Parameters
    ----------
    output_dir : str, optional
        Output directory for results.
    logger : logging.Logger, optional
        Logger instance.
    """
    if output_dir is None:
        output_dir = os.path.join("results", "synthetic_test")
    
    if logger:
        logger.info("=" * 80)
        logger.info("RUNNING SYNTHETIC CELL MATCHING TEST")
        logger.info("=" * 80)
    
    # Create synthetic data
    test_data = create_synthetic_test_data(output_dir, logger=logger)
    atomx_df = test_data["atomx_df"]
    proseg_df = test_data["proseg_df"]
    atomx_polys = test_data["atomx_polygons"]
    proseg_polys = test_data["proseg_polygons"]
    
    # Run matching
    if logger:
        logger.info("\n--- Running Cell Matching ---")
    
    results = cells.match_cells(
        atomx_cells=atomx_df,
        proseg_cells=proseg_df,
        atomx_coords={"x": "x", "y": "y"},
        proseg_coords={"x": "x", "y": "y"},
        atomx_cell_id="cell_id",
        proseg_cell_id="cell_id",
        atomx_polygons=atomx_polys,
        proseg_polygons=proseg_polys,
        logger=logger,
    )
    
    # Print summary
    if logger:
        logger.info("\n--- Correspondence Summary ---")
        logger.info(f"Total correspondences: {len(results)}")
        
        if "relationship_type" in results.columns:
            logger.info("\nRelationship Type Breakdown:")
            for rel_type, count in results["relationship_type"].value_counts().items():
                logger.info(f"  {rel_type}: {count}")
            
            one_to_one = results[results["relationship_type"] == "one_to_one"]
            if len(one_to_one) > 0 and "iou" in one_to_one.columns:
                logger.info(f"\nOne-to-One Statistics:")
                logger.info(f"  Count: {len(one_to_one)}")
                logger.info(f"  Median IoU: {one_to_one['iou'].median():.4f}")
                logger.info(f"  Mean IoU: {one_to_one['iou'].mean():.4f}")
                logger.info(f"  Min IoU: {one_to_one['iou'].min():.4f}")
                logger.info(f"  Max IoU: {one_to_one['iou'].max():.4f}")
    
    # Create visualizations
    if logger:
        logger.info("\n--- Creating Visualizations ---")
    
    fig_path = os.path.join(output_dir, "correspondence_figure.png")
    cell_correspondence.create_correspondence_figure(
        atomx_polys, proseg_polys, results,
        output_path=fig_path,
        logger=logger
    )
    
    summary_path = os.path.join(output_dir, "relationship_summary.png")
    cell_correspondence.create_relationship_summary_figure(
        results,
        output_path=summary_path,
        logger=logger
    )
    
    # Save results
    if logger:
        logger.info("\n--- Saving Results ---")
    
    csv_path = os.path.join(output_dir, "synthetic_matches.csv.gz")
    results.to_csv(csv_path, index=False, compression="gzip")
    if logger:
        logger.info(f"Saved CSV to: {csv_path}")
    
    parquet_path = os.path.join(output_dir, "synthetic_matches.parquet")
    results.to_parquet(parquet_path, index=False)
    if logger:
        logger.info(f"Saved Parquet to: {parquet_path}")
    
    if logger:
        logger.info("=" * 80)
        logger.info("SYNTHETIC TEST COMPLETE")
        logger.info("=" * 80)
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run synthetic cell matching test")
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/synthetic_test",
        help="Output directory",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    
    args = parser.parse_args()
    
    # Setup logging
    os.makedirs(args.output_dir, exist_ok=True)
    logger = logging_util.setup_logger(
        "synthetic_test",
        output_dir=os.path.join(args.output_dir, "logs"),
        verbose=args.verbose,
    )
    
    run_synthetic_test(output_dir=args.output_dir, logger=logger)
