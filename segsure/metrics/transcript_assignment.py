"""Transcript assignment disagreement metrics for SegSure."""

from typing import Dict, Optional

import numpy as np
import pandas as pd


def compute_transcript_overlap_metrics(
    assignment_comparisons: pd.DataFrame,
    logger=None
) -> pd.DataFrame:
    """Compute transcript overlap metrics from assignment comparisons.
    
    Parameters
    ----------
    assignment_comparisons : pd.DataFrame
        Results from compare_transcript_assignments.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Transcript overlap metrics.
    """
    if logger:
        logger.info("Computing transcript overlap metrics")
    
    results = assignment_comparisons.copy()
    
    # Compute additional metrics
    results["discord_rate"] = 1.0 - results["jaccard_index"]
    
    # Asymmetric discord (what fraction of AtoMx transcripts are in Proseg)
    total_tx = results["atomx_tx_count"] + results["proseg_tx_count"]
    results["symmetric_discord"] = (
        (results["atomx_tx_count"] + results["proseg_tx_count"] - 
         2 * results["overlap_tx_count"]) / 
        total_tx.replace(0, 1)
    )
    
    # Precision and recall
    results["atomx_recall"] = results["overlap_tx_count"] / (
        results["atomx_tx_count"].replace(0, 1)
    )
    results["proseg_recall"] = results["overlap_tx_count"] / (
        results["proseg_tx_count"].replace(0, 1)
    )
    
    if logger:
        logger.info(f"Computed metrics for {len(results)} cell pairs")
    
    return results


def compute_discord_severity(
    assignment_comparisons: pd.DataFrame,
    logger=None
) -> pd.DataFrame:
    """Classify discord severity (high, medium, low, none).
    
    Parameters
    ----------
    assignment_comparisons : pd.DataFrame
        Results from compute_transcript_overlap_metrics.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Assignment comparisons with severity classification.
    """
    if logger:
        logger.info("Classifying discord severity")
    
    results = assignment_comparisons.copy()
    
    # Define severity thresholds based on Jaccard index
    def classify_severity(jaccard):
        if jaccard >= 0.9:
            return "none"
        elif jaccard >= 0.7:
            return "low"
        elif jaccard >= 0.5:
            return "medium"
        else:
            return "high"
    
    results["discord_severity"] = results["jaccard_index"].apply(classify_severity)
    
    if logger:
        severity_counts = results["discord_severity"].value_counts()
        logger.info(f"Discord severity distribution:\n{severity_counts}")
    
    return results


def compute_gene_specific_discord(
    atomx_transcripts: pd.DataFrame,
    proseg_transcripts: pd.DataFrame,
    atomx_gene_col: str,
    proseg_gene_col: str,
    atomx_cell_col: str,
    proseg_cell_col: str,
    logger=None
) -> pd.DataFrame:
    """Compute discord metrics for each gene across cells.
    
    Parameters
    ----------
    atomx_transcripts : pd.DataFrame
        AtoMx transcript data.
    proseg_transcripts : pd.DataFrame
        Proseg transcript data.
    atomx_gene_col : str
        AtoMx gene column name.
    proseg_gene_col : str
        Proseg gene column name.
    atomx_cell_col : str
        AtoMx cell column name.
    proseg_cell_col : str
        Proseg cell column name.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Per-gene discord metrics.
    """
    if logger:
        logger.info("Computing gene-specific discord metrics")
    
    # Get all genes
    all_genes = set(atomx_transcripts[atomx_gene_col].unique()) | \
                set(proseg_transcripts[proseg_gene_col].unique())
    
    results = []
    
    for gene in all_genes:
        atomx_count = len(
            atomx_transcripts[atomx_transcripts[atomx_gene_col] == gene]
        )
        proseg_count = len(
            proseg_transcripts[proseg_transcripts[proseg_gene_col] == gene]
        )
        
        # Count assigned cells
        atomx_cells = set(
            atomx_transcripts[atomx_transcripts[atomx_gene_col] == gene][
                atomx_cell_col
            ].unique()
        )
        proseg_cells = set(
            proseg_transcripts[proseg_transcripts[proseg_gene_col] == gene][
                proseg_cell_col
            ].unique()
        )
        
        overlap_cells = atomx_cells & proseg_cells
        
        results.append({
            "gene": gene,
            "atomx_transcript_count": atomx_count,
            "proseg_transcript_count": proseg_count,
            "atomx_cell_count": len(atomx_cells),
            "proseg_cell_count": len(proseg_cells),
            "overlap_cell_count": len(overlap_cells),
            "asymmetry": abs(atomx_count - proseg_count) / (atomx_count + proseg_count + 1),
        })
    
    return pd.DataFrame(results) if results else pd.DataFrame()
