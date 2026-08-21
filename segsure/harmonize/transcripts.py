"""Transcript assignment harmonization for SegSure."""

from typing import Dict, Optional

import pandas as pd


def get_transcript_id_column(df: pd.DataFrame, logger=None) -> Optional[str]:
    """Infer transcript ID column name from DataFrame.
    
    Looks for columns with common transcript ID naming patterns.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to inspect.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    str or None
        Inferred transcript ID column name.
    """
    transcript_names = ["transcript_id", "transcriptID", "feature", "gene", "gene_name"]
    
    for col in df.columns:
        if col.lower() in [n.lower() for n in transcript_names]:
            if logger:
                logger.info(f"Inferred transcript ID column: {col}")
            return col
    
    if logger:
        logger.warning("Could not infer transcript ID column")
    
    return None


def get_cell_assignment_column(df: pd.DataFrame, logger=None) -> Optional[str]:
    """Infer cell assignment column name from DataFrame.
    
    Looks for columns indicating which cell a transcript belongs to.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to inspect.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    str or None
        Inferred cell assignment column name.
    """
    cell_names = ["cell", "cell_id", "cellID", "assigned_cell", "nucleus_id"]
    
    for col in df.columns:
        if col.lower() in [n.lower() for n in cell_names]:
            if logger:
                logger.info(f"Inferred cell assignment column: {col}")
            return col
    
    if logger:
        logger.warning("Could not infer cell assignment column")
    
    return None


def harmonize_transcript_assignments(
    atomx_transcripts: pd.DataFrame,
    proseg_transcripts: pd.DataFrame,
    atomx_gene_col: str,
    proseg_gene_col: str,
    atomx_cell_col: str,
    proseg_cell_col: str,
    logger=None
) -> Dict[str, pd.DataFrame]:
    """Harmonize transcript assignments between AtoMx and Proseg.
    
    Parameters
    ----------
    atomx_transcripts : pd.DataFrame
        AtoMx transcript data with cell assignments.
    proseg_transcripts : pd.DataFrame
        Proseg transcript data with cell assignments.
    atomx_gene_col : str
        AtoMx gene/transcript column name.
    proseg_gene_col : str
        Proseg gene/transcript column name.
    atomx_cell_col : str
        AtoMx cell assignment column name.
    proseg_cell_col : str
        Proseg cell assignment column name.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    dict
        Dictionary with harmonized transcript assignments.
    """
    if logger:
        logger.info("Harmonizing transcript assignments")
    
    # Get gene name mappings
    atomx_genes = set(atomx_transcripts[atomx_gene_col].unique())
    proseg_genes = set(proseg_transcripts[proseg_gene_col].unique())
    
    common_genes = atomx_genes & proseg_genes
    atomx_only = atomx_genes - proseg_genes
    proseg_only = proseg_genes - atomx_genes
    
    if logger:
        logger.info(f"Common genes: {len(common_genes)}")
        logger.info(f"AtoMx-only genes: {len(atomx_only)}")
        logger.info(f"Proseg-only genes: {len(proseg_only)}")
    
    result = {
        "atomx_transcripts": atomx_transcripts,
        "proseg_transcripts": proseg_transcripts,
        "common_genes": common_genes,
        "atomx_only_genes": atomx_only,
        "proseg_only_genes": proseg_only,
    }
    
    return result


def compare_transcript_assignments(
    atomx_transcripts: pd.DataFrame,
    proseg_transcripts: pd.DataFrame,
    cell_matches: pd.DataFrame,
    atomx_cell_col: str,
    proseg_cell_col: str,
    atomx_tx_col: str,
    proseg_tx_col: str,
    logger=None
) -> pd.DataFrame:
    """Compare transcript assignments for matched cells.
    
    Parameters
    ----------
    atomx_transcripts : pd.DataFrame
        AtoMx transcript data.
    proseg_transcripts : pd.DataFrame
        Proseg transcript data.
    cell_matches : pd.DataFrame
        Cell matching results (atomx_cell, proseg_cell).
    atomx_cell_col : str
        AtoMx cell assignment column name.
    proseg_cell_col : str
        Proseg cell assignment column name.
    atomx_tx_col : str
        AtoMx transcript column name.
    proseg_tx_col : str
        Proseg transcript column name.
    logger : logging.Logger, optional
        Logger instance.
    
    Returns
    -------
    pd.DataFrame
        Comparison results with overlap and discord information.
    """
    if logger:
        logger.info("Comparing transcript assignments")
    
    comparisons = []
    
    for _, match_row in cell_matches.iterrows():
        atomx_id = match_row["atomx_cell"]
        proseg_id = match_row["proseg_cell"]
        
        atomx_txs = set(
            atomx_transcripts[
                atomx_transcripts[atomx_cell_col] == atomx_id
            ][atomx_tx_col].unique()
        )
        proseg_txs = set(
            proseg_transcripts[
                proseg_transcripts[proseg_cell_col] == proseg_id
            ][proseg_tx_col].unique()
        )
        
        overlap = atomx_txs & proseg_txs
        atomx_only = atomx_txs - proseg_txs
        proseg_only = proseg_txs - atomx_txs
        
        comparisons.append({
            "atomx_cell": atomx_id,
            "proseg_cell": proseg_id,
            "atomx_tx_count": len(atomx_txs),
            "proseg_tx_count": len(proseg_txs),
            "overlap_tx_count": len(overlap),
            "atomx_only_tx_count": len(atomx_only),
            "proseg_only_tx_count": len(proseg_only),
            "jaccard_index": len(overlap) / len(atomx_txs | proseg_txs)
            if (atomx_txs or proseg_txs)
            else 0.0,
        })
    
    return pd.DataFrame(comparisons) if comparisons else pd.DataFrame()
