"""Export functionality for SegSure disagreement analysis results."""

import json
import os
from typing import Dict, List, Optional

import pandas as pd


class DisagreementExporter:
    """Export disagreement analysis results in multiple formats."""

    def __init__(self, output_dir: str, logger=None):
        """Initialize exporter.

        Parameters
        ----------
        output_dir : str
            Output directory for results.
        logger : logging.Logger, optional
            Logger instance.
        """
        self.output_dir = output_dir
        self.logger = logger

        os.makedirs(output_dir, exist_ok=True)

    def export_molecules_parquet(
        self, molecules_df: pd.DataFrame, filename: str = "molecules_analysis.parquet"
    ) -> str:
        """Export molecule-level results to Parquet format.

        Preserves all fields:
        - transcript_id (or molecule_id)
        - gene
        - FOV
        - x, y, z coordinates
        - AtoMx cell assignment
        - Proseg cell assignment
        - match_method (stable_id or spatial_proximity)
        - status (classification)
        - spatial_distance (if spatial match)

        Parameters
        ----------
        molecules_df : pd.DataFrame
            Molecule-level data with all fields.
        filename : str, optional
            Output filename.

        Returns
        -------
        str
            Path to exported file.
        """
        output_path = os.path.join(self.output_dir, filename)

        molecules_df.to_parquet(output_path, index=False, compression="snappy")

        if self.logger:
            self.logger.info(
                f"Exported {len(molecules_df)} molecules to Parquet: {output_path}"
            )

        return output_path

    def export_molecules_csv(
        self, molecules_df: pd.DataFrame, filename: str = "molecules_analysis.csv"
    ) -> str:
        """Export molecule-level results to CSV format.

        Parameters
        ----------
        molecules_df : pd.DataFrame
            Molecule-level data.
        filename : str, optional
            Output filename.

        Returns
        -------
        str
            Path to exported file.
        """
        output_path = os.path.join(self.output_dir, filename)

        molecules_df.to_csv(output_path, index=False)

        if self.logger:
            self.logger.info(f"Exported molecules to CSV: {output_path}")

        return output_path

    def export_cell_level_summaries(
        self,
        cell_summaries_df: pd.DataFrame,
        filename: str = "cell_level_summaries.csv",
    ) -> str:
        """Export cell-level summary statistics.

        Columns:
        - cell_id
        - n_molecules
        - n_consistent
        - n_changed
        - n_unresolved
        - fraction_changed
        - per-status counts (n_same_matched_cell, n_changed_neighbor, etc.)

        Parameters
        ----------
        cell_summaries_df : pd.DataFrame
            Cell-level summary data.
        filename : str, optional
            Output filename.

        Returns
        -------
        str
            Path to exported file.
        """
        output_path = os.path.join(self.output_dir, filename)

        cell_summaries_df.to_csv(output_path, index=False)

        if self.logger:
            self.logger.info(f"Exported cell-level summaries to CSV: {output_path}")

        return output_path

    def export_geometry_metrics(
        self, geometry_df: pd.DataFrame, filename: str = "geometry_metrics.csv"
    ) -> str:
        """Export cell-pair geometry disagreement metrics.

        Columns:
        - atomx_cell, proseg_cell
        - iou (Intersection over Union)
        - coverage_atomx, coverage_proseg
        - area_difference, area_ratio
        - atomx_area, proseg_area
        - boundary_hausdorff, mean_boundary_distance

        Parameters
        ----------
        geometry_df : pd.DataFrame
            Geometry metrics by cell pair.
        filename : str, optional
            Output filename.

        Returns
        -------
        str
            Path to exported file.
        """
        output_path = os.path.join(self.output_dir, filename)

        geometry_df.to_csv(output_path, index=False)

        if self.logger:
            self.logger.info(f"Exported geometry metrics to CSV: {output_path}")

        return output_path

    def export_transcript_metrics(
        self, transcript_df: pd.DataFrame, filename: str = "transcript_metrics.csv"
    ) -> str:
        """Export cell-pair transcript assignment disagreement metrics.

        Columns:
        - atomx_cell, proseg_cell
        - atomx_tx_count, proseg_tx_count
        - overlap_tx_count
        - jaccard_index
        - discord_rate
        - (optional) discord_severity

        Parameters
        ----------
        transcript_df : pd.DataFrame
            Transcript metrics by cell pair.
        filename : str, optional
            Output filename.

        Returns
        -------
        str
            Path to exported file.
        """
        output_path = os.path.join(self.output_dir, filename)

        transcript_df.to_csv(output_path, index=False)

        if self.logger:
            self.logger.info(f"Exported transcript metrics to CSV: {output_path}")

        return output_path

    def export_neighborhood_metrics(
        self,
        neighborhood_df: pd.DataFrame,
        filename: str = "neighborhood_metrics.csv",
    ) -> str:
        """Export cell-pair neighborhood topology disagreement metrics.

        Columns:
        - atomx_cell, proseg_cell
        - n_neighbors_atomx, n_neighbors_proseg
        - shared_neighbors, gained_neighbors, lost_neighbors
        - neighbor_jaccard
        - local_split_density, local_merge_density
        - high_discord_neighbors, neighbor_discord_density

        Parameters
        ----------
        neighborhood_df : pd.DataFrame
            Neighborhood metrics by cell pair.
        filename : str, optional
            Output filename.

        Returns
        -------
        str
            Path to exported file.
        """
        output_path = os.path.join(self.output_dir, filename)

        neighborhood_df.to_csv(output_path, index=False)

        if self.logger:
            self.logger.info(f"Exported neighborhood metrics to CSV: {output_path}")

        return output_path

    def export_ambiguous_matches(
        self, ambiguous_list: List[Dict], filename: str = "ambiguous_matches.json"
    ) -> str:
        """Export ambiguous molecule matches for manual review.

        Parameters
        ----------
        ambiguous_list : list
            List of ambiguous match records.
        filename : str, optional
            Output filename.

        Returns
        -------
        str
            Path to exported file.
        """
        output_path = os.path.join(self.output_dir, filename)

        with open(output_path, "w") as f:
            json.dump(ambiguous_list, f, indent=2, default=str)

        if self.logger:
            self.logger.info(
                f"Exported {len(ambiguous_list)} ambiguous matches to: {output_path}"
            )

        return output_path

    def export_analysis_summary(
        self,
        summary: Dict,
        filename: str = "analysis_summary.json",
    ) -> str:
        """Export comprehensive analysis summary.

        Parameters
        ----------
        summary : dict
            Summary dictionary with keys:
            - fov_id
            - n_molecules_matched
            - n_molecules_unresolved
            - fraction_consistent
            - fraction_changed
            - geometry_summary
            - neighborhood_summary
            - molecule_method
            - spatial_tolerance
            - Generated outputs
        filename : str, optional
            Output filename.

        Returns
        -------
        str
            Path to exported file.
        """
        output_path = os.path.join(self.output_dir, filename)

        with open(output_path, "w") as f:
            json.dump(summary, f, indent=2, default=str)

        if self.logger:
            self.logger.info(f"Exported analysis summary to: {output_path}")

        return output_path
