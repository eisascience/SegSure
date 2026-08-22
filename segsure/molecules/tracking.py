"""Molecule-level transcript tracking for SegSure."""

from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


class MoleculeTracker:
    """Track and link molecule assignments across segmentation methods."""

    # Molecule status categories
    STATUS_SAME_MATCHED_CELL = "same_matched_cell"
    STATUS_CHANGED_NEIGHBOR = "changed_neighbor"
    STATUS_SPLIT_RELATED = "split_related"
    STATUS_MERGE_RELATED = "merge_related"
    STATUS_ATOMX_ASSIGNED_PROSEG_UNASSIGNED = "atomx_assigned_proseg_unassigned"
    STATUS_ATOMX_UNASSIGNED_PROSEG_ASSIGNED = "atomx_unassigned_proseg_assigned"
    STATUS_CHANGED_UNRELATED = "changed_unrelated"
    STATUS_UNRESOLVED = "unresolved"

    VALID_STATUSES = {
        STATUS_SAME_MATCHED_CELL,
        STATUS_CHANGED_NEIGHBOR,
        STATUS_SPLIT_RELATED,
        STATUS_MERGE_RELATED,
        STATUS_ATOMX_ASSIGNED_PROSEG_UNASSIGNED,
        STATUS_ATOMX_UNASSIGNED_PROSEG_ASSIGNED,
        STATUS_CHANGED_UNRELATED,
        STATUS_UNRESOLVED,
    }

    def __init__(
        self,
        atomx_molecules: pd.DataFrame,
        proseg_molecules: pd.DataFrame,
        cell_matches: pd.DataFrame,
        spatial_tolerance: float = 1.5,
        logger=None,
    ):
        """Initialize molecule tracker.

        Parameters
        ----------
        atomx_molecules : pd.DataFrame
            AtoMx transcript/molecule data with columns:
            - transcript_id or similar: unique molecule/transcript identifier
            - gene or gene_name: gene identity
            - fov: field of view ID
            - x, y, z: spatial coordinates
            - cell_id or assigned_cell: cell assignment
        proseg_molecules : pd.DataFrame
            Proseg transcript/molecule data (same schema).
        cell_matches : pd.DataFrame
            Cell matching results with atomx_cell and proseg_cell columns.
        spatial_tolerance : float, optional
            Maximum distance (microns) for spatial coordinate matching.
            Default is 1.5.
        logger : logging.Logger, optional
            Logger instance.
        """
        self.atomx_molecules = atomx_molecules.copy()
        self.proseg_molecules = proseg_molecules.copy()
        self.cell_matches = cell_matches.copy()
        self.spatial_tolerance = spatial_tolerance
        self.logger = logger

        # Detect column names
        self._detect_columns()

        # Build cell correspondence map
        self.cell_correspondence = self._build_cell_correspondence_map()

        # Track matches and ambiguities
        self.matched_molecules = []
        self.ambiguous_matches = []
        self.unmatched_atomx = []
        self.unmatched_proseg = []

    def _detect_columns(self):
        """Infer column names for molecule data."""
        # Transcript/molecule ID columns
        tx_patterns = [
            "transcript_id",
            "transcriptID",
            "molecule_id",
            "moleculeID",
            "feature",
        ]
        self.atomx_tx_col = self._find_column(self.atomx_molecules, tx_patterns)
        self.proseg_tx_col = self._find_column(self.proseg_molecules, tx_patterns)

        # Gene columns
        gene_patterns = ["gene", "gene_name", "gene_names", "feature_name"]
        self.atomx_gene_col = self._find_column(self.atomx_molecules, gene_patterns)
        self.proseg_gene_col = self._find_column(self.proseg_molecules, gene_patterns)

        # FOV columns
        fov_patterns = ["fov", "fov_id", "field_of_view"]
        self.atomx_fov_col = self._find_column(self.atomx_molecules, fov_patterns)
        self.proseg_fov_col = self._find_column(self.proseg_molecules, fov_patterns)

        # Coordinate columns
        coord_patterns_x = ["x", "x_coord"]
        coord_patterns_y = ["y", "y_coord"]
        coord_patterns_z = ["z", "z_coord"]

        self.atomx_x_col = self._find_column(self.atomx_molecules, coord_patterns_x)
        self.atomx_y_col = self._find_column(self.atomx_molecules, coord_patterns_y)
        self.atomx_z_col = self._find_column(self.atomx_molecules, coord_patterns_z)

        self.proseg_x_col = self._find_column(self.proseg_molecules, coord_patterns_x)
        self.proseg_y_col = self._find_column(self.proseg_molecules, coord_patterns_y)
        self.proseg_z_col = self._find_column(self.proseg_molecules, coord_patterns_z)

        # Cell assignment columns
        cell_patterns = ["cell_id", "cellID", "assigned_cell", "nucleus_id", "cell"]
        self.atomx_cell_col = self._find_column(self.atomx_molecules, cell_patterns)
        self.proseg_cell_col = self._find_column(self.proseg_molecules, cell_patterns)

        if self.logger:
            self.logger.info(
                f"Detected columns: "
                f"AtoMx tx={self.atomx_tx_col}, Proseg tx={self.proseg_tx_col}"
            )

    @staticmethod
    def _find_column(df: pd.DataFrame, patterns: List[str]) -> Optional[str]:
        """Find column matching any pattern (case-insensitive)."""
        for col in df.columns:
            for pattern in patterns:
                if col.lower() == pattern.lower():
                    return col
        return None

    def _build_cell_correspondence_map(self) -> Dict:
        """Build bidirectional cell correspondence mapping from cell_matches.

        Returns
        -------
        dict
            Mapping: atomx_cell_id -> proseg_cell_id and vice versa.
        """
        correspondence = {"atomx_to_proseg": {}, "proseg_to_atomx": {}}

        for _, row in self.cell_matches.iterrows():
            atomx_id = row.get("atomx_cell")
            proseg_id = row.get("proseg_cell")

            if atomx_id is not None and proseg_id is not None:
                correspondence["atomx_to_proseg"][atomx_id] = proseg_id
                correspondence["proseg_to_atomx"][proseg_id] = atomx_id

        return correspondence

    def match_molecules_by_id(self) -> pd.DataFrame:
        """Match molecules using stable IDs if available in both datasets.

        Returns
        -------
        pd.DataFrame
            Matched molecules with fields:
            - tx_id, gene, fov, x, y, z
            - atomx_cell, proseg_cell
            - match_method='stable_id'
        """
        if not self.atomx_tx_col or not self.proseg_tx_col:
            return pd.DataFrame()

        # Find common transcript IDs
        atomx_ids = set(self.atomx_molecules[self.atomx_tx_col].unique())
        proseg_ids = set(self.proseg_molecules[self.proseg_tx_col].unique())
        common_ids = atomx_ids & proseg_ids

        if not common_ids:
            return pd.DataFrame()

        matches = []
        for tx_id in common_ids:
            atomx_rows = self.atomx_molecules[
                self.atomx_molecules[self.atomx_tx_col] == tx_id
            ]
            proseg_rows = self.proseg_molecules[
                self.proseg_molecules[self.proseg_tx_col] == tx_id
            ]

            # For now, assume 1-1 per FOV; handle multiples conservatively
            if len(atomx_rows) == 1 and len(proseg_rows) == 1:
                a_row = atomx_rows.iloc[0]
                p_row = proseg_rows.iloc[0]

                matches.append(
                    {
                        "tx_id": tx_id,
                        "gene": a_row.get(self.atomx_gene_col, "unknown"),
                        "fov": a_row.get(self.atomx_fov_col, "unknown"),
                        "x": (
                            a_row.get(self.atomx_x_col)
                            if self.atomx_x_col
                            else p_row.get(self.proseg_x_col)
                        ),
                        "y": (
                            a_row.get(self.atomx_y_col)
                            if self.atomx_y_col
                            else p_row.get(self.proseg_y_col)
                        ),
                        "z": (
                            a_row.get(self.atomx_z_col)
                            if self.atomx_z_col
                            else p_row.get(self.proseg_z_col)
                        ),
                        "atomx_cell": a_row.get(self.atomx_cell_col),
                        "proseg_cell": p_row.get(self.proseg_cell_col),
                        "match_method": "stable_id",
                        "atomx_x": a_row.get(self.atomx_x_col),
                        "atomx_y": a_row.get(self.atomx_y_col),
                        "atomx_z": a_row.get(self.atomx_z_col),
                        "proseg_x": p_row.get(self.proseg_x_col),
                        "proseg_y": p_row.get(self.proseg_y_col),
                        "proseg_z": p_row.get(self.proseg_z_col),
                    }
                )

        result_df = pd.DataFrame(matches)
        if self.logger:
            self.logger.info(f"Matched {len(result_df)} molecules by stable ID")

        return result_df

    def match_molecules_by_spatial_proximity(
        self, exclude_matched_ids: Optional[Set[str]] = None
    ) -> Tuple[pd.DataFrame, List[Dict]]:
        """Match unmatched molecules using gene identity + spatial coordinates.

        Parameters
        ----------
        exclude_matched_ids : set, optional
            Set of transcript IDs already matched by stable ID.

        Returns
        -------
        tuple
            - pd.DataFrame: Spatially matched molecules
            - list: Ambiguous matches (multiple candidates within tolerance)
        """
        if exclude_matched_ids is None:
            exclude_matched_ids = set()

        # Filter to unmatched molecules
        atomx_unmatched = self.atomx_molecules[
            ~self.atomx_molecules[self.atomx_tx_col].isin(exclude_matched_ids)
        ].copy()
        proseg_unmatched = self.proseg_molecules[
            ~self.proseg_molecules[self.proseg_tx_col].isin(exclude_matched_ids)
        ].copy()

        if atomx_unmatched.empty or proseg_unmatched.empty:
            return pd.DataFrame(), []

        matches = []
        ambiguous = []
        matched_proseg_indices = set()

        # Group by gene for efficiency
        for gene in atomx_unmatched[self.atomx_gene_col].unique():
            atomx_gene_group = atomx_unmatched[
                atomx_unmatched[self.atomx_gene_col] == gene
            ]
            proseg_gene_group = proseg_unmatched[
                proseg_unmatched[self.proseg_gene_col] == gene
            ]

            if atomx_gene_group.empty or proseg_gene_group.empty:
                continue

            # Extract coordinates
            atomx_coords = atomx_gene_group[
                [self.atomx_x_col, self.atomx_y_col]
            ].values
            proseg_coords = proseg_gene_group[
                [self.proseg_x_col, self.proseg_y_col]
            ].values

            # Compute pairwise distances
            distances = cdist(atomx_coords, proseg_coords, metric="euclidean")

            # Find matches within tolerance
            for a_idx, (a_dist_row) in enumerate(distances):
                candidates = np.where(a_dist_row <= self.spatial_tolerance)[0]

                if len(candidates) == 0:
                    # No match found
                    self.unmatched_atomx.append(
                        atomx_gene_group.iloc[a_idx].to_dict()
                    )
                elif len(candidates) == 1:
                    # Unambiguous match
                    p_idx = candidates[0]
                    if p_idx not in matched_proseg_indices:
                        a_row = atomx_gene_group.iloc[a_idx]
                        p_row = proseg_gene_group.iloc[p_idx]

                        matches.append(
                            {
                                "tx_id": (
                                    a_row.get(self.atomx_tx_col, "unknown")
                                    if self.atomx_tx_col
                                    else f"atomx_{a_idx}"
                                ),
                                "gene": gene,
                                "fov": a_row.get(self.atomx_fov_col, "unknown"),
                                "x": a_row.get(self.atomx_x_col),
                                "y": a_row.get(self.atomx_y_col),
                                "z": a_row.get(self.atomx_z_col),
                                "atomx_cell": a_row.get(self.atomx_cell_col),
                                "proseg_cell": p_row.get(self.proseg_cell_col),
                                "spatial_distance": a_dist_row[p_idx],
                                "match_method": "spatial_proximity",
                                "atomx_x": a_row.get(self.atomx_x_col),
                                "atomx_y": a_row.get(self.atomx_y_col),
                                "atomx_z": a_row.get(self.atomx_z_col),
                                "proseg_x": p_row.get(self.proseg_x_col),
                                "proseg_y": p_row.get(self.proseg_y_col),
                                "proseg_z": p_row.get(self.proseg_z_col),
                            }
                        )
                        matched_proseg_indices.add(p_idx)
                else:
                    # Ambiguous: multiple candidates
                    a_row = atomx_gene_group.iloc[a_idx]
                    ambig_candidates = [
                        {
                            "proseg_idx": p_idx,
                            "proseg_cell": proseg_gene_group.iloc[p_idx].get(
                                self.proseg_cell_col
                            ),
                            "distance": a_dist_row[p_idx],
                        }
                        for p_idx in candidates
                        if p_idx not in matched_proseg_indices
                    ]

                    if ambig_candidates:
                        ambiguous.append(
                            {
                                "atomx_tx": a_row.get(self.atomx_tx_col),
                                "atomx_cell": a_row.get(self.atomx_cell_col),
                                "gene": gene,
                                "x": a_row.get(self.atomx_x_col),
                                "y": a_row.get(self.atomx_y_col),
                                "candidates": ambig_candidates,
                            }
                        )

        result_df = pd.DataFrame(matches)
        if self.logger:
            self.logger.info(
                f"Matched {len(result_df)} molecules by spatial proximity"
            )
            self.logger.warning(f"Found {len(ambiguous)} ambiguous matches")

        return result_df, ambiguous

    def classify_molecules(
        self, matched_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Classify each matched molecule into one of 8 status categories.

        Parameters
        ----------
        matched_df : pd.DataFrame
            Output from match_molecules_by_id or match_molecules_by_spatial_proximity.

        Returns
        -------
        pd.DataFrame
            Input dataframe with added 'status' column.
        """
        results = matched_df.copy()
        statuses = []

        for _, row in results.iterrows():
            atomx_cell = row.get("atomx_cell")
            proseg_cell = row.get("proseg_cell")

            # Determine if cells correspond via validated cell matching
            same_cell_correspondence = (
                self.cell_correspondence["atomx_to_proseg"].get(atomx_cell)
                == proseg_cell
            )

            # Check if cells are neighbors
            is_neighbor = False
            if atomx_cell and proseg_cell:
                # Determine neighbors (cells within small distance)
                # This requires cell coordinate data; for now, check via cell_matches
                is_neighbor = self._are_cells_adjacent(atomx_cell, proseg_cell)

            # Classify status
            if pd.isna(atomx_cell) and not pd.isna(proseg_cell):
                status = self.STATUS_ATOMX_UNASSIGNED_PROSEG_ASSIGNED
            elif not pd.isna(atomx_cell) and pd.isna(proseg_cell):
                status = self.STATUS_ATOMX_ASSIGNED_PROSEG_UNASSIGNED
            elif same_cell_correspondence:
                status = self.STATUS_SAME_MATCHED_CELL
            elif is_neighbor:
                status = self.STATUS_CHANGED_NEIGHBOR
            elif self._is_split_related(atomx_cell):
                status = self.STATUS_SPLIT_RELATED
            elif self._is_merge_related(proseg_cell):
                status = self.STATUS_MERGE_RELATED
            else:
                status = self.STATUS_CHANGED_UNRELATED

            statuses.append(status)

        results["status"] = statuses
        return results

    def _are_cells_adjacent(self, atomx_cell, proseg_cell) -> bool:
        """Check if two cells are adjacent (neighbors).

        Currently uses cell matching data; could be enhanced with
        actual coordinate-based distance checking.
        """
        # Check if either cell has a relationship to the other
        # This is a simplification; full implementation would use spatial distance
        for _, match_row in self.cell_matches.iterrows():
            if (
                match_row.get("atomx_cell") == atomx_cell
                and match_row.get("proseg_cell") != proseg_cell
            ):
                # atomx_cell has a primary match to a different proseg cell
                return False

        return False  # Conservative: assume not adjacent if no direct evidence

    def _is_split_related(self, atomx_cell) -> bool:
        """Check if atomx_cell is involved in a split relationship."""
        # A split occurs when one AtoMx cell maps to multiple Proseg cells
        matches_for_cell = self.cell_matches[
            self.cell_matches["atomx_cell"] == atomx_cell
        ]
        return len(matches_for_cell) > 1

    def _is_merge_related(self, proseg_cell) -> bool:
        """Check if proseg_cell is involved in a merge relationship."""
        # A merge occurs when multiple AtoMx cells map to one Proseg cell
        matches_for_cell = self.cell_matches[
            self.cell_matches["proseg_cell"] == proseg_cell
        ]
        return len(matches_for_cell) > 1

    def process_all_molecules(self) -> Tuple[pd.DataFrame, List[Dict]]:
        """Perform complete molecule matching and classification.

        Returns
        -------
        tuple
            - pd.DataFrame: All matched molecules with status classifications
            - list: Ambiguous match candidates to review manually
        """
        # Match by stable ID
        stable_matched = self.match_molecules_by_id()
        stable_ids = (
            set(stable_matched["tx_id"].unique())
            if not stable_matched.empty
            else set()
        )

        # Match by spatial proximity
        spatial_matched, ambiguous = self.match_molecules_by_spatial_proximity(
            exclude_matched_ids=stable_ids
        )

        # Combine
        all_matched = pd.concat(
            [stable_matched, spatial_matched], ignore_index=True
        )

        # Classify
        if not all_matched.empty:
            all_matched = self.classify_molecules(all_matched)

        if self.logger:
            if not all_matched.empty:
                status_counts = all_matched["status"].value_counts()
                self.logger.info(f"Molecule status distribution:\n{status_counts}")

        return all_matched, ambiguous
