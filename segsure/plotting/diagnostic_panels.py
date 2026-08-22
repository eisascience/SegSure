"""Advanced diagnostic neighborhood panels for SegSure."""

import os
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

try:
    from shapely.geometry import Polygon
    HAS_SHAPELY = True
except ImportError:
    HAS_SHAPELY = False


class NeighborhoodDiagnosticPanel:
    """Create comprehensive 8-panel diagnostic views of neighborhoods."""

    def __init__(
        self,
        atomx_molecules: pd.DataFrame,
        proseg_molecules: pd.DataFrame,
        atomx_polygons: Optional[pd.DataFrame] = None,
        proseg_polygons: Optional[pd.DataFrame] = None,
        cell_matches: Optional[pd.DataFrame] = None,
        neighborhood_metrics: Optional[pd.DataFrame] = None,
        logger=None,
    ):
        """Initialize diagnostic panel generator.

        Parameters
        ----------
        atomx_molecules : pd.DataFrame
            AtoMx transcript/molecule data.
        proseg_molecules : pd.DataFrame
            Proseg transcript/molecule data.
        atomx_polygons : pd.DataFrame, optional
            AtoMx polygon vertices.
        proseg_polygons : pd.DataFrame, optional
            Proseg polygon vertices.
        cell_matches : pd.DataFrame, optional
            Cell matching results.
        neighborhood_metrics : pd.DataFrame, optional
            Neighborhood-level metrics.
        logger : logging.Logger, optional
            Logger instance.
        """
        self.atomx_molecules = atomx_molecules
        self.proseg_molecules = proseg_molecules
        self.atomx_polygons = atomx_polygons
        self.proseg_polygons = proseg_polygons
        self.cell_matches = cell_matches or pd.DataFrame()
        self.neighborhood_metrics = neighborhood_metrics or pd.DataFrame()
        self.logger = logger

        # Build data structures
        self._build_molecule_index()
        self._build_polygon_dict()
        self._build_cell_correspondence()

    def _build_molecule_index(self):
        """Build indexes for fast molecule lookup by cell."""
        self.atomx_molecules_by_cell = {}
        for _, row in self.atomx_molecules.iterrows():
            cell_id = row.get("atomx_cell") or row.get("cell_id")
            if cell_id:
                if cell_id not in self.atomx_molecules_by_cell:
                    self.atomx_molecules_by_cell[cell_id] = []
                self.atomx_molecules_by_cell[cell_id].append(row)

        self.proseg_molecules_by_cell = {}
        for _, row in self.proseg_molecules.iterrows():
            cell_id = row.get("proseg_cell") or row.get("cell_id")
            if cell_id:
                if cell_id not in self.proseg_molecules_by_cell:
                    self.proseg_molecules_by_cell[cell_id] = []
                self.proseg_molecules_by_cell[cell_id].append(row)

    def _build_polygon_dict(self):
        """Build polygon dictionaries for fast lookup."""
        self.atomx_poly_dict = {}
        self.proseg_poly_dict = {}

        if HAS_SHAPELY:
            if self.atomx_polygons is not None:
                self.atomx_poly_dict = self._build_poly_dict_from_df(
                    self.atomx_polygons
                )
            if self.proseg_polygons is not None:
                self.proseg_poly_dict = self._build_poly_dict_from_df(
                    self.proseg_polygons
                )

    @staticmethod
    def _build_poly_dict_from_df(
        polygon_df: pd.DataFrame,
    ) -> Dict[str, Optional[Polygon]]:
        """Build polygon dict from dataframe."""
        poly_dict = {}

        # Infer columns
        cell_col = None
        for col in ["cell_id", "cellID", "cell", "nucleus_id"]:
            if col in polygon_df.columns:
                cell_col = col
                break

        x_col = None
        y_col = None
        for col in polygon_df.columns:
            if "x" in col.lower() and x_col is None:
                x_col = col
            elif "y" in col.lower() and y_col is None:
                y_col = col

        if cell_col and x_col and y_col:
            for cell_id in polygon_df[cell_col].unique():
                vertices = list(
                    zip(
                        polygon_df[polygon_df[cell_col] == cell_id][x_col],
                        polygon_df[polygon_df[cell_col] == cell_id][y_col],
                    )
                )
                if len(vertices) >= 3:
                    try:
                        poly_dict[cell_id] = Polygon(vertices)
                    except Exception:
                        pass

        return poly_dict

    def _build_cell_correspondence(self):
        """Build bidirectional cell correspondence."""
        self.atomx_to_proseg = {}
        for _, row in self.cell_matches.iterrows():
            atomx_id = row.get("atomx_cell")
            proseg_id = row.get("proseg_cell")
            if atomx_id and proseg_id:
                self.atomx_to_proseg[atomx_id] = proseg_id

    def create_diagnostic_panel(
        self,
        atomx_cell: str,
        panel_width: float = 100.0,
        figsize: Tuple[int, int] = (16, 12),
    ) -> plt.Figure:
        """Create 8-panel diagnostic view for a neighborhood.

        Panels:
        A. AtoMx boundaries
        B. Proseg boundaries
        C. Boundary overlay
        D. Correspondence classes
        E. Molecule assignment status
        F. Fraction transcript assignment changed
        G. Geometry disagreement
        H. Neighborhood disagreement

        Parameters
        ----------
        atomx_cell : str
            Central AtoMx cell ID.
        panel_width : float, optional
            Width/height of panel region in microns.
        figsize : tuple, optional
            Figure size (inches).

        Returns
        -------
        plt.Figure
            Diagnostic figure.
        """
        proseg_cell = self.atomx_to_proseg.get(atomx_cell)

        # Get central coordinates
        atomx_mols = self.atomx_molecules_by_cell.get(atomx_cell, [])
        if not atomx_mols:
            return None

        center_x = np.mean([m.get("x", 0) for m in atomx_mols])
        center_y = np.mean([m.get("y", 0) for m in atomx_mols])

        # Get neighborhood bounds
        x_min, x_max = center_x - panel_width / 2, center_x + panel_width / 2
        y_min, y_max = center_y - panel_width / 2, center_y + panel_width / 2

        # Filter molecules to neighborhood
        atomx_neighborhood = self._filter_molecules_to_region(
            self.atomx_molecules, x_min, x_max, y_min, y_max
        )
        proseg_neighborhood = self._filter_molecules_to_region(
            self.proseg_molecules, x_min, x_max, y_min, y_max
        )

        # Create figure
        fig, axes = plt.subplots(2, 4, figsize=figsize)
        fig.suptitle(
            f"Neighborhood Diagnostic: AtoMx cell {atomx_cell} vs Proseg cell {proseg_cell}",
            fontsize=14,
            fontweight="bold",
        )

        # Panel A: AtoMx boundaries
        self._plot_atomx_boundaries(
            axes[0, 0], atomx_neighborhood, x_min, x_max, y_min, y_max
        )

        # Panel B: Proseg boundaries
        self._plot_proseg_boundaries(
            axes[0, 1], proseg_neighborhood, x_min, x_max, y_min, y_max
        )

        # Panel C: Boundary overlay
        self._plot_boundary_overlay(
            axes[0, 2],
            atomx_neighborhood,
            proseg_neighborhood,
            x_min,
            x_max,
            y_min,
            y_max,
        )

        # Panel D: Correspondence classes
        self._plot_correspondence_classes(
            axes[0, 3], atomx_neighborhood, x_min, x_max, y_min, y_max
        )

        # Panel E: Molecule assignment status
        self._plot_molecule_status(
            axes[1, 0], atomx_neighborhood, x_min, x_max, y_min, y_max
        )

        # Panel F: Fraction changed
        self._plot_fraction_changed(
            axes[1, 1], atomx_neighborhood, x_min, x_max, y_min, y_max
        )

        # Panel G: Geometry disagreement
        self._plot_geometry_disagreement(
            axes[1, 2], atomx_neighborhood, x_min, x_max, y_min, y_max
        )

        # Panel H: Neighborhood metrics
        self._plot_neighborhood_metrics(
            axes[1, 3], atomx_cell, proseg_cell
        )

        plt.tight_layout()
        return fig

    def _filter_molecules_to_region(
        self,
        molecules: pd.DataFrame,
        x_min: float,
        x_max: float,
        y_min: float,
        y_max: float,
    ) -> pd.DataFrame:
        """Filter molecules to spatial region."""
        x_col = None
        y_col = None

        for col in molecules.columns:
            if "x" in col.lower() and x_col is None:
                x_col = col
            elif "y" in col.lower() and y_col is None:
                y_col = col

        if not x_col or not y_col:
            return molecules.copy()

        return molecules[
            (molecules[x_col] >= x_min)
            & (molecules[x_col] <= x_max)
            & (molecules[y_col] >= y_min)
            & (molecules[y_col] <= y_max)
        ].copy()

    def _plot_atomx_boundaries(self, ax, molecules, x_min, x_max, y_min, y_max):
        """Plot A: AtoMx boundaries."""
        ax.set_title("A. AtoMx Boundaries", fontweight="bold")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")

        # Plot polygons if available
        if HAS_SHAPELY and self.atomx_poly_dict:
            for cell_id, poly in self.atomx_poly_dict.items():
                if poly and poly.is_valid:
                    x, y = poly.exterior.xy
                    ax.plot(x, y, "b-", linewidth=2, alpha=0.7)

        # Plot molecules
        for _, mol in molecules.iterrows():
            x, y = mol.get("x", 0), mol.get("y", 0)
            ax.scatter(x, y, c="blue", s=20, alpha=0.5)

        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")

    def _plot_proseg_boundaries(self, ax, molecules, x_min, x_max, y_min, y_max):
        """Plot B: Proseg boundaries."""
        ax.set_title("B. Proseg Boundaries", fontweight="bold")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")

        # Plot polygons if available
        if HAS_SHAPELY and self.proseg_poly_dict:
            for cell_id, poly in self.proseg_poly_dict.items():
                if poly and poly.is_valid:
                    x, y = poly.exterior.xy
                    ax.plot(x, y, "r-", linewidth=2, alpha=0.7)

        # Plot molecules
        for _, mol in molecules.iterrows():
            x, y = mol.get("x", 0), mol.get("y", 0)
            ax.scatter(x, y, c="red", s=20, alpha=0.5)

        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")

    def _plot_boundary_overlay(
        self, ax, atomx_mols, proseg_mols, x_min, x_max, y_min, y_max
    ):
        """Plot C: Boundary overlay."""
        ax.set_title("C. Boundary Overlay", fontweight="bold")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")

        # AtoMx in blue
        if HAS_SHAPELY and self.atomx_poly_dict:
            for cell_id, poly in self.atomx_poly_dict.items():
                if poly and poly.is_valid:
                    x, y = poly.exterior.xy
                    ax.plot(x, y, "b-", linewidth=2, alpha=0.7, label="AtoMx")

        # Proseg in red
        if HAS_SHAPELY and self.proseg_poly_dict:
            for cell_id, poly in self.proseg_poly_dict.items():
                if poly and poly.is_valid:
                    x, y = poly.exterior.xy
                    ax.plot(x, y, "r-", linewidth=2, alpha=0.7, label="Proseg")

        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")
        ax.legend(loc="upper right", fontsize=8)

    def _plot_correspondence_classes(
        self, ax, molecules, x_min, x_max, y_min, y_max
    ):
        """Plot D: Correspondence classes."""
        ax.set_title("D. Correspondence Classes", fontweight="bold")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")

        for _, mol in molecules.iterrows():
            x = mol.get("x", 0)
            y = mol.get("y", 0)
            cell_id = mol.get("atomx_cell") or mol.get("cell_id")

            # Color by correspondence
            if cell_id in self.atomx_to_proseg:
                color = "green"  # Has correspondence
                label = "Matched"
            else:
                color = "orange"  # No correspondence
                label = "Unmatched"

            ax.scatter(x, y, c=color, s=30, alpha=0.7)

        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")

    def _plot_molecule_status(
        self, ax, molecules, x_min, x_max, y_min, y_max
    ):
        """Plot E: Molecule assignment status."""
        ax.set_title("E. Molecule Assignment Status", fontweight="bold")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")

        # Use status column if available, otherwise use cell correspondence
        status_colors = {
            "same_matched_cell": "green",
            "changed_neighbor": "yellow",
            "split_related": "orange",
            "merge_related": "purple",
            "atomx_assigned_proseg_unassigned": "red",
            "atomx_unassigned_proseg_assigned": "blue",
            "changed_unrelated": "brown",
            "unresolved": "gray",
        }

        for _, mol in molecules.iterrows():
            x = mol.get("x", 0)
            y = mol.get("y", 0)
            status = mol.get("status", "unknown")
            color = status_colors.get(status, "gray")

            ax.scatter(x, y, c=color, s=30, alpha=0.7)

        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")

    def _plot_fraction_changed(
        self, ax, molecules, x_min, x_max, y_min, y_max
    ):
        """Plot F: Fraction transcript assignment changed."""
        ax.set_title("F. Fraction Transcript Changed", fontweight="bold")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")

        # This would be populated from transcript_assignment data
        # For now, show placeholder
        ax.text(
            0.5,
            0.5,
            "Fraction Changed\n(Per-cell)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
            color="gray",
        )

        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")

    def _plot_geometry_disagreement(
        self, ax, molecules, x_min, x_max, y_min, y_max
    ):
        """Plot G: Geometry disagreement."""
        ax.set_title("G. Geometry Disagreement (IoU)", fontweight="bold")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")

        # This would be populated from geometry metrics
        # For now, show placeholder
        ax.text(
            0.5,
            0.5,
            "Geometry Metrics\n(IoU, Area Ratio)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
            color="gray",
        )

        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")

    def _plot_neighborhood_metrics(self, ax, atomx_cell, proseg_cell):
        """Plot H: Neighborhood metrics."""
        ax.set_title("H. Neighborhood Metrics", fontweight="bold")
        ax.axis("off")

        # Get metrics if available
        metrics_text = f"""
        AtoMx Cell: {atomx_cell}
        Proseg Cell: {proseg_cell}
        
        Neighborhood Summary:
        - Shared neighbors: N/A
        - Gained neighbors: N/A
        - Lost neighbors: N/A
        - Neighbor Jaccard: N/A
        - Split density: N/A
        - Merge density: N/A
        """

        ax.text(
            0.05,
            0.95,
            metrics_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            family="monospace",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

    def save_panel(self, fig: plt.Figure, output_path: str):
        """Save diagnostic panel to file."""
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        if self.logger:
            self.logger.info(f"Saved diagnostic panel to {output_path}")
