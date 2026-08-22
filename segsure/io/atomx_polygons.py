"""Specialized loader for AtoMx polygon/segmentation data."""

import os
from typing import Dict, List, Optional

import pandas as pd


class AtoMxPolygonsLoader:
    """Loader for AtoMx polygon/segmentation data."""

    def __init__(self, root_dir: str, logger=None):
        """Initialize AtoMx polygons loader.

        Parameters
        ----------
        root_dir : str
            Root directory containing AtoMx export files.
        logger : logging.Logger, optional
            Logger instance.
        """
        self.root_dir = root_dir
        self.logger = logger

    def _log(self, msg: str, level: str = "info"):
        """Helper for logging."""
        if self.logger:
            getattr(self.logger, level)(msg)

    def load_polygons(
        self,
        polygons_file: str,
        fov: Optional[str] = None,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load polygon/segmentation data.

        Parameters
        ----------
        polygons_file : str
            Polygons filename.
        fov : str, optional
            FOV to filter by.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Polygon vertex data.
        """
        filepath = os.path.join(self.root_dir, polygons_file)
        if not os.path.isfile(filepath):
            self._log(f"Polygons file not found: {filepath}", "error")
            return None

        self._log(f"Loading polygons from: {polygons_file}")
        polygons = pd.read_csv(filepath, nrows=nrows)

        if fov and "fov" in polygons.columns:
            polygons = polygons[polygons["fov"] == fov]
            self._log(f"Filtered to FOV {fov}: {len(polygons)} polygon vertices")

        self._log(f"Loaded {len(polygons)} polygon vertices")
        return polygons

    def get_unique_cells(
        self,
        polygons_file: str,
        fov: Optional[str] = None,
    ) -> Optional[List[str]]:
        """Get list of unique cell identifiers in polygons.

        Parameters
        ----------
        polygons_file : str
            Polygons filename.
        fov : str, optional
            FOV to filter by.

        Returns
        -------
        list or None
            Unique cell identifiers.
        """
        polygons = self.load_polygons(polygons_file, fov=fov)
        if polygons is None:
            return None

        # Look for cell ID columns
        cell_cols = [
            col for col in polygons.columns
            if any(c in col.lower() for c in ["cell", "nucleus", "segment", "id"])
        ]

        if not cell_cols:
            self._log("No cell ID columns found in polygons", "warning")
            return None

        cell_col = cell_cols[0]
        unique_cells = sorted(polygons[cell_col].unique().tolist())
        self._log(f"Found {len(unique_cells)} unique cells")
        return unique_cells

    def load_cell_boundary(
        self,
        polygons_file: str,
        cell_id: str,
        fov: Optional[str] = None,
    ) -> Optional[pd.DataFrame]:
        """Load polygon vertices for a specific cell.

        Parameters
        ----------
        polygons_file : str
            Polygons filename.
        cell_id : str
            Cell identifier.
        fov : str, optional
            FOV to filter by.

        Returns
        -------
        pd.DataFrame or None
            Polygon vertices for the cell.
        """
        polygons = self.load_polygons(polygons_file, fov=fov)
        if polygons is None:
            return None

        # Find cell ID column
        cell_cols = [
            col for col in polygons.columns
            if any(c in col.lower() for c in ["cell", "nucleus", "segment"])
        ]

        if not cell_cols:
            self._log("No cell ID columns found", "warning")
            return None

        cell_col = cell_cols[0]
        cell_polygons = polygons[polygons[cell_col] == cell_id]

        if cell_polygons.empty:
            self._log(f"No polygons found for cell {cell_id}", "warning")
            return None

        self._log(f"Found {len(cell_polygons)} vertices for cell {cell_id}")
        return cell_polygons

    def load_polygon_coordinates(
        self,
        polygons_file: str,
        fov: Optional[str] = None,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load polygon vertex coordinates only.

        Parameters
        ----------
        polygons_file : str
            Polygons filename.
        fov : str, optional
            FOV to filter by.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Polygon coordinates.
        """
        polygons = self.load_polygons(polygons_file, fov=fov, nrows=nrows)
        if polygons is None:
            return None

        # Find coordinate columns
        coord_cols = [
            col for col in polygons.columns
            if any(c in col.lower() for c in ["x", "y", "z", "coord"])
        ]

        if not coord_cols:
            self._log("No coordinate columns found", "warning")
            return None

        self._log(f"Found {len(coord_cols)} coordinate columns")
        return polygons[coord_cols]
