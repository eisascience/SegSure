"""Specialized loader for AtoMx cell-level data."""

import os
from typing import Optional

import pandas as pd


class AtoMxCellsLoader:
    """Loader for AtoMx cell-level data (metadata + expression)."""

    def __init__(self, root_dir: str, logger=None):
        """Initialize AtoMx cells loader.

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

    def load_cells(
        self,
        metadata_file: str,
        expression_file: str,
        fov: Optional[str] = None,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load cell metadata and expression.

        Parameters
        ----------
        metadata_file : str
            Cell metadata filename.
        expression_file : str
            Expression matrix filename.
        fov : str, optional
            FOV to filter by.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Cell metadata with optional expression columns.
        """
        # Load metadata
        meta_path = os.path.join(self.root_dir, metadata_file)
        if not os.path.isfile(meta_path):
            self._log(f"Metadata file not found: {meta_path}", "error")
            return None

        self._log(f"Loading cell metadata from: {metadata_file}")
        metadata = pd.read_csv(meta_path, index_col=0, nrows=nrows)

        # Filter by FOV if requested
        if fov and "fov" in metadata.columns:
            metadata = metadata[metadata["fov"] == fov]
            self._log(f"Filtered to FOV {fov}: {len(metadata)} cells")

        self._log(f"Loaded {len(metadata)} cells")
        return metadata

    def load_centroids(
        self,
        metadata_file: str,
        fov: Optional[str] = None,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load cell centroids/coordinates.

        Parameters
        ----------
        metadata_file : str
            Cell metadata filename (contains centroids).
        fov : str, optional
            FOV to filter by.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Cell centroids (cell_id x coordinates).
        """
        metadata = self.load_cells(
            metadata_file, "", fov=fov, nrows=nrows
        )
        if metadata is None:
            return None

        # Look for coordinate columns
        coord_cols = []
        for col in metadata.columns:
            if any(c in col.lower() for c in ["cx", "cy", "cz", "x_centroid", "y_centroid"]):
                coord_cols.append(col)

        if not coord_cols:
            self._log("No centroid columns found in metadata", "warning")
            return None

        self._log(f"Found {len(coord_cols)} coordinate columns: {coord_cols}")
        return metadata[coord_cols]
