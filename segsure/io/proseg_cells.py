"""Specialized loader for Proseg cell-level data."""

import os
from typing import Optional

import anndata as ad
import pandas as pd


class ProsegCellsLoader:
    """Loader for Proseg cell-level data."""

    def __init__(self, root_dir: str, logger=None):
        """Initialize Proseg cells loader.

        Parameters
        ----------
        root_dir : str
            Root directory containing Proseg output.
        logger : logging.Logger, optional
            Logger instance.
        """
        self.root_dir = root_dir
        self.logger = logger

    def _log(self, msg: str, level: str = "info"):
        """Helper for logging."""
        if self.logger:
            getattr(self.logger, level)(msg)

    def load_cells_from_zarr(
        self,
        zarr_path: str,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load cell metadata from Proseg Zarr output.

        Parameters
        ----------
        zarr_path : str
            Zarr directory path.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Cell metadata.
        """
        full_path = os.path.join(self.root_dir, zarr_path)
        if not os.path.isdir(full_path):
            self._log(f"Zarr directory not found: {full_path}", "error")
            return None

        try:
            self._log(f"Loading cells from Zarr: {zarr_path}")
            adata = ad.read_zarr(full_path)

            # Extract obs (observations = cells)
            cells = adata.obs.copy()
            if nrows:
                cells = cells.iloc[:nrows]

            self._log(f"Loaded {len(cells)} cells")
            return cells
        except Exception as e:
            self._log(f"Error loading cells from Zarr: {e}", "error")
            return None

    def load_cells_from_h5ad(
        self,
        h5ad_path: str,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load cell metadata from Proseg H5AD output.

        Parameters
        ----------
        h5ad_path : str
            H5AD file path.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Cell metadata.
        """
        full_path = os.path.join(self.root_dir, h5ad_path)
        if not os.path.isfile(full_path):
            self._log(f"H5AD file not found: {full_path}", "error")
            return None

        try:
            self._log(f"Loading cells from H5AD: {h5ad_path}")
            adata = ad.read_h5ad(full_path)

            # Extract obs (observations = cells)
            cells = adata.obs.copy()
            if nrows:
                cells = cells.iloc[:nrows]

            self._log(f"Loaded {len(cells)} cells")
            return cells
        except Exception as e:
            self._log(f"Error loading cells from H5AD: {e}", "error")
            return None

    def load_centroids_from_zarr(
        self,
        zarr_path: str,
    ) -> Optional[pd.DataFrame]:
        """Load cell centroids from Proseg Zarr obsm.

        Parameters
        ----------
        zarr_path : str
            Zarr directory path.

        Returns
        -------
        pd.DataFrame or None
            Cell centroids with shape (n_cells, n_coords).
        """
        full_path = os.path.join(self.root_dir, zarr_path)
        if not os.path.isdir(full_path):
            self._log(f"Zarr directory not found: {full_path}", "error")
            return None

        try:
            self._log(f"Loading centroids from Zarr: {zarr_path}")
            adata = ad.read_zarr(full_path)

            if "centroids" in adata.obsm:
                centroids = pd.DataFrame(
                    adata.obsm["centroids"],
                    index=adata.obs_names,
                    columns=[f"coord_{i}" for i in range(adata.obsm["centroids"].shape[1])]
                )
                self._log(f"Loaded {len(centroids)} centroids")
                return centroids
            elif "spatial" in adata.obsm:
                spatial = pd.DataFrame(
                    adata.obsm["spatial"],
                    index=adata.obs_names,
                    columns=[f"spatial_{i}" for i in range(adata.obsm["spatial"].shape[1])]
                )
                self._log(f"Loaded {len(spatial)} spatial coordinates")
                return spatial
            else:
                self._log("No centroids or spatial coordinates found", "warning")
                return None
        except Exception as e:
            self._log(f"Error loading centroids: {e}", "error")
            return None
