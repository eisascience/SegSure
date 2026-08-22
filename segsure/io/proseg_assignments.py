"""Specialized loaders for Proseg polygons and transcript assignments."""

import json
import os
from typing import Dict, List, Optional

import anndata as ad
import pandas as pd
import zarr


class ProsegPolygonsLoader:
    """Loader for Proseg polygon/segmentation data."""

    def __init__(self, root_dir: str, logger=None):
        """Initialize Proseg polygons loader.

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

    def inspect_polygon_storage(self, zarr_path: str) -> Optional[Dict]:
        """Inspect where polygons are stored in Zarr.

        Parameters
        ----------
        zarr_path : str
            Zarr directory path.

        Returns
        -------
        dict or None
            Information about polygon storage.
        """
        full_path = os.path.join(self.root_dir, zarr_path)
        if not os.path.isdir(full_path):
            return None

        try:
            root = zarr.open_group(full_path, mode='r')
            storage_info = {
                "found_in": [],
                "details": {}
            }

            # Check obsm for polygon data
            if 'obsm' in root:
                obsm_group = root['obsm']
                polygon_keys = [
                    k for k in obsm_group.keys()
                    if 'polygon' in k.lower() or 'boundary' in k.lower()
                ]
                if polygon_keys:
                    storage_info["found_in"].append("obsm")
                    storage_info["details"]["obsm_keys"] = polygon_keys

            # Check uns for polygon data
            if 'uns' in root:
                uns_group = root['uns']
                polygon_keys = [
                    k for k in uns_group.keys()
                    if 'polygon' in k.lower() or 'boundary' in k.lower()
                ]
                if polygon_keys:
                    storage_info["found_in"].append("uns")
                    storage_info["details"]["uns_keys"] = polygon_keys

            # Check for raw polygon data in root
            polygon_arrays = [
                k for k in root.array_keys()
                if 'polygon' in k.lower() or 'boundary' in k.lower()
            ]
            if polygon_arrays:
                storage_info["found_in"].append("root_arrays")
                storage_info["details"]["array_keys"] = polygon_arrays

            return storage_info
        except Exception as e:
            self._log(f"Error inspecting polygon storage: {e}", "error")
            return None

    def load_polygons_from_uns(
        self, zarr_path: str
    ) -> Optional[Dict]:
        """Load polygons from uns in Zarr.

        Parameters
        ----------
        zarr_path : str
            Zarr directory path.

        Returns
        -------
        dict or None
            Polygon data from uns.
        """
        full_path = os.path.join(self.root_dir, zarr_path)
        if not os.path.isdir(full_path):
            self._log(f"Zarr directory not found", "error")
            return None

        try:
            adata = ad.read_zarr(full_path)
            if adata.uns:
                polygon_data = {}
                for key in adata.uns.keys():
                    if 'polygon' in key.lower() or 'boundary' in key.lower():
                        polygon_data[key] = adata.uns[key]
                        self._log(f"Found polygon data in uns: {key}")
                return polygon_data if polygon_data else None
            return None
        except Exception as e:
            self._log(f"Error loading polygons from uns: {e}", "error")
            return None


class ProsegTranscriptAssignmentsLoader:
    """Loader for Proseg transcript-to-cell assignments."""

    def __init__(self, root_dir: str, logger=None):
        """Initialize Proseg transcript assignments loader.

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

    def inspect_transcript_data_availability(
        self, zarr_path: str
    ) -> Optional[Dict]:
        """Inspect what transcript/molecule data is available.

        Parameters
        ----------
        zarr_path : str
            Zarr directory path.

        Returns
        -------
        dict or None
            Information about transcript data availability.
        """
        full_path = os.path.join(self.root_dir, zarr_path)
        if not os.path.isdir(full_path):
            return None

        try:
            adata = ad.read_zarr(full_path)
            info = {
                "obs_shape": adata.shape,
                "obs_is_molecules": False,  # Assume obs are cells by default
                "available_layers": list(adata.layers.keys()) if adata.layers else [],
                "available_obsm": list(adata.obsm.keys()) if adata.obsm else [],
                "transcript_assignment_columns": [],
                "molecule_id_columns": [],
            }

            # Look for transcript-related columns in obs
            obs_cols = adata.obs.columns
            for col in obs_cols:
                col_lower = col.lower()
                if any(t in col_lower for t in ["transcript", "molecule", "assignment"]):
                    info["transcript_assignment_columns"].append(col)
                if "id" in col_lower and "transcript" in col_lower:
                    info["molecule_id_columns"].append(col)

            return info
        except Exception as e:
            self._log(f"Error inspecting transcript data: {e}", "error")
            return None

    def get_molecule_ids(self, zarr_path: str) -> Optional[List[str]]:
        """Get unique molecule/transcript identifiers if available.

        Parameters
        ----------
        zarr_path : str
            Zarr directory path.

        Returns
        -------
        list or None
            Unique molecule IDs.
        """
        full_path = os.path.join(self.root_dir, zarr_path)
        if not os.path.isdir(full_path):
            return None

        try:
            adata = ad.read_zarr(full_path)

            # Look for molecule ID column
            for col in adata.obs.columns:
                if "id" in col.lower() and "transcript" in col.lower():
                    molecule_ids = adata.obs[col].unique().tolist()
                    self._log(f"Found {len(molecule_ids)} unique molecule IDs in '{col}'")
                    return molecule_ids

            self._log("No molecule ID column found", "warning")
            return None
        except Exception as e:
            self._log(f"Error getting molecule IDs: {e}", "error")
            return None

    def load_transcript_assignments(
        self, zarr_path: str, nrows: Optional[int] = None
    ) -> Optional[pd.DataFrame]:
        """Load transcript-to-cell assignments.

        Parameters
        ----------
        zarr_path : str
            Zarr directory path.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Transcript assignments.
        """
        full_path = os.path.join(self.root_dir, zarr_path)
        if not os.path.isdir(full_path):
            self._log(f"Zarr directory not found", "error")
            return None

        try:
            adata = ad.read_zarr(full_path)

            # Check if obs contains transcript/molecule data
            obs_data = adata.obs.copy()
            if nrows:
                obs_data = obs_data.iloc[:nrows]

            # Extract assignment-related columns
            assignment_cols = [
                col for col in obs_data.columns
                if any(a in col.lower() for a in ["cell", "assignment", "transcript", "molecule"])
            ]

            if assignment_cols:
                self._log(f"Found {len(assignment_cols)} assignment columns")
                return obs_data[assignment_cols]

            self._log("No assignment columns found", "warning")
            return None
        except Exception as e:
            self._log(f"Error loading assignments: {e}", "error")
            return None
