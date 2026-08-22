"""Proseg data loader for SegSure."""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import zarr


class ProsegLoader:
    """Loader for Proseg output data."""
    
    def __init__(self, root_dir: str, logger=None):
        """Initialize Proseg loader.
        
        Parameters
        ----------
        root_dir : str
            Root directory containing Proseg output files.
        logger : logging.Logger, optional
            Logger instance.
        """
        self.root_dir = root_dir
        self.logger = logger
    
    def _log(self, msg: str, level: str = "info"):
        """Helper for logging."""
        if self.logger:
            getattr(self.logger, level)(msg)
    
    def inspect_zarr_structure(self, zarr_path: str, max_depth: int = 5) -> Optional[Dict]:
        """Recursively inspect Zarr directory structure.
        
        Parameters
        ----------
        zarr_path : str
            Zarr directory path (relative or absolute).
        max_depth : int
            Maximum recursion depth.
        
        Returns
        -------
        dict or None
            Structure info including arrays, groups, and metadata.
        """
        if zarr_path.startswith('/'):
            full_path = zarr_path
        else:
            full_path = os.path.join(self.root_dir, zarr_path)
        
        if not os.path.isdir(full_path):
            self._log(f"Zarr directory not found: {full_path}", "warning")
            return None
        
        try:
            root = zarr.open_group(full_path, mode='r')
            structure = self._inspect_zarr_group(root, max_depth=max_depth)
            return structure
        except Exception as e:
            self._log(f"Error inspecting Zarr structure: {e}", "error")
            return None
    
    def _inspect_zarr_group(self, group, current_depth: int = 0, max_depth: int = 5) -> Dict:
        """Helper to recursively inspect Zarr groups."""
        if current_depth > max_depth:
            return {"truncated": True}
        
        info = {
            "type": "group",
            "path": str(group.path),
            "attrs": dict(group.attrs),
            "arrays": {},
            "groups": {},
        }
        
        # Inspect arrays
        for key in group.array_keys():
            try:
                arr = group[key]
                info["arrays"][key] = {
                    "shape": arr.shape,
                    "dtype": str(arr.dtype),
                    "chunks": arr.chunks,
                    "attrs": dict(arr.attrs),
                }
            except Exception as e:
                info["arrays"][key] = {"error": str(e)}
        
        # Inspect subgroups
        for key in group.group_keys():
            try:
                subgroup = group[key]
                info["groups"][key] = self._inspect_zarr_group(
                    subgroup, 
                    current_depth=current_depth + 1, 
                    max_depth=max_depth
                )
            except Exception as e:
                info["groups"][key] = {"error": str(e)}
        
        return info
    
    def get_zarr_stats(self, zarr_path: str) -> Optional[Dict]:
        """Get basic statistics about Zarr store without loading full data.
        
        Parameters
        ----------
        zarr_path : str
            Zarr directory path.
        
        Returns
        -------
        dict or None
            Statistics including size, array names, coordinate info.
        """
        if zarr_path.startswith('/'):
            full_path = zarr_path
        else:
            full_path = os.path.join(self.root_dir, zarr_path)
        
        if not os.path.isdir(full_path):
            return None
        
        try:
            # Check for adata.h5ad in Zarr
            h5ad_path = os.path.join(full_path, "adata.h5ad")
            size_mb = 0
            if os.path.exists(h5ad_path):
                size_mb += os.path.getsize(h5ad_path) / (1024 ** 2)
            
            # Size of entire directory
            for root, dirs, files in os.walk(full_path):
                for f in files:
                    size_mb += os.path.getsize(os.path.join(root, f)) / (1024 ** 2)
            
            root = zarr.open_group(full_path, mode='r')
            stats = {
                "path": full_path,
                "exists": True,
                "size_mb": size_mb,
                "root_attrs": dict(root.attrs),
                "arrays": list(root.array_keys()),
                "groups": list(root.group_keys()),
            }
            return stats
        except Exception as e:
            self._log(f"Error getting Zarr stats: {e}", "error")
            return None
    
    
    def load_h5ad(self, filename: str) -> Optional[ad.AnnData]:
        """Load Proseg H5AD output.
        
        Parameters
        ----------
        filename : str
            H5AD filename.
        
        Returns
        -------
        anndata.AnnData or None
            Proseg AnnData object.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            self._log(f"H5AD file not found: {filepath}", "error")
            return None
        
        self._log(f"Loading H5AD: {filepath}")
        
        adata = ad.read_h5ad(filepath)
        self._log(f"H5AD shape: {adata.shape}")
        self._log(f"H5AD obs columns: {list(adata.obs.columns)}")
        self._log(f"H5AD var columns (first 10): {list(adata.var.columns[:10])}")
        if adata.obsm:
            self._log(f"H5AD obsm: {list(adata.obsm.keys())}")
        if adata.uns:
            self._log(f"H5AD uns keys: {list(adata.uns.keys())[:5]}")
        
        return adata
    
    def load_zarr(self, filename: str) -> Optional[ad.AnnData]:
        """Load Proseg Zarr output.
        
        Parameters
        ----------
        filename : str
            Zarr directory name.
        
        Returns
        -------
        anndata.AnnData or None
            Proseg AnnData object.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isdir(filepath):
            self._log(f"Zarr directory not found: {filepath}", "error")
            return None
        
        self._log(f"Loading Zarr: {filepath}")
        
        adata = ad.read_zarr(filepath)
        self._log(f"Zarr shape: {adata.shape}")
        self._log(f"Zarr obs columns: {list(adata.obs.columns)}")
        self._log(f"Zarr var columns (first 10): {list(adata.var.columns[:10])}")
        if adata.obsm:
            self._log(f"Zarr obsm: {list(adata.obsm.keys())}")
        if adata.uns:
            self._log(f"Zarr uns keys: {list(adata.uns.keys())[:5]}")
        
        return adata
    
    def load_seurat(self, filename: str) -> Optional[object]:
        """Load Proseg Seurat output.
        
        Note: This requires R and rpy2 integration, which is not implemented here.
        Use H5AD or Zarr output for Python-side work.
        
        Parameters
        ----------
        filename : str
            Seurat RDS filename.
        
        Returns
        -------
        None
            Not implemented.
        """
        self._log(
            "Seurat RDS loading not implemented. Use H5AD or Zarr output.",
            "warning"
        )
        return None
