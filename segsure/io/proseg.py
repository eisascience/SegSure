"""Proseg data loader for SegSure."""

import os
from typing import Dict, Optional

import anndata as ad


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
            if self.logger:
                self.logger.error(f"H5AD file not found: {filepath}")
            return None
        
        if self.logger:
            self.logger.info(f"Loading H5AD: {filepath}")
        
        adata = ad.read_h5ad(filepath)
        if self.logger:
            self.logger.info(f"H5AD shape: {adata.shape}")
            self.logger.info(f"H5AD obs columns: {list(adata.obs.columns)}")
            self.logger.info(f"H5AD var columns: {list(adata.var.columns)}")
            if adata.obsm:
                self.logger.info(f"H5AD obsm: {list(adata.obsm.keys())}")
            if adata.uns:
                self.logger.info(f"H5AD uns keys: {list(adata.uns.keys())}")
        
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
            if self.logger:
                self.logger.error(f"Zarr directory not found: {filepath}")
            return None
        
        if self.logger:
            self.logger.info(f"Loading Zarr: {filepath}")
        
        adata = ad.read_zarr(filepath)
        if self.logger:
            self.logger.info(f"Zarr shape: {adata.shape}")
            self.logger.info(f"Zarr obs columns: {list(adata.obs.columns)}")
            self.logger.info(f"Zarr var columns: {list(adata.var.columns)}")
            if adata.obsm:
                self.logger.info(f"Zarr obsm: {list(adata.obsm.keys())}")
            if adata.uns:
                self.logger.info(f"Zarr uns keys: {list(adata.uns.keys())}")
        
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
        if self.logger:
            self.logger.warning(
                "Seurat RDS loading not implemented. Use H5AD or Zarr output."
            )
        return None
