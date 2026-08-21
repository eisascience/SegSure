"""AtoMx data loader for SegSure."""

import os
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd


class AtoMxLoader:
    """Loader for AtoMx-exported CosMx data."""
    
    def __init__(self, root_dir: str, logger=None):
        """Initialize AtoMx loader.
        
        Parameters
        ----------
        root_dir : str
            Root directory containing AtoMx export files.
        logger : logging.Logger, optional
            Logger instance.
        """
        self.root_dir = root_dir
        self.logger = logger
    
    def load_expression_matrix(self, filename: str) -> Optional[pd.DataFrame]:
        """Load gene expression matrix.
        
        Parameters
        ----------
        filename : str
            Expression matrix filename.
        
        Returns
        -------
        pd.DataFrame or None
            Expression matrix (cells x genes).
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            if self.logger:
                self.logger.error(f"Expression file not found: {filepath}")
            return None
        
        if self.logger:
            self.logger.info(f"Loading expression matrix: {filepath}")
        
        df = pd.read_csv(filepath, index_col=0)
        if self.logger:
            self.logger.info(f"Expression matrix shape: {df.shape}")
        
        return df
    
    def load_fov_positions(self, filename: str) -> Optional[pd.DataFrame]:
        """Load FOV positions.
        
        Parameters
        ----------
        filename : str
            FOV positions filename.
        
        Returns
        -------
        pd.DataFrame or None
            FOV position data.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            if self.logger:
                self.logger.error(f"FOV positions file not found: {filepath}")
            return None
        
        if self.logger:
            self.logger.info(f"Loading FOV positions: {filepath}")
        
        df = pd.read_csv(filepath)
        if self.logger:
            self.logger.info(f"FOV positions shape: {df.shape}")
        
        return df
    
    def load_metadata(self, filename: str) -> Optional[pd.DataFrame]:
        """Load cell metadata.
        
        Parameters
        ----------
        filename : str
            Metadata filename.
        
        Returns
        -------
        pd.DataFrame or None
            Cell metadata.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            if self.logger:
                self.logger.error(f"Metadata file not found: {filepath}")
            return None
        
        if self.logger:
            self.logger.info(f"Loading metadata: {filepath}")
        
        df = pd.read_csv(filepath, index_col=0)
        if self.logger:
            self.logger.info(f"Metadata shape: {df.shape}")
        
        return df
    
    def load_transcripts(self, filename: str) -> Optional[pd.DataFrame]:
        """Load transcript data.
        
        Parameters
        ----------
        filename : str
            Transcript filename.
        
        Returns
        -------
        pd.DataFrame or None
            Transcript data.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            if self.logger:
                self.logger.error(f"Transcript file not found: {filepath}")
            return None
        
        if self.logger:
            self.logger.info(f"Loading transcripts: {filepath}")
        
        df = pd.read_csv(filepath)
        if self.logger:
            self.logger.info(f"Transcript data shape: {df.shape}")
        
        return df
    
    def load_polygons(self, filename: str) -> Optional[pd.DataFrame]:
        """Load polygon/segmentation data.
        
        Parameters
        ----------
        filename : str
            Polygons filename.
        
        Returns
        -------
        pd.DataFrame or None
            Polygon vertex data.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            if self.logger:
                self.logger.error(f"Polygons file not found: {filepath}")
            return None
        
        if self.logger:
            self.logger.info(f"Loading polygons: {filepath}")
        
        df = pd.read_csv(filepath)
        if self.logger:
            self.logger.info(f"Polygon data shape: {df.shape}")
        
        return df
    
    def load_all(
        self,
        expression: str,
        fov_positions: str,
        metadata: str,
        transcripts: str,
        polygons: str,
    ) -> Dict[str, Optional[pd.DataFrame]]:
        """Load all AtoMx files.
        
        Parameters
        ----------
        expression : str
            Expression matrix filename.
        fov_positions : str
            FOV positions filename.
        metadata : str
            Metadata filename.
        transcripts : str
            Transcript filename.
        polygons : str
            Polygons filename.
        
        Returns
        -------
        dict
            Dictionary with all loaded dataframes.
        """
        data = {
            "expression": self.load_expression_matrix(expression),
            "fov_positions": self.load_fov_positions(fov_positions),
            "metadata": self.load_metadata(metadata),
            "transcripts": self.load_transcripts(transcripts),
            "polygons": self.load_polygons(polygons),
        }
        return data
