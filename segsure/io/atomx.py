"""AtoMx data loader for SegSure."""

import gzip
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
    
    def _log(self, msg: str, level: str = "info"):
        """Helper for logging."""
        if self.logger:
            getattr(self.logger, level)(msg)
    
    def _open_file(self, filepath: str):
        """Open a file, handling gzip compression if needed."""
        if filepath.endswith('.gz'):
            return gzip.open(filepath, 'rt')
        return open(filepath, 'r')
    
    def inspect_file_schema(self, filename: str) -> Optional[Dict]:
        """Inspect file schema without loading full data.
        
        Parameters
        ----------
        filename : str
            File to inspect.
        
        Returns
        -------
        dict or None
            Schema info: {'path', 'exists', 'size_mb', 'columns', 'n_rows_estimate', 'dtypes'}
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            self._log(f"File not found: {filepath}", "warning")
            return None
        
        schema = {
            'path': filepath,
            'exists': True,
            'size_mb': os.path.getsize(filepath) / (1024 ** 2),
            'filename': filename,
        }
        
        try:
            # Read just the header and first few rows to infer schema
            with self._open_file(filepath) as f:
                header = f.readline().strip().split(',')
                schema['columns'] = header
                schema['n_columns'] = len(header)
                
                # Count rows by reading the file
                f.seek(0)
                n_rows = sum(1 for _ in f) - 1  # -1 for header
                schema['n_rows_estimate'] = n_rows
        except Exception as e:
            self._log(f"Error inspecting {filename}: {e}", "error")
            return None
        
        return schema
    
    def load_expression_matrix(
        self, 
        filename: str, 
        nrows: Optional[int] = None,
        fov: Optional[str] = None,
    ) -> Optional[pd.DataFrame]:
        """Load gene expression matrix.
        
        Parameters
        ----------
        filename : str
            Expression matrix filename.
        nrows : int, optional
            Maximum number of rows to load.
        fov : str, optional
            FOV to filter by (requires 'fov' column in metadata).
        
        Returns
        -------
        pd.DataFrame or None
            Expression matrix (cells x genes).
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            self._log(f"Expression file not found: {filepath}", "error")
            return None
        
        self._log(f"Loading expression matrix: {filepath}")
        
        df = pd.read_csv(filepath, index_col=0, nrows=nrows)
        self._log(f"Expression matrix shape: {df.shape}")
        
        return df
    
    def load_fov_positions(
        self, 
        filename: str,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load FOV positions.
        
        Parameters
        ----------
        filename : str
            FOV positions filename.
        nrows : int, optional
            Maximum number of rows to load.
        
        Returns
        -------
        pd.DataFrame or None
            FOV position data.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            self._log(f"FOV positions file not found: {filepath}", "error")
            return None
        
        self._log(f"Loading FOV positions: {filepath}")
        
        df = pd.read_csv(filepath, nrows=nrows)
        self._log(f"FOV positions shape: {df.shape}")
        
        return df
    
    def load_metadata(
        self, 
        filename: str,
        nrows: Optional[int] = None,
        fov: Optional[str] = None,
    ) -> Optional[pd.DataFrame]:
        """Load cell metadata.
        
        Parameters
        ----------
        filename : str
            Metadata filename.
        nrows : int, optional
            Maximum number of rows to load.
        fov : str, optional
            FOV to filter by (requires 'fov' column in metadata).
        
        Returns
        -------
        pd.DataFrame or None
            Cell metadata.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            self._log(f"Metadata file not found: {filepath}", "error")
            return None
        
        self._log(f"Loading metadata: {filepath}")
        
        df = pd.read_csv(filepath, index_col=0, nrows=nrows)
        
        if fov:
            if 'fov' in df.columns:
                df = df[df['fov'] == fov]
                self._log(f"Filtered to FOV {fov}: {df.shape[0]} cells")
            else:
                self._log(f"Warning: 'fov' column not found in metadata", "warning")
        
        self._log(f"Metadata shape: {df.shape}")
        
        return df
    
    def load_transcripts(
        self, 
        filename: str,
        nrows: Optional[int] = None,
        fov: Optional[str] = None,
    ) -> Optional[pd.DataFrame]:
        """Load transcript data.
        
        Parameters
        ----------
        filename : str
            Transcript filename.
        nrows : int, optional
            Maximum number of rows to load.
        fov : str, optional
            FOV to filter by (requires 'fov' column in transcripts).
        
        Returns
        -------
        pd.DataFrame or None
            Transcript data.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            self._log(f"Transcript file not found: {filepath}", "error")
            return None
        
        self._log(f"Loading transcripts: {filepath}")
        
        df = pd.read_csv(filepath, nrows=nrows)
        
        if fov:
            if 'fov' in df.columns:
                df = df[df['fov'] == fov]
                self._log(f"Filtered to FOV {fov}: {df.shape[0]} transcripts")
            else:
                self._log(f"Warning: 'fov' column not found in transcripts", "warning")
        
        self._log(f"Transcript data shape: {df.shape}")
        
        return df
    
    def load_polygons(
        self, 
        filename: str,
        nrows: Optional[int] = None,
        fov: Optional[str] = None,
    ) -> Optional[pd.DataFrame]:
        """Load polygon/segmentation data.
        
        Parameters
        ----------
        filename : str
            Polygons filename.
        nrows : int, optional
            Maximum number of rows to load.
        fov : str, optional
            FOV to filter by (requires 'fov' column in polygons).
        
        Returns
        -------
        pd.DataFrame or None
            Polygon vertex data.
        """
        filepath = os.path.join(self.root_dir, filename)
        if not os.path.isfile(filepath):
            self._log(f"Polygons file not found: {filepath}", "error")
            return None
        
        self._log(f"Loading polygons: {filepath}")
        
        df = pd.read_csv(filepath, nrows=nrows)
        
        if fov:
            if 'fov' in df.columns:
                df = df[df['fov'] == fov]
                self._log(f"Filtered to FOV {fov}: {df.shape[0]} polygon vertices")
            else:
                self._log(f"Warning: 'fov' column not found in polygons", "warning")
        
        self._log(f"Polygon data shape: {df.shape}")
        
        return df
    
    def load_all(
        self,
        expression: str,
        fov_positions: str,
        metadata: str,
        transcripts: str,
        polygons: str,
        nrows: Optional[int] = None,
        fov: Optional[str] = None,
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
        nrows : int, optional
            Maximum number of rows to load.
        fov : str, optional
            FOV to filter by.
        
        Returns
        -------
        dict
            Dictionary with all loaded dataframes.
        """
        data = {
            "expression": self.load_expression_matrix(expression, nrows=nrows, fov=fov),
            "fov_positions": self.load_fov_positions(fov_positions, nrows=nrows),
            "metadata": self.load_metadata(metadata, nrows=nrows, fov=fov),
            "transcripts": self.load_transcripts(transcripts, nrows=nrows, fov=fov),
            "polygons": self.load_polygons(polygons, nrows=nrows, fov=fov),
        }
        return data
    
    def get_available_fovs(self, fov_positions_file: str) -> Optional[List[str]]:
        """Get list of available FOVs from FOV positions file.
        
        Parameters
        ----------
        fov_positions_file : str
            FOV positions filename.
        
        Returns
        -------
        list or None
            List of unique FOV identifiers.
        """
        fov_df = self.load_fov_positions(fov_positions_file)
        if fov_df is None:
            return None
        
        # Try common column names for FOV identifier
        fov_cols = ['fov', 'fov_id', 'fov_name', 'slide_id']
        for col in fov_cols:
            if col in fov_df.columns:
                return sorted(fov_df[col].unique().tolist())
        
        self._log("Could not find FOV column in FOV positions", "warning")
        return None
    
    def inspect_transcript_columns(self, transcripts_file: str) -> Optional[Dict]:
        """Inspect transcript data columns without loading all data.
        
        Parameters
        ----------
        transcripts_file : str
            Transcript filename.
        
        Returns
        -------
        dict or None
            Schema info including columns, dtypes, and coordinate columns.
        """
        schema = self.inspect_file_schema(transcripts_file)
        if schema is None:
            return None
        
        # Try to identify coordinate columns
        coord_info = {}
        coord_cols = ['x', 'y', 'z', 'x_coord', 'y_coord', 'z_coord', 'cx', 'cy']
        for col in schema.get('columns', []):
            if any(coord in col.lower() for coord in coord_cols):
                coord_info[col.lower()] = col
        
        schema['coordinate_columns'] = coord_info
        
        # Try to identify gene column
        gene_cols = ['gene', 'target', 'gene_name']
        for col in schema.get('columns', []):
            if any(g in col.lower() for g in gene_cols):
                schema['gene_column'] = col
                break
        
        # Try to identify cell assignment column
        cell_cols = ['cell_id', 'cellid', 'cell', 'nucleus_id', 'segment_id']
        for col in schema.get('columns', []):
            if any(c in col.lower() for c in cell_cols):
                schema['cell_id_column'] = col
                break
        
        return schema
