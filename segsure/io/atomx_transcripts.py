"""Specialized loader for AtoMx transcript-level data."""

import os
from typing import Dict, List, Optional, Tuple

import pandas as pd


class AtoMxTranscriptsLoader:
    """Loader for AtoMx transcript-level data."""

    def __init__(self, root_dir: str, logger=None):
        """Initialize AtoMx transcripts loader.

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

    def load_transcripts(
        self,
        transcripts_file: str,
        fov: Optional[str] = None,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load transcript data.

        Parameters
        ----------
        transcripts_file : str
            Transcript filename.
        fov : str, optional
            FOV to filter by.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Transcript data.
        """
        filepath = os.path.join(self.root_dir, transcripts_file)
        if not os.path.isfile(filepath):
            self._log(f"Transcript file not found: {filepath}", "error")
            return None

        self._log(f"Loading transcripts from: {transcripts_file}")
        transcripts = pd.read_csv(filepath, nrows=nrows)

        if fov and "fov" in transcripts.columns:
            transcripts = transcripts[transcripts["fov"] == fov]
            self._log(f"Filtered to FOV {fov}: {len(transcripts)} transcripts")

        self._log(f"Loaded {len(transcripts)} transcripts")
        return transcripts

    def get_transcript_columns(
        self, transcripts_file: str
    ) -> Optional[Dict[str, List[str]]]:
        """Get categorized column names for transcripts.

        Parameters
        ----------
        transcripts_file : str
            Transcript filename.

        Returns
        -------
        dict or None
            Categorized columns: {'coordinates', 'genes', 'assignments', 'metadata'}.
        """
        filepath = os.path.join(self.root_dir, transcripts_file)
        if not os.path.isfile(filepath):
            return None

        # Read just the header
        with open(filepath, "r") as f:
            header = f.readline().strip().split(",")

        columns = {
            "coordinates": [],
            "genes": [],
            "assignments": [],
            "metadata": [],
            "fov": [],
        }

        for col in header:
            col_lower = col.lower()
            if any(c in col_lower for c in ["x", "y", "z", "coord"]):
                columns["coordinates"].append(col)
            elif any(g in col_lower for g in ["gene", "target"]):
                columns["genes"].append(col)
            elif any(a in col_lower for a in ["cell", "nucleus", "segment"]):
                columns["assignments"].append(col)
            elif "fov" in col_lower:
                columns["fov"].append(col)
            else:
                columns["metadata"].append(col)

        self._log(f"Identified column categories: {columns}")
        return columns

    def load_transcript_coordinates(
        self,
        transcripts_file: str,
        fov: Optional[str] = None,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load transcript coordinates only.

        Parameters
        ----------
        transcripts_file : str
            Transcript filename.
        fov : str, optional
            FOV to filter by.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Transcript coordinates.
        """
        transcripts = self.load_transcripts(
            transcripts_file, fov=fov, nrows=nrows
        )
        if transcripts is None:
            return None

        # Find coordinate columns
        coord_cols = [
            col for col in transcripts.columns
            if any(c in col.lower() for c in ["x", "y", "z", "coord"])
        ]

        if not coord_cols:
            self._log("No coordinate columns found", "warning")
            return None

        self._log(f"Found {len(coord_cols)} coordinate columns")
        return transcripts[coord_cols]

    def load_transcript_genes(
        self,
        transcripts_file: str,
        fov: Optional[str] = None,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load transcript gene identities.

        Parameters
        ----------
        transcripts_file : str
            Transcript filename.
        fov : str, optional
            FOV to filter by.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Transcript gene data.
        """
        transcripts = self.load_transcripts(
            transcripts_file, fov=fov, nrows=nrows
        )
        if transcripts is None:
            return None

        # Find gene columns
        gene_cols = [
            col for col in transcripts.columns
            if any(g in col.lower() for g in ["gene", "target"])
        ]

        if not gene_cols:
            self._log("No gene columns found", "warning")
            return None

        self._log(f"Found {len(gene_cols)} gene columns")
        return transcripts[gene_cols]

    def load_transcript_assignments(
        self,
        transcripts_file: str,
        fov: Optional[str] = None,
        nrows: Optional[int] = None,
    ) -> Optional[pd.DataFrame]:
        """Load transcript-to-cell assignments.

        Parameters
        ----------
        transcripts_file : str
            Transcript filename.
        fov : str, optional
            FOV to filter by.
        nrows : int, optional
            Maximum number of rows to load.

        Returns
        -------
        pd.DataFrame or None
            Transcript-to-cell assignments.
        """
        transcripts = self.load_transcripts(
            transcripts_file, fov=fov, nrows=nrows
        )
        if transcripts is None:
            return None

        # Find assignment columns
        assign_cols = [
            col for col in transcripts.columns
            if any(a in col.lower() for a in ["cell", "nucleus", "segment"])
        ]

        if not assign_cols:
            self._log("No assignment columns found", "warning")
            return None

        self._log(f"Found {len(assign_cols)} assignment columns")
        return transcripts[assign_cols]
