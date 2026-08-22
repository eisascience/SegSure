"""Tests for SegSure I/O modules."""

import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from segsure.io import atomx


class TestAtoMxLoader:
    """Test AtoMx loader."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory with test data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create dummy CSV files
            test_data = {
                "expr": pd.DataFrame(
                    {"gene1": [1, 2], "gene2": [3, 4]},
                    index=["cell1", "cell2"]
                ),
                "meta": pd.DataFrame(
                    {"x": [10, 20], "y": [30, 40]},
                    index=["cell1", "cell2"]
                ),
                "tx": pd.DataFrame({
                    "transcript_id": ["tx1", "tx2", "tx3"],
                    "gene": ["gene1", "gene2", "gene1"],
                    "x": [15, 25, 35],
                    "y": [35, 45, 55],
                    "cell_id": ["cell1", "cell2", "cell1"],
                }),
            }
            
            test_data["expr"].to_csv(os.path.join(tmpdir, "expr.csv.gz"))
            test_data["meta"].to_csv(os.path.join(tmpdir, "meta.csv.gz"))
            test_data["tx"].to_csv(os.path.join(tmpdir, "tx.csv.gz"), index=False)
            
            yield tmpdir
    
    def test_load_expression_matrix(self, temp_dir):
        """Test loading expression matrix."""
        loader = atomx.AtoMxLoader(temp_dir)
        df = loader.load_expression_matrix("expr.csv.gz")
        
        assert df is not None
        assert df.shape == (2, 2)
        assert "gene1" in df.columns
    
    def test_load_metadata(self, temp_dir):
        """Test loading metadata."""
        loader = atomx.AtoMxLoader(temp_dir)
        df = loader.load_metadata("meta.csv.gz")
        
        assert df is not None
        assert df.shape == (2, 2)
        assert "x" in df.columns
        assert "y" in df.columns
    
    def test_load_transcripts(self, temp_dir):
        """Test loading transcripts."""
        loader = atomx.AtoMxLoader(temp_dir)
        df = loader.load_transcripts("tx.csv.gz")
        
        assert df is not None
        assert df.shape[0] == 3
        assert "transcript_id" in df.columns
    
    def test_file_not_found(self, temp_dir):
        """Test loading non-existent file."""
        loader = atomx.AtoMxLoader(temp_dir)
        df = loader.load_expression_matrix("nonexistent.csv.gz")
        
        assert df is None
