"""Tests for cell matching modules."""

import pandas as pd
import pytest

from segsure.harmonize import cells


class TestCellMatching:
    """Test cell matching functions."""
    
    @pytest.fixture
    def sample_data(self):
        """Create sample cell data for testing."""
        atomx_cells = pd.DataFrame({
            "cell_id": ["a1", "a2", "a3"],
            "x": [10.0, 20.0, 30.0],
            "y": [10.0, 20.0, 30.0],
        })
        
        proseg_cells = pd.DataFrame({
            "cell_id": ["p1", "p2", "p3"],
            "x": [10.5, 20.5, 30.5],
            "y": [10.5, 20.5, 30.5],
        })
        
        return {"atomx": atomx_cells, "proseg": proseg_cells}
    
    def test_match_cells(self, sample_data):
        """Test basic cell matching."""
        matches = cells.match_cells(
            atomx_cells=sample_data["atomx"],
            proseg_cells=sample_data["proseg"],
            atomx_coords={"x": "x", "y": "y"},
            proseg_coords={"x": "x", "y": "y"},
            atomx_cell_id="cell_id",
            proseg_cell_id="cell_id",
            distance_threshold=1.0,
        )
        
        assert len(matches) == 3
        assert "atomx_cell" in matches.columns
        assert "proseg_cell" in matches.columns
        assert "distance" in matches.columns
    
    def test_identify_split_merges(self, sample_data):
        """Test split/merge identification."""
        matches = pd.DataFrame({
            "atomx_cell": ["a1", "a2", "a2"],
            "proseg_cell": ["p1", "p2", "p3"],
        })
        
        result = cells.identify_split_merges(matches)
        
        assert "splits" in result
        assert "merges" in result
        assert len(result["splits"]) == 1  # a2 is split
