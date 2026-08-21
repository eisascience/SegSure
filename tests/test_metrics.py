"""Tests for metrics computation."""

import pandas as pd
import pytest

from segsure.metrics import transcript_assignment


class TestMetrics:
    """Test metrics computation functions."""
    
    @pytest.fixture
    def sample_comparisons(self):
        """Create sample transcript comparison data."""
        return pd.DataFrame({
            "atomx_cell": ["c1", "c2", "c3"],
            "proseg_cell": ["p1", "p2", "p3"],
            "atomx_tx_count": [10, 15, 20],
            "proseg_tx_count": [9, 15, 20],
            "overlap_tx_count": [8, 14, 19],
            "jaccard_index": [0.8, 0.93, 0.95],
        })
    
    def test_compute_overlap_metrics(self, sample_comparisons):
        """Test transcript overlap metric computation."""
        result = transcript_assignment.compute_transcript_overlap_metrics(
            assignment_comparisons=sample_comparisons
        )
        
        assert len(result) == 3
        assert "discord_rate" in result.columns
        assert "atomx_recall" in result.columns
        assert "proseg_recall" in result.columns
    
    def test_compute_discord_severity(self, sample_comparisons):
        """Test discord severity classification."""
        result = transcript_assignment.compute_discord_severity(
            assignment_comparisons=sample_comparisons
        )
        
        assert "discord_severity" in result.columns
        assert all(s in ["none", "low", "medium", "high"] for s in result["discord_severity"])
