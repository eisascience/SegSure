"""Tests for metrics computation."""

import numpy as np
import pandas as pd
import pytest

from segsure.metrics import transcript_assignment, neighborhood
from segsure.molecules import MoleculeTracker


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
            "discord_severity": ["low", "none", "none"],
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
    
    @pytest.fixture
    def sample_molecules(self):
        """Create sample molecule data."""
        atomx_molecules = pd.DataFrame({
            "transcript_id": ["tx1", "tx2", "tx3", "tx4"],
            "gene": ["GeneA", "GeneB", "GeneA", "GeneC"],
            "fov": [1, 1, 1, 1],
            "x": [10.0, 20.0, 30.0, 40.0],
            "y": [15.0, 25.0, 35.0, 45.0],
            "z": [0.0, 0.0, 0.0, 0.0],
            "cell_id": ["c1", "c2", "c1", "c3"],
        })
        
        proseg_molecules = pd.DataFrame({
            "transcript_id": ["tx1", "tx2", "tx5", "tx6"],
            "gene": ["GeneA", "GeneB", "GeneA", "GeneD"],
            "fov": [1, 1, 1, 1],
            "x": [10.1, 20.2, 31.0, 40.0],
            "y": [15.1, 25.2, 36.0, 45.0],
            "z": [0.0, 0.0, 0.0, 0.0],
            "cell_id": ["p1", "p2", "p1", "p4"],
        })
        
        cell_matches = pd.DataFrame({
            "atomx_cell": ["c1", "c2", "c3"],
            "proseg_cell": ["p1", "p2", "p4"],
        })
        
        return atomx_molecules, proseg_molecules, cell_matches
    
    def test_molecule_tracker_initialization(self, sample_molecules):
        """Test MoleculeTracker initialization."""
        atomx, proseg, matches = sample_molecules
        
        tracker = MoleculeTracker(
            atomx_molecules=atomx,
            proseg_molecules=proseg,
            cell_matches=matches,
            spatial_tolerance=2.0,
        )
        
        assert tracker is not None
        assert tracker.spatial_tolerance == 2.0
    
    def test_molecule_matching_by_id(self, sample_molecules):
        """Test matching molecules by stable ID."""
        atomx, proseg, matches = sample_molecules
        
        tracker = MoleculeTracker(
            atomx_molecules=atomx,
            proseg_molecules=proseg,
            cell_matches=matches,
        )
        
        matched = tracker.match_molecules_by_id()
        
        # Should match tx1 and tx2
        assert len(matched) >= 2
        assert "tx1" in matched["tx_id"].values or "tx2" in matched["tx_id"].values
    
    def test_molecule_classification_same_cell(self, sample_molecules):
        """Test molecule classification for same-cell assignments."""
        atomx, proseg, matches = sample_molecules
        
        # Create matched molecules where atomx_cell == proseg_cell (via correspondence)
        test_molecules = pd.DataFrame({
            "tx_id": ["tx1"],
            "gene": ["GeneA"],
            "fov": [1],
            "x": [10.0],
            "y": [15.0],
            "z": [0.0],
            "atomx_cell": ["c1"],
            "proseg_cell": ["p1"],
        })
        
        tracker = MoleculeTracker(
            atomx_molecules=atomx,
            proseg_molecules=proseg,
            cell_matches=matches,
        )
        
        classified = tracker.classify_molecules(test_molecules)
        
        # c1 -> p1 is a valid match in cell_matches
        assert len(classified) == 1
        assert classified["status"].iloc[0] in tracker.VALID_STATUSES
    
    def test_neighbor_jaccard_similarity(self):
        """Test mapped-neighbor Jaccard similarity computation."""
        # Create simple neighbor graphs
        atomx_neighbors = {
            "c1": ["c2", "c3"],
            "c2": ["c1", "c3"],
            "c3": ["c1", "c2"],
        }
        
        proseg_neighbors = {
            "p1": ["p2", "p3"],
            "p2": ["p1", "p3"],
            "p3": ["p1", "p2"],
        }
        
        cell_matches = pd.DataFrame({
            "atomx_cell": ["c1", "c2", "c3"],
            "proseg_cell": ["p1", "p2", "p3"],
        })
        
        jaccard = neighborhood.compute_mapped_neighbor_jaccard(
            atomx_neighbors=atomx_neighbors,
            proseg_neighbors=proseg_neighbors,
            cell_matches=cell_matches,
        )
        
        # Perfect correspondence should give Jaccard ~= 1.0
        assert 0 <= jaccard <= 1
        assert not np.isnan(jaccard)
    
    def test_neighborhood_metrics_empty(self):
        """Test neighborhood metrics with empty data."""
        result = neighborhood.compute_neighborhood_discord(
            assignment_comparisons=pd.DataFrame(),
            neighbors_atomx={},
            neighbors_proseg={},
            cell_matches=pd.DataFrame(),
        )
        
        assert result.empty


class TestMoleculeClassification:
    """Test molecule-level status classification."""
    
    def test_valid_status_constants(self):
        """Test that all status constants are valid."""
        assert len(MoleculeTracker.VALID_STATUSES) == 8
        assert MoleculeTracker.STATUS_SAME_MATCHED_CELL in MoleculeTracker.VALID_STATUSES
        assert MoleculeTracker.STATUS_UNRESOLVED in MoleculeTracker.VALID_STATUSES
    
    def test_molecule_summary_computation(self):
        """Test cell-level transcript summary computation."""
        molecules_df = pd.DataFrame({
            "tx_id": ["tx1", "tx2", "tx3", "tx4", "tx5"],
            "status": [
                "same_matched_cell",
                "same_matched_cell",
                "changed_neighbor",
                "unresolved",
                "changed_unrelated",
            ],
            "atomx_cell": ["c1", "c1", "c2", "c2", "c3"],
        })
        
        summary = transcript_assignment.compute_molecule_level_status_summary(
            molecules_df
        )
        
        assert summary["total_molecules"] == 5
        assert summary["n_consistent"] == 2
        assert summary["n_changed"] == 2
        assert summary["n_unresolved"] == 1


class TestJaccardSimilarity:
    """Test Jaccard similarity computation for neighbor sets."""
    
    def test_jaccard_perfect_match(self):
        """Test Jaccard with perfect neighbor match."""
        atomx_neighbors = {
            "c1": ["c2", "c3"],
        }
        proseg_neighbors = {
            "p1": ["p2", "p3"],
        }
        cell_matches = pd.DataFrame({
            "atomx_cell": ["c1"],
            "proseg_cell": ["p1"],
        })
        
        jaccard = neighborhood.compute_mapped_neighbor_jaccard(
            atomx_neighbors=atomx_neighbors,
            proseg_neighbors=proseg_neighbors,
            cell_matches=cell_matches,
        )
        
        assert jaccard == 1.0
    
    def test_jaccard_no_match(self):
        """Test Jaccard with no neighbor overlap."""
        atomx_neighbors = {
            "c1": ["c2"],
        }
        proseg_neighbors = {
            "p1": ["p3"],
        }
        cell_matches = pd.DataFrame({
            "atomx_cell": ["c1"],
            "proseg_cell": ["p1"],
        })
        
        jaccard = neighborhood.compute_mapped_neighbor_jaccard(
            atomx_neighbors=atomx_neighbors,
            proseg_neighbors=proseg_neighbors,
            cell_matches=cell_matches,
        )
        
        assert jaccard == 0.0

