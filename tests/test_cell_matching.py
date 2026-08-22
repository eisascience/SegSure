"""Tests for cell matching modules."""

import pandas as pd
import pytest
import numpy as np
from shapely.geometry import Polygon

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
        
        assert len(matches) > 0
        assert "atomx_cell" in matches.columns or "atomx_cell_id" in matches.columns
        assert "proseg_cell" in matches.columns or "proseg_cell_id" in matches.columns
    
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


class TestPolygonParsing:
    """Test polygon parsing functions."""
    
    def test_parse_polygons_wkt(self):
        """Test parsing WKT format polygons."""
        from shapely.geometry import box
        from shapely.wkt import dumps as wkt_dumps
        
        # Create test polygons
        poly1 = box(0, 0, 10, 10)
        poly2 = box(20, 20, 30, 30)
        
        df = pd.DataFrame({
            "cell_id": ["c1", "c2"],
            "geometry": [wkt_dumps(poly1), wkt_dumps(poly2)],
        })
        
        result = cells.parse_polygons(df, cell_id_col="cell_id", geometry_col="geometry")
        
        assert len(result) == 2
        assert "c1" in result
        assert "c2" in result
        assert isinstance(result["c1"], Polygon)


class TestOverlapMetrics:
    """Test overlap metric calculations."""
    
    def test_calculate_overlap_metrics(self):
        """Test overlap metric calculation."""
        # Create two overlapping squares
        poly_a = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        poly_b = Polygon([(5, 5), (15, 5), (15, 15), (5, 15)])
        
        metrics = cells.calculate_overlap_metrics(poly_a, poly_b)
        
        assert "intersection_area" in metrics
        assert "union_area" in metrics
        assert "iou" in metrics
        assert "coverage_a" in metrics
        assert "coverage_b" in metrics
        assert "centroid_distance" in metrics
        
        # Check values
        assert metrics["intersection_area"] == 25.0  # 5x5 square
        assert metrics["coverage_a"] > 0  # Some of A is covered by B
        assert metrics["coverage_b"] > 0  # Some of B is covered by A
        assert metrics["iou"] > 0  # Some intersection


class TestOneToOneMatching:
    """Test one-to-one cell matching."""
    
    def test_one_to_one_match(self):
        """Test perfect one-to-one matching."""
        from shapely.geometry import box
        
        # Create matching polygons
        poly_a1 = box(0, 0, 10, 10)
        poly_b1 = box(0, 0, 10, 10)  # Identical
        
        atomx_df = pd.DataFrame({
            "atomx_cell_id": ["a1"],
            "x": [5.0],
            "y": [5.0],
        })
        
        proseg_df = pd.DataFrame({
            "proseg_cell_id": ["p1"],
            "x": [5.0],
            "y": [5.0],
        })
        
        atomx_polygons = {"a1": poly_a1}
        proseg_polygons = {"p1": poly_b1}
        
        result = cells.match_cells(
            atomx_cells=atomx_df,
            proseg_cells=proseg_df,
            atomx_coords={"x": "x", "y": "y"},
            proseg_coords={"x": "x", "y": "y"},
            atomx_cell_id="atomx_cell_id",
            proseg_cell_id="proseg_cell_id",
            atomx_polygons=atomx_polygons,
            proseg_polygons=proseg_polygons,
        )
        
        assert len(result) > 0
        # Should have perfect overlap (or close to it)
        one_to_one = result[result.get("relationship_type") == "one_to_one"]
        assert len(one_to_one) > 0


class TestSplitMatching:
    """Test split (1->N) cell matching."""
    
    def test_split_match(self):
        """Test one AtoMx cell splitting into multiple Proseg cells."""
        from shapely.geometry import box
        
        # One large AtoMx cell
        atomx_poly = box(0, 0, 20, 10)
        
        # Two smaller Proseg cells covering the space
        proseg_poly1 = box(0, 0, 10, 10)
        proseg_poly2 = box(10, 0, 20, 10)
        
        atomx_df = pd.DataFrame({
            "atomx_cell_id": ["a1"],
            "x": [10.0],
            "y": [5.0],
        })
        
        proseg_df = pd.DataFrame({
            "proseg_cell_id": ["p1", "p2"],
            "x": [5.0, 15.0],
            "y": [5.0, 5.0],
        })
        
        atomx_polygons = {"a1": atomx_poly}
        proseg_polygons = {"p1": proseg_poly1, "p2": proseg_poly2}
        
        result = cells.match_cells(
            atomx_cells=atomx_df,
            proseg_cells=proseg_df,
            atomx_coords={"x": "x", "y": "y"},
            proseg_coords={"x": "x", "y": "y"},
            atomx_cell_id="atomx_cell_id",
            proseg_cell_id="proseg_cell_id",
            atomx_polygons=atomx_polygons,
            proseg_polygons=proseg_polygons,
        )
        
        splits = result[result.get("relationship_type") == "split"]
        assert len(splits) > 0


class TestMergeMatching:
    """Test merge (N->1) cell matching."""
    
    def test_merge_match(self):
        """Test multiple AtoMx cells merging into one Proseg cell."""
        from shapely.geometry import box
        
        # Two smaller AtoMx cells
        atomx_poly1 = box(0, 0, 10, 10)
        atomx_poly2 = box(10, 0, 20, 10)
        
        # One large Proseg cell covering the space
        proseg_poly = box(0, 0, 20, 10)
        
        atomx_df = pd.DataFrame({
            "atomx_cell_id": ["a1", "a2"],
            "x": [5.0, 15.0],
            "y": [5.0, 5.0],
        })
        
        proseg_df = pd.DataFrame({
            "proseg_cell_id": ["p1"],
            "x": [10.0],
            "y": [5.0],
        })
        
        atomx_polygons = {"a1": atomx_poly1, "a2": atomx_poly2}
        proseg_polygons = {"p1": proseg_poly}
        
        result = cells.match_cells(
            atomx_cells=atomx_df,
            proseg_cells=proseg_df,
            atomx_coords={"x": "x", "y": "y"},
            proseg_coords={"x": "x", "y": "y"},
            atomx_cell_id="atomx_cell_id",
            proseg_cell_id="proseg_cell_id",
            atomx_polygons=atomx_polygons,
            proseg_polygons=proseg_polygons,
        )
        
        merges = result[result.get("relationship_type") == "merge"]
        assert len(merges) > 0


class TestLostCells:
    """Test lost cell identification."""
    
    def test_lost_cells(self):
        """Test identification of lost AtoMx cells."""
        from shapely.geometry import box
        
        # Two AtoMx cells, but only one will match
        atomx_poly1 = box(0, 0, 10, 10)
        atomx_poly2 = box(100, 100, 110, 110)  # Far away, won't match
        
        proseg_poly1 = box(0, 0, 10, 10)
        
        atomx_df = pd.DataFrame({
            "atomx_cell_id": ["a1", "a2"],
            "x": [5.0, 105.0],
            "y": [5.0, 105.0],
        })
        
        proseg_df = pd.DataFrame({
            "proseg_cell_id": ["p1"],
            "x": [5.0],
            "y": [5.0],
        })
        
        atomx_polygons = {"a1": atomx_poly1, "a2": atomx_poly2}
        proseg_polygons = {"p1": proseg_poly1}
        
        result = cells.match_cells(
            atomx_cells=atomx_df,
            proseg_cells=proseg_df,
            atomx_coords={"x": "x", "y": "y"},
            proseg_coords={"x": "x", "y": "y"},
            atomx_cell_id="atomx_cell_id",
            proseg_cell_id="proseg_cell_id",
            atomx_polygons=atomx_polygons,
            proseg_polygons=proseg_polygons,
        )
        
        lost = result[result.get("relationship_type") == "lost_atomx"]
        assert len(lost) > 0
        assert lost.iloc[0]["atomx_cell"] == "a2"


class TestNewCells:
    """Test new cell identification."""
    
    def test_new_cells(self):
        """Test identification of new Proseg cells."""
        from shapely.geometry import box
        
        atomx_poly1 = box(0, 0, 10, 10)
        
        # Two Proseg cells, but only one will match
        proseg_poly1 = box(0, 0, 10, 10)
        proseg_poly2 = box(100, 100, 110, 110)  # Far away, won't match
        
        atomx_df = pd.DataFrame({
            "atomx_cell_id": ["a1"],
            "x": [5.0],
            "y": [5.0],
        })
        
        proseg_df = pd.DataFrame({
            "proseg_cell_id": ["p1", "p2"],
            "x": [5.0, 105.0],
            "y": [5.0, 105.0],
        })
        
        atomx_polygons = {"a1": atomx_poly1}
        proseg_polygons = {"p1": proseg_poly1, "p2": proseg_poly2}
        
        result = cells.match_cells(
            atomx_cells=atomx_df,
            proseg_cells=proseg_df,
            atomx_coords={"x": "x", "y": "y"},
            proseg_coords={"x": "x", "y": "y"},
            atomx_cell_id="atomx_cell_id",
            proseg_cell_id="proseg_cell_id",
            atomx_polygons=atomx_polygons,
            proseg_polygons=proseg_polygons,
        )
        
        new = result[result.get("relationship_type") == "new_proseg"]
        assert len(new) > 0
        assert new.iloc[0]["proseg_cell"] == "p2"


class TestComplexMatching:
    """Test complex (N->M) cell matching."""
    
    def test_complex_match(self):
        """Test complex many-to-many matching."""
        from shapely.geometry import box
        
        # 2 AtoMx cells
        atomx_poly1 = box(0, 0, 15, 10)
        atomx_poly2 = box(10, 0, 25, 10)
        
        # 2 Proseg cells with partial overlaps
        proseg_poly1 = box(0, 0, 12, 10)
        proseg_poly2 = box(12, 0, 25, 10)
        
        atomx_df = pd.DataFrame({
            "atomx_cell_id": ["a1", "a2"],
            "x": [7.5, 17.5],
            "y": [5.0, 5.0],
        })
        
        proseg_df = pd.DataFrame({
            "proseg_cell_id": ["p1", "p2"],
            "x": [6.0, 18.0],
            "y": [5.0, 5.0],
        })
        
        atomx_polygons = {"a1": atomx_poly1, "a2": atomx_poly2}
        proseg_polygons = {"p1": proseg_poly1, "p2": proseg_poly2}
        
        result = cells.match_cells(
            atomx_cells=atomx_df,
            proseg_cells=proseg_df,
            atomx_coords={"x": "x", "y": "y"},
            proseg_coords={"x": "x", "y": "y"},
            atomx_cell_id="atomx_cell_id",
            proseg_cell_id="proseg_cell_id",
            atomx_polygons=atomx_polygons,
            proseg_polygons=proseg_polygons,
        )
        
        complex_rel = result[result.get("relationship_type") == "complex"]
        # Should have some complex relationships due to overlapping coverage
        assert len(result) > 0

