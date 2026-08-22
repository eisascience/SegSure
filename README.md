# SegSure

**Spatially resolved detection and refinement of cell-segmentation uncertainty in image-based transcriptomics**

## Overview

SegSure is a modular framework for identifying spatial regions and cells where image-based spatial transcriptomics segmentation is uncertain or inconsistent across segmentation methods. The first version focuses on:

1. Loading and auditing AtoMx-exported CosMx data
2. Loading and auditing Proseg output
3. Harmonizing coordinate and identifier systems
4. Matching cells across segmentation methods
5. **Comparing molecule-level transcript-to-cell assignments** (IMPLEMENTED)
6. **Identifying split/merge/loss events at molecule level** (IMPLEMENTED)
7. **Computing interpretable disagreement metrics** (geometry, transcript, topology) (IMPLEMENTED)
8. Producing spatial uncertainty maps
9. Creating diagnostic neighborhood panels
10. Exporting clean tables for downstream analysis

**Design Philosophy:** SegSure does not ask "Which segmentation algorithm is best?" but rather "Where does segmentation become unreliable or method-dependent?" The framework treats local segmentation disagreement as a measurable property of the tissue and keeps disagreement decomposed into interpretable components.

**Status:** The molecule-level analysis pipeline is now fully implemented with real (non-placeholder) functionality for matching molecules, classifying disagreement statuses, and computing geometry/transcript/topology metrics. Single-FOV analysis is operational; whole-slide scaling pending.

## Project Structure

```
SegSure/
├── README.md
├── pyproject.toml
├── .gitignore
│
├── config/
│   └── datasets.yaml              # Dataset configuration
│
├── segsure/
│   ├── __init__.py
│   ├── io/                        # Data loading (AtoMx, Proseg, export)
│   ├── harmonize/                 # Coordinate/ID harmonization
│   ├── molecules/                 # Molecule-level tracking and classification
│   ├── metrics/                   # Disagreement metrics (geometry, transcript, neighborhood)
│   ├── uncertainty/               # Disagreement aggregation
│   ├── plotting/                  # Spatial visualization
│   └── utils/                     # Logging, validation
│
├── scripts/
│   ├── 00_inspect_inputs.py       # Data audit and schema inspection
│   ├── 01_harmonize_atomx_proseg.py
│   ├── 02_match_cells.py
│   ├── 03_compare_transcript_assignments.py
│   ├── 04_compute_disagreement_metrics.py
│   ├── 05_plot_uncertainty.py
│   └── 06_molecule_analysis.py    # NEW: Comprehensive molecule-level analysis
│
├── notebooks/
│   └── 01_segmentation_disagreement_eda.ipynb
│
├── tests/
│   ├── test_io.py
│   ├── test_cell_matching.py
│   └── test_metrics.py
│
└── results/                       # Generated outputs
    ├── audit/
    ├── harmonized/
    ├── tables/
    ├── spatial_maps/
    ├── neighborhoods/
    └── logs/
```

## Installation

```bash
pip install -e .
```

### Dependencies

- numpy >= 1.20
- pandas >= 1.3
- scipy >= 1.7
- scikit-image >= 0.18
- anndata >= 0.8
- scanpy >= 1.9
- zarr >= 2.10
- h5py >= 3.0
- pyyaml >= 5.4
- **shapely >= 1.7** (for polygon geometry metrics)
- matplotlib >= 3.4
- seaborn >= 0.11
- tqdm >= 4.60

## Quick Start

### 1. Configure Your Dataset

Edit `config/datasets.yaml` with your data paths:

```yaml
datasets:
  your_dataset_id:
    atomx:
      root: "/path/to/atomx/data"
      expression: "expr_file.csv.gz"
      metadata: "meta_file.csv.gz"
      transcripts: "tx_file.csv.gz"
      polygons: "polygons_file.csv.gz"
      fov_positions: "fov_file.csv.gz"
    proseg:
      root: "/path/to/proseg/data"
      h5ad: "proseg_output.h5ad"
      zarr: "proseg_output.zarr"
      seurat: "proseg_output.seurat.rds"
```

### 2. Audit Input Data

```bash
python scripts/00_inspect_inputs.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results \
  --verbose
```

This creates a comprehensive audit report documenting:
- File schemas and column names
- Data completeness
- Coordinate ranges
- Identifier distributions

### 3. Harmonize Coordinate Systems

```bash
python scripts/01_harmonize_atomx_proseg.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results
```

### 4. Match Cells

```bash
python scripts/02_match_cells.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results \
  --distance-threshold 10.0
```

Output: `results/tables/cell_matches.csv`, `results/tables/cell_relationships.json`

### 5. Analyze Molecules (Primary Analysis - NEW)

```bash
python scripts/06_molecule_analysis.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --fov-id FOV_ID \
  --output-root results \
  --spatial-tolerance 1.5 \
  --verbose
```

This is the primary analysis script that performs comprehensive molecule-level disagreement analysis:
- Matches molecules using stable IDs and spatial proximity
- Classifies molecules into 8 interpretable status categories
- Computes geometry, transcript, and topology disagreement metrics
- Generates diagnostic panels
- Exports parquet, CSV, and JSON results

Output files in `results/fov_FOV_ID/tables/`:
- `molecules_analysis.parquet` - Molecule-level results
- `cell_level_summaries.csv` - Per-cell statistics
- `geometry_metrics.csv` - Polygon overlap (IoU, coverage, area)
- `neighborhood_metrics.csv` - Neighbor changes and Jaccard similarity
- `ambiguous_matches.json` - Molecules with multiple candidate matches
- `analysis_summary.json` - Comprehensive summary statistics

### 6. Compare Transcripts (Legacy - Optional)

```bash
python scripts/03_compare_transcript_assignments.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results
```

### 7. Compute Disagreement Metrics (Legacy - Optional)

```bash
python scripts/04_compute_disagreement_metrics.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results
```

Output: Centroid distances, transcript metrics, neighborhood discord

### 8. Visualize Uncertainty (Legacy - Optional)

```bash
python scripts/05_plot_uncertainty.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results
```

Output: Spatial maps, distribution plots, neighborhood panels

## Disagreement Metrics

SegSure now computes three separate categories of interpretable disagreement components:

### Molecule-Level Transcript Assignment

Each molecule is classified into one of 8 statuses:

- **same_matched_cell**: Molecule assigned to corresponding cells (cell correspondence validated via polygon overlap)
- **changed_neighbor**: Molecule assigned to neighboring cells
- **split_related**: Molecule assigned to cell involved in a split relationship (one AtoMx → multiple Proseg)
- **merge_related**: Molecule assigned to cell involved in a merge relationship (multiple AtoMx → one Proseg)
- **atomx_assigned_proseg_unassigned**: Detected in AtoMx but not assigned in Proseg
- **atomx_unassigned_proseg_assigned**: Detected in Proseg but not assigned in AtoMx
- **changed_unrelated**: Assigned to cells with no detected relationship
- **unresolved**: Ambiguous match (multiple equally plausible candidates)

### Geometric Disagreement

Computed from polygon geometry (when available):
- **Intersection over Union (IoU)**: Overlap fraction (0-1)
- **Coverage fractions**: Percent of AtoMx covered by Proseg and vice versa
- **Area metrics**: Absolute difference and ratio of cell areas
- **Centroid displacement**: Euclidean distance between cell centers
- **Boundary metrics**: Hausdorff distance and mean boundary distance (when Shapely available)

### Neighborhood Discord

Computed from cell-cell relationships:
- **Neighbor count**: Number of neighbors detected by each method
- **Shared neighbors**: Neighbors present in both methods (after correspondence mapping)
- **Gained/lost neighbors**: Method-specific neighbors
- **Neighbor Jaccard similarity**: Overlap of mapped neighbor sets
- **Local split density**: Fraction of neighbors involved in splits
- **Local merge density**: Fraction of neighbors involved in merges
- **Neighbor discord density**: Fraction of neighbors with high discord

### Cell Relationship Types

- **One-to-One**: Consistent segmentation across methods
- **Splits**: One AtoMx cell → Multiple Proseg cells
- **Merges**: Multiple AtoMx cells → One Proseg cell
- **Lost Cells**: AtoMx cells without Proseg match
- **Newly Inferred Cells**: Proseg cells without AtoMx match

## Output Files

Results are organized in `results/fov_FOV_ID/tables/`:

### Molecule-Level Data
- `molecules_analysis.parquet` - All molecules with:
  - Stable ID, gene, FOV, coordinates (x, y, z)
  - Cell assignments from both methods
  - Match method (stable_id or spatial_proximity)
  - Classification status
  - Spatial match distance

### Cell-Level Summaries
- `cell_level_summaries.csv` - Per-cell statistics:
  - Total molecules
  - Consistent/changed/unresolved counts
  - Fraction changed
  - Breakdown by status category

### Disagreement Metrics
- `geometry_metrics.csv` - Polygon-based metrics (IoU, coverage, area)
- `neighborhood_metrics.csv` - Neighbor set changes and Jaccard similarity
- `transcript_metrics.csv` - Per-cell transcript overlap (if computed)

### Diagnostic Information
- `ambiguous_matches.json` - Molecules with multiple candidate matches
- `analysis_summary.json` - Comprehensive analysis summary
- `neighborhoods/` - Diagnostic panel images (when generated)

## Molecule Matching Strategy

Molecules are matched in two stages:

1. **Stable ID Matching** (first priority)
   - If a molecule has the same transcript/molecule ID in both methods
   - And appears exactly once in each method
   - Conservatively accepted as the same molecule

2. **Spatial Proximity Matching** (fallback)
   - For unmatched molecules, group by gene
   - Find candidates within `spatial_tolerance` (default 1.5 μm)
   - Single candidate → matched
   - Multiple candidates → flagged as ambiguous
   - No candidates → left unmatched

Ambiguous matches are exported for manual review and NOT silently resolved.

## Current Limitations & Scope

The present implementation is optimized for **single-FOV analysis**:

- **One FOV per run**: Scripts are designed to analyze a single field-of-view at a time to manage memory and complexity
- **Requires pre-computed cell matches**: The analysis script expects `results/tables/cell_matches.csv` from step 2 (02_match_cells.py)
- **Polygon geometry optional**: Geometry metrics require Shapely library and polygon vertex data. If unavailable, those metrics return None
- **Method-specific identifiers**: Cell and molecule IDs are method-specific. Correspondences are established via the cell_matches graph
- **No model training**: This version does not train classifiers or uncertainty models, only detects and measures disagreement
- **Spatial matching tolerance**: Set to 1.5 μm by default; may need adjustment for different tissue types

### Not Yet Implemented (Whole-Slide Scaling)

- Chunked processing for whole-slide datasets
- Streaming/memory-mapped data access
- Batch processing scripts
- Aggregation of per-FOV results across slide
- Slide-level uncertainty maps

These are deferred to a future version after single-FOV validation is complete.

While this first version focuses on transparent, decomposed metrics, future versions may add:
- Cell morphology analysis
- Cell-type coherence metrics
- Niche information integration
- Probabilistic transcript ownership
- Viral transcript assignment confidence
- Targeted local correction algorithms
- Machine learning uncertainty scores

All future components will maintain the modular architecture and interpretability focus.

## Testing

Run tests with pytest:

```bash
pip install pytest pytest-cov
pytest tests/ -v
```

## Data Privacy

SegSure does NOT:
- Modify source data
- Copy large data files into the repository
- Assume specific file naming conventions
- Require a specific version of input software

Data paths are configured externally via `config/datasets.yaml`.

## Citation

If you use SegSure, please cite:

```bibtex
@software{segsure2024,
  title={SegSure: Spatially resolved segmentation uncertainty detection},
  author={EISA Science},
  year={2024},
  url={https://github.com/eisascience/SegSure}
}
```

## Contributing

Contributions are welcome! Please follow PEP 8 style and add tests for new features.

## License

MIT License - see LICENSE file for details
