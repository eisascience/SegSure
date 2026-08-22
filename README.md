# SegSure

**Spatially resolved detection and refinement of cell-segmentation uncertainty in image-based transcriptomics**

## Overview

SegSure is a modular framework for identifying spatial regions and cells where image-based spatial transcriptomics segmentation is uncertain or inconsistent across segmentation methods. The first version focuses on:

1. Loading and auditing AtoMx-exported CosMx data
2. Loading and auditing Proseg output
3. Harmonizing coordinate and identifier systems
4. Matching cells across segmentation methods
5. Comparing transcript-to-cell assignments
6. Identifying split/merge/loss events
7. Computing interpretable disagreement metrics
8. Producing spatial uncertainty maps
9. Creating diagnostic neighborhood panels
10. Exporting clean tables for downstream analysis

**Design Philosophy:** SegSure does not ask "Which segmentation algorithm is best?" but rather "Where does segmentation become unreliable or method-dependent?" The framework treats local segmentation disagreement as a measurable property of the tissue and keeps disagreement decomposed into interpretable components.

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
│   ├── io/                        # Data loading (AtoMx, Proseg)
│   ├── harmonize/                 # Coordinate/ID harmonization
│   ├── metrics/                   # Disagreement metrics
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
│   └── 05_plot_uncertainty.py
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
- shapely >= 1.7
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

### 5. Compare Transcripts

```bash
python scripts/03_compare_transcript_assignments.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results
```

### 6. Compute Disagreement Metrics

```bash
python scripts/04_compute_disagreement_metrics.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results
```

Output: Centroid distances, transcript metrics, neighborhood discord

### 7. Visualize Uncertainty

```bash
python scripts/05_plot_uncertainty.py \
  --config config/datasets.yaml \
  --dataset-id your_dataset_id \
  --output-root results
```

Output: Spatial maps, distribution plots, neighborhood panels

## Disagreement Metrics

SegSure computes interpretable disagreement components:

### Transcript Assignment Discord
- **Jaccard Index**: Overlap of transcript assignments (0-1)
- **Discord Rate**: 1 - Jaccard Index
- **Asymmetric Discord**: Directional discordance
- **Precision/Recall**: Per-method assignment quality

### Geometric Disagreement
- **Centroid Distance**: Distance between segmented cell centers
- **Polygon Overlap**: Intersection over Union (if geometries available)
- **Boundary Disagreement**: Hausdorff-like boundary metrics

### Neighborhood Discord
- **Neighbor Discord Density**: Fraction of neighbors with high discord
- **Local Disagreement Clustering**: Spatial distribution of discord

### Cell Relationship Types
- **One-to-One**: Consistent segmentation across methods
- **Splits**: One AtoMx cell → Multiple Proseg cells
- **Merges**: Multiple AtoMx cells → One Proseg cell
- **Lost Cells**: AtoMx cells without Proseg match
- **Newly Inferred Cells**: Proseg cells without AtoMx match

## Output Files

Results are organized in `results/`:

### `audit/`
- `audit_report_*.json` - Comprehensive data audit

### `harmonized/`
- `harmonization_info.json` - Coordinate/ID mapping info
- `atomx_metadata_harmonized.csv` - Harmonized cell metadata

### `tables/`
- `cell_matches.csv` - Matched cell pairs with distances
- `cell_relationships.json` - Split/merge/loss classifications
- `transcript_comparisons.csv` - Per-cell transcript agreement
- `centroid_distances.csv` - Geometric disagreement
- `transcript_metrics.csv` - Transcript discord metrics
- `neighborhood_metrics.csv` - Neighborhood disagreement

### `spatial_maps/`
- `uncertainty_heatmap.png` - Spatial discord map
- `discord_distribution.png` - Distribution plots
- `neighborhoods/` - High-discord diagnostic panels

## Future Extensions

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
