# Cancer-transcriptomics-analysis-project-
# Cancer Transcriptomics Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![R 4.0+](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

## Overview

A comprehensive, production-ready bioinformatics pipeline for analyzing cancer transcriptomics data with clinical impact. This project implements state-of-the-art methods for differential gene expression analysis, survival prediction, and biomarker discovery using RNA-seq data from The Cancer Genome Atlas (TCGA) and other public repositories.

### Clinical Impact

- **Biomarker Discovery**: Identifies prognostic and predictive gene signatures
- **Patient Stratification**: Classifies patients into risk groups for personalized treatment
- **Survival Prediction**: Builds robust models for outcome prediction
- **Pathway Analysis**: Reveals dysregulated biological processes in cancer
- **Tumor Heterogeneity**: Characterizes spatial and cellular diversity within tumors

### Key Features

- Multi-method differential expression analysis (DESeq2, edgeR, limma-voom)
- Comprehensive quality control and normalization
- Survival analysis with Cox regression and Kaplan-Meier curves
- Machine learning-based risk stratification
- Pathway enrichment and gene ontology analysis
- Interactive visualizations and clinical reports
- Reproducible workflow with Docker support
- Extensive documentation and citations

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Data Sources](#data-sources)
- [Analysis Modules](#analysis-modules)
- [Results & Interpretation](#results--interpretation)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Prerequisites

- Python 3.8 or higher
- R 4.0 or higher
- conda or mamba (recommended)
- 16GB RAM minimum (32GB recommended)
- 50GB disk space

### Setup Environment

```bash
# Clone the repository
git clone https://github.com/yourusername/cancer-transcriptomics-analysis.git
cd cancer-transcriptomics-analysis

# Create conda environment
conda env create -f environment.yml
conda activate cancer-transcriptomics

# Install R dependencies
Rscript scripts/install_r_packages.R

# Verify installation
python scripts/verify_installation.py
```

### Docker Installation (Recommended)

```bash
docker pull yourusername/cancer-transcriptomics:latest
docker run -v $(pwd):/workspace -it yourusername/cancer-transcriptomics:latest
```

## Quick Start

### Example 1: Analyze TCGA Breast Cancer Data

```bash
# Download TCGA-BRCA data
python scripts/download_tcga_data.py --project TCGA-BRCA --data-type HTSeq-Counts

# Run complete analysis pipeline
python pipeline.py \
    --config configs/brca_analysis.yaml \
    --output results/brca_analysis \
    --threads 8

# Generate clinical report
python scripts/generate_report.py --input results/brca_analysis
```

### Example 2: Custom RNA-seq Analysis

```python
from cancer_transcriptomics import Pipeline, DataLoader

# Load your RNA-seq count data
data = DataLoader.from_counts_matrix(
    counts_file="data/counts.csv",
    metadata_file="data/clinical_metadata.csv"
)

# Initialize pipeline
pipeline = Pipeline(
    methods=['deseq2', 'edger', 'limma'],
    fdr_threshold=0.05,
    log2fc_threshold=1.5
)

# Run analysis
results = pipeline.run(data)

# Generate visualizations
pipeline.plot_results(results, output_dir="figures/")
```

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────┐
│                      DATA ACQUISITION                        │
│  • TCGA Download  • GEO Import  • Custom Data Upload        │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│                   QUALITY CONTROL                            │
│  • Raw Read QC    • Alignment Stats  • Sample Filtering     │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│              PREPROCESSING & NORMALIZATION                   │
│  • TMM/DESeq2 Normalization  • Batch Effect Removal         │
│  • Low Expression Filtering  • Log Transformation           │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│           DIFFERENTIAL EXPRESSION ANALYSIS                   │
│  • DESeq2 (NB model)  • edgeR (QL F-test)                   │
│  • limma-voom (Linear model)  • Consensus DEGs              │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│            FUNCTIONAL ENRICHMENT ANALYSIS                    │
│  • GO Enrichment  • KEGG Pathways  • Reactome               │
│  • Gene Set Enrichment Analysis (GSEA)                      │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│              SURVIVAL & CLINICAL ANALYSIS                    │
│  • Cox Proportional Hazards  • Kaplan-Meier Curves          │
│  • Risk Score Calculation  • Time-dependent ROC              │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│           MACHINE LEARNING CLASSIFICATION                    │
│  • LASSO Regression  • Random Forest  • SVM                 │
│  • Cross-validation  • Feature Selection                    │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│            VISUALIZATION & REPORTING                         │
│  • Volcano Plots  • Heatmaps  • Survival Curves             │
│  • Interactive Dashboard  • Clinical Report PDF             │
└─────────────────────────────────────────────────────────────┘
```

## Data Sources

### TCGA (The Cancer Genome Atlas)

The pipeline supports 33 TCGA cancer types with >10,000 patient samples:

- **Access**: Via GDC Data Portal or TCGAbiolinks R package
- **Data Types**: RNA-seq (HTSeq counts), Clinical data, Survival information
- **Reference**: [Cancer Genome Atlas Research Network, 2013]

### GEO (Gene Expression Omnibus)

- **Access**: Via GEOquery R package
- **Datasets**: Custom cancer expression studies
- **Formats**: Series Matrix, Raw counts

### Custom Data

Supports standard bioinformatics formats:
- Count matrices (CSV, TSV, Excel)
- Clinical metadata
- FASTQ files (with preprocessing module)

## Analysis Modules

### 1. Differential Expression Analysis

Implements three complementary methods following best practices:

#### DESeq2 (Love et al., 2014)
- Negative binomial model for count data
- Optimal for small sample sizes
- Handles batch effects effectively
- **Use case**: Standard RNA-seq DE analysis

#### edgeR (Robinson et al., 2010)
- Quasi-likelihood F-test
- Excellent for low-count genes
- Flexible dispersion estimation
- **Use case**: Low-expression transcripts

#### limma-voom (Law et al., 2014)
- Linear modeling with precision weights
- Fast computation on large datasets
- Robust to outliers
- **Use case**: Large cohorts, complex designs

**Consensus Approach**: Genes identified by ≥2 methods provide high-confidence DEG lists.

### 2. Survival Analysis

Kaplan-Meier survival curves and Cox proportional hazards regression:

```python
# Example survival analysis
from cancer_transcriptomics.survival import SurvivalAnalyzer

analyzer = SurvivalAnalyzer(
    expression_data=gene_expression,
    clinical_data=clinical_info,
    time_col='overall_survival',
    event_col='deceased'
)

# Single gene analysis
results = analyzer.analyze_gene('TP53')

# Multi-gene signature
signature_genes = ['BRCA1', 'BRCA2', 'TP53', 'ATM']
risk_model = analyzer.build_risk_model(signature_genes)

# Stratify patients
high_risk, low_risk = analyzer.stratify_patients(risk_model)
```

### 3. Pathway Enrichment

Identifies dysregulated biological pathways:

- **Gene Ontology (GO)**: Biological processes, molecular functions, cellular components
- **KEGG Pathways**: Cancer-related signaling cascades
- **Reactome**: Detailed pathway annotations
- **Hallmark Gene Sets**: Cancer hallmark processes

### 4. Machine Learning Models

Risk stratification using supervised learning:

- **LASSO Cox Regression**: Feature selection for survival models
- **Random Forest**: Classification of cancer subtypes
- **Support Vector Machines**: High-dimensional classification
- **Cross-validation**: Rigorous model evaluation

### 5. Tumor Microenvironment Analysis

Deconvolution of immune cell populations:

- **CIBERSORT**: Estimates 22 immune cell types
- **xCell**: 64 cell types including immune and stromal
- **TIMER**: Tumor-immune interaction analysis

## Results & Interpretation

### Output Directory Structure

```
results/
├── quality_control/
│   ├── sample_qc_report.html
│   ├── pca_plot.png
│   └── correlation_heatmap.png
├── differential_expression/
│   ├── deseq2_results.csv
│   ├── edgeR_results.csv
│   ├── limma_results.csv
│   ├── consensus_degs.csv
│   └── volcano_plots/
├── functional_enrichment/
│   ├── go_enrichment.csv
│   ├── kegg_pathways.csv
│   └── gsea_results/
├── survival_analysis/
│   ├── kaplan_meier_curves/
│   ├── cox_regression_results.csv
│   └── risk_stratification.csv
├── machine_learning/
│   ├── model_performance.csv
│   ├── roc_curves.png
│   └── feature_importance.csv
├── visualizations/
│   ├── heatmaps/
│   ├── interactive_plots/
│   └── network_diagrams/
└── clinical_report.pdf
```

### Key Metrics

- **Differentially Expressed Genes**: FDR < 0.05, |log2FC| > 1.5
- **Survival Significance**: Log-rank p < 0.05
- **Model Performance**: C-index > 0.70 for clinical utility
- **Pathway Enrichment**: FDR < 0.05

## Scientific References

This pipeline implements methods from the following key publications:

1. **DESeq2**: Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.

2. **edgeR**: Robinson, M.D., McCarthy, D.J., Smyth, G.K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139-140.

3. **limma-voom**: Law, C.W., Chen, Y., Shi, W., Smyth, G.K. (2014). voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology*, 15(2), R29.

4. **TCGA**: Cancer Genome Atlas Research Network (2013). The Cancer Genome Atlas Pan-Cancer analysis project. *Nature Genetics*, 45(10), 1113-1120.

5. **Spatial Transcriptomics**: Williams, C.G., Lee, H.J., Asatsuma, T., et al. (2022). An introduction to spatial transcriptomics for biomedical research. *Genome Medicine*, 14(1), 68.

6. **Survival Analysis**: Therneau, T.M., Grambsch, P.M. (2000). Modeling Survival Data: Extending the Cox Model. Springer.

7. **GSEA**: Subramanian, A., et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545-15550.

## Configuration

### Example Configuration File (YAML)

```yaml
# configs/cancer_analysis.yaml

project:
  name: "BRCA_Analysis"
  description: "Breast cancer transcriptomics study"
  output_dir: "results/brca"

data:
  source: "TCGA"
  project_id: "TCGA-BRCA"
  data_type: "HTSeq-Counts"
  
quality_control:
  min_counts_per_gene: 10
  min_samples_per_gene: 3
  variance_stabilization: true
  batch_correction: "ComBat"

differential_expression:
  methods: ["deseq2", "edgeR", "limma"]
  fdr_threshold: 0.05
  log2fc_threshold: 1.5
  comparison: "Tumor vs Normal"

survival_analysis:
  time_column: "overall_survival"
  event_column: "deceased"
  stratification_method: "median"
  
enrichment:
  databases: ["GO_BP", "KEGG", "REACTOME"]
  organism: "Homo sapiens"
  fdr_cutoff: 0.05

machine_learning:
  algorithms: ["lasso", "random_forest"]
  cv_folds: 5
  test_size: 0.2
  random_seed: 42

visualization:
  interactive: true
  export_formats: ["png", "pdf", "svg"]
```

## Clinical Interpretation Guide

### Biomarker Discovery Workflow

1. **Identify Consensus DEGs**: Genes detected by multiple methods
2. **Validate Clinical Relevance**: Check survival association
3. **Functional Context**: Examine pathway enrichment
4. **Literature Mining**: Cross-reference with known cancer genes
5. **Independent Validation**: Test in external cohorts

### Risk Stratification

High-risk patients identified by the model may benefit from:
- Aggressive treatment regimens
- Close monitoring schedules
- Clinical trial enrollment

Low-risk patients may be candidates for:
- De-escalated therapy
- Surveillance protocols

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
# Fork and clone the repository
git clone https://github.com/yourusername/cancer-transcriptomics-analysis.git

# Create a feature branch
git checkout -b feature/your-feature-name

# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/

# Submit a pull request
```

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run specific test module
pytest tests/test_differential_expression.py

# Generate coverage report
pytest --cov=cancer_transcriptomics tests/
```

## Docker Usage

```bash
# Build Docker image
docker build -t cancer-transcriptomics:latest .

# Run analysis in container
docker run -v $(pwd)/data:/data \
           -v $(pwd)/results:/results \
           cancer-transcriptomics:latest \
           python pipeline.py --config /data/config.yaml
```

## Troubleshooting

### Common Issues

1. **Memory Error**: Increase RAM or use downsampling options
2. **R Package Installation Fails**: Check R version and BiocManager
3. **TCGA Download Timeout**: Use `--retry` flag or manual download

See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for detailed solutions.

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{cancer_transcriptomics_2024,
  title = {Cancer Transcriptomics Analysis Pipeline},
  author = {Dorra Rjaibi},
  year = {2025},
  url = {[https://github.com/](https://github.com/dorra28/Cancer-transcriptomics-analysis-project-/edit/main/README.md)},
  version = {1.0.0}
}
```

## License

This project is licensed under the MIT License 

## Acknowledgments

- The Cancer Genome Atlas (TCGA) Research Network
- Bioconductor Community
- All contributors and users

## Contact

- **Issues**: [GitHub Issues](https://github.com/dorra28)
- **Email**: dorra.rjaibi@pasteur.utm.tn
  

---

**Disclaimer**: This software is for research purposes only and should not be used for clinical decision-making without appropriate validation and regulatory approval.
