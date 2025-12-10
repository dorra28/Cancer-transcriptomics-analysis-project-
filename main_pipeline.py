#!/usr/bin/env python3
"""
Cancer Transcriptomics Analysis Pipeline
Main orchestration script for comprehensive cancer RNA-seq analysis

Author: Dorra Rjaibi
Date: 2025
License: MIT
"""

import os
import sys
import yaml
import logging
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple

# Import custom modules
from cancer_transcriptomics.data_loader import TCGADataLoader, GEODataLoader
from cancer_transcriptomics.quality_control import QualityController
from cancer_transcriptomics.normalization import DataNormalizer
from cancer_transcriptomics.differential_expression import (
    DESeq2Analyzer, EdgeRAnalyzer, LimmaAnalyzer, ConsensusAnalyzer
)
from cancer_transcriptomics.survival_analysis import SurvivalAnalyzer
from cancer_transcriptomics.enrichment import EnrichmentAnalyzer
from cancer_transcriptomics.machine_learning import RiskStratifier
from cancer_transcriptomics.visualization import Visualizer
from cancer_transcriptomics.report_generator import ClinicalReportGenerator
from cancer_transcriptomics.utils import setup_logging, validate_config


class CancerTranscriptomicsPipeline:
    """
    Main pipeline class for comprehensive cancer transcriptomics analysis.
    
    This pipeline implements state-of-the-art methods for:
    - Multi-method differential expression analysis
    - Survival prediction and risk stratification
    - Pathway enrichment analysis
    - Biomarker discovery
    - Clinical report generation
    
    References:
        - Love et al. (2014) Genome Biology - DESeq2
        - Robinson et al. (2010) Bioinformatics - edgeR
        - Law et al. (2014) Genome Biology - limma-voom
    """
    
    def __init__(self, config_path: str):
        """
        Initialize the pipeline with configuration.
        
        Args:
            config_path: Path to YAML configuration file
        """
        self.config = self._load_config(config_path)
        self.output_dir = Path(self.config['project']['output_dir'])
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        log_file = self.output_dir / 'pipeline.log'
        self.logger = setup_logging(log_file)
        self.logger.info(f"Pipeline initialized: {self.config['project']['name']}")
        
        # Initialize components
        self.data_loader = None
        self.qc_controller = QualityController()
        self.normalizer = DataNormalizer()
        self.visualizer = Visualizer(self.output_dir / 'visualizations')
        
    def _load_config(self, config_path: str) -> Dict:
        """Load and validate configuration file."""
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        validate_config(config)
        return config
    
    def run(self):
        """Execute the complete analysis pipeline."""
        self.logger.info("=" * 80)
        self.logger.info("CANCER TRANSCRIPTOMICS ANALYSIS PIPELINE")
        self.logger.info("=" * 80)
        
        try:
            # Step 1: Data Acquisition
            self.logger.info("\n[STEP 1] DATA ACQUISITION")
            counts_data, metadata = self._load_data()
            
            # Step 2: Quality Control
            self.logger.info("\n[STEP 2] QUALITY CONTROL")
            qc_results = self._quality_control(counts_data, metadata)
            
            # Step 3: Normalization & Preprocessing
            self.logger.info("\n[STEP 3] NORMALIZATION & PREPROCESSING")
            normalized_data = self._normalize_data(
                qc_results['filtered_counts'],
                qc_results['filtered_metadata']
            )
            
            # Step 4: Differential Expression Analysis
            self.logger.info("\n[STEP 4] DIFFERENTIAL EXPRESSION ANALYSIS")
            de_results = self._differential_expression(
                normalized_data,
                qc_results['filtered_metadata']
            )
            
            # Step 5: Functional Enrichment
            self.logger.info("\n[STEP 5] FUNCTIONAL ENRICHMENT ANALYSIS")
            enrichment_results = self._enrichment_analysis(de_results)
            
            # Step 6: Survival Analysis
            self.logger.info("\n[STEP 6] SURVIVAL ANALYSIS")
            survival_results = self._survival_analysis(
                normalized_data,
                qc_results['filtered_metadata'],
                de_results['consensus_degs']
            )
            
            # Step 7: Machine Learning Classification
            self.logger.info("\n[STEP 7] MACHINE LEARNING & RISK STRATIFICATION")
            ml_results = self._machine_learning(
                normalized_data,
                qc_results['filtered_metadata'],
                de_results['consensus_degs']
            )
            
            # Step 8: Visualization
            self.logger.info("\n[STEP 8] GENERATING VISUALIZATIONS")
            self._generate_visualizations(
                normalized_data,
                de_results,
                survival_results,
                enrichment_results
            )
            
            # Step 9: Clinical Report
            self.logger.info("\n[STEP 9] GENERATING CLINICAL REPORT")
            self._generate_clinical_report(
                de_results,
                survival_results,
                enrichment_results,
                ml_results
            )
            
            self.logger.info("\n" + "=" * 80)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info(f"Results saved to: {self.output_dir}")
            self.logger.info("=" * 80)
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
            raise
    
    def _load_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Load RNA-seq counts and clinical metadata.
        
        Returns:
            Tuple of (counts_matrix, clinical_metadata)
        """
        data_config = self.config['data']
        
        if data_config['source'] == 'TCGA':
            self.logger.info(f"Loading TCGA data: {data_config['project_id']}")
            loader = TCGADataLoader(
                project_id=data_config['project_id'],
                data_type=data_config['data_type']
            )
            counts, metadata = loader.load()
            
        elif data_config['source'] == 'GEO':
            self.logger.info(f"Loading GEO data: {data_config['series_id']}")
            loader = GEODataLoader(series_id=data_config['series_id'])
            counts, metadata = loader.load()
            
        elif data_config['source'] == 'custom':
            self.logger.info("Loading custom data")
            counts = pd.read_csv(data_config['counts_file'], index_col=0)
            metadata = pd.read_csv(data_config['metadata_file'], index_col=0)
        
        else:
            raise ValueError(f"Unknown data source: {data_config['source']}")
        
        self.logger.info(f"Loaded: {counts.shape[0]} genes, {counts.shape[1]} samples")
        return counts, metadata
    
    def _quality_control(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame
    ) -> Dict:
        """
        Perform comprehensive quality control.
        
        Implements:
        - Low-count gene filtering
        - Sample outlier detection
        - Batch effect visualization
        - PCA analysis
        
        Args:
            counts: Raw count matrix
            metadata: Clinical metadata
            
        Returns:
            Dictionary with filtered data and QC metrics
        """
        qc_config = self.config['quality_control']
        qc_dir = self.output_dir / 'quality_control'
        qc_dir.mkdir(exist_ok=True)
        
        # Filter low-count genes
        self.logger.info("Filtering low-expression genes...")
        filtered_counts = self.qc_controller.filter_low_counts(
            counts,
            min_count=qc_config['min_counts_per_gene'],
            min_samples=qc_config['min_samples_per_gene']
        )
        self.logger.info(f"Retained {filtered_counts.shape[0]} genes after filtering")
        
        # Detect outlier samples
        self.logger.info("Detecting outlier samples...")
        outliers = self.qc_controller.detect_outliers(filtered_counts)
        if outliers:
            self.logger.warning(f"Detected {len(outliers)} outlier samples: {outliers}")
            # Remove outliers
            filtered_counts = filtered_counts.drop(columns=outliers)
            metadata = metadata.drop(index=outliers)
        
        # Generate QC visualizations
        self.logger.info("Generating QC plots...")
        
        # Sample correlation heatmap
        self.qc_controller.plot_sample_correlation(
            filtered_counts,
            output_path=qc_dir / 'sample_correlation.png'
        )
        
        # PCA plot
        self.qc_controller.plot_pca(
            filtered_counts,
            metadata,
            output_path=qc_dir / 'pca_plot.png'
        )
        
        # Distribution plots
        self.qc_controller.plot_count_distribution(
            filtered_counts,
            output_path=qc_dir / 'count_distribution.png'
        )
        
        # Generate QC report
        qc_report = self.qc_controller.generate_report(
            counts, filtered_counts, outliers
        )
        qc_report.to_html(qc_dir / 'qc_report.html')
        
        return {
            'filtered_counts': filtered_counts,
            'filtered_metadata': metadata,
            'outliers': outliers,
            'qc_metrics': qc_report
        }
    
    def _normalize_data(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Normalize count data and correct for batch effects.
        
        Implements:
        - Variance stabilizing transformation (VST)
        - TMM normalization
        - Batch effect correction (ComBat)
        
        Args:
            counts: Filtered count matrix
            metadata: Clinical metadata
            
        Returns:
            Normalized expression matrix
        """
        qc_config = self.config['quality_control']
        norm_dir = self.output_dir / 'normalization'
        norm_dir.mkdir(exist_ok=True)
        
        # Variance stabilization
        if qc_config.get('variance_stabilization', True):
            self.logger.info("Applying variance stabilization...")
            normalized = self.normalizer.vst_transform(counts)
        else:
            self.logger.info("Applying log2 normalization...")
            normalized = self.normalizer.log2_transform(counts)
        
        # Batch correction
        if 'batch' in metadata.columns and qc_config.get('batch_correction'):
            self.logger.info(f"Correcting batch effects using {qc_config['batch_correction']}...")
            normalized = self.normalizer.correct_batch_effects(
                normalized,
                metadata['batch'],
                method=qc_config['batch_correction']
            )
        
        # Save normalized data
        normalized.to_csv(norm_dir / 'normalized_expression.csv')
        self.logger.info(f"Normalized data saved: {normalized.shape}")
        
        return normalized
    
    def _differential_expression(
        self,
        normalized_data: pd.DataFrame,
        metadata: pd.DataFrame
    ) -> Dict:
        """
        Multi-method differential expression analysis.
        
        Implements three complementary approaches:
        1. DESeq2: Negative binomial model
        2. edgeR: Quasi-likelihood F-test
        3. limma-voom: Linear modeling
        
        Consensus DEGs are identified across methods for robustness.
        
        References:
            - Love et al. (2014) Genome Biology
            - Robinson et al. (2010) Bioinformatics
            - Law et al. (2014) Genome Biology
        
        Args:
            normalized_data: Normalized expression matrix
            metadata: Clinical metadata with group information
            
        Returns:
            Dictionary containing results from all methods and consensus
        """
        de_config = self.config['differential_expression']
        de_dir = self.output_dir / 'differential_expression'
        de_dir.mkdir(exist_ok=True)
        
        methods = de_config['methods']
        fdr_threshold = de_config['fdr_threshold']
        log2fc_threshold = de_config['log2fc_threshold']
        
        results = {}
        
        # DESeq2 Analysis
        if 'deseq2' in methods:
            self.logger.info("Running DESeq2 analysis...")
            deseq2 = DESeq2Analyzer(
                fdr_threshold=fdr_threshold,
                log2fc_threshold=log2fc_threshold
            )
            results['deseq2'] = deseq2.analyze(normalized_data, metadata)
            results['deseq2'].to_csv(de_dir / 'deseq2_results.csv')
            self.logger.info(f"DESeq2: {len(results['deseq2'])} DEGs identified")
        
        # edgeR Analysis
        if 'edgeR' in methods:
            self.logger.info("Running edgeR analysis...")
            edger = EdgeRAnalyzer(
                fdr_threshold=fdr_threshold,
                log2fc_threshold=log2fc_threshold
            )
            results['edgeR'] = edger.analyze(normalized_data, metadata)
            results['edgeR'].to_csv(de_dir / 'edgeR_results.csv')
            self.logger.info(f"edgeR: {len(results['edgeR'])} DEGs identified")
        
        # limma-voom Analysis
        if 'limma' in methods:
            self.logger.info("Running limma-voom analysis...")
            limma = LimmaAnalyzer(
                fdr_threshold=fdr_threshold,
                log2fc_threshold=log2fc_threshold
            )
            results['limma'] = limma.analyze(normalized_data, metadata)
            results['limma'].to_csv(de_dir / 'limma_results.csv')
            self.logger.info(f"limma: {len(results['limma'])} DEGs identified")
        
        # Consensus Analysis
        self.logger.info("Identifying consensus DEGs...")
        consensus = ConsensusAnalyzer()
        consensus_degs = consensus.find_consensus(
            [results[m] for m in methods if m in results],
            min_methods=2
        )
        results['consensus_degs'] = consensus_degs
        consensus_degs.to_csv(de_dir / 'consensus_degs.csv')
        self.logger.info(f"Consensus: {len(consensus_degs)} high-confidence DEGs")
        
        # Generate Venn diagram
        self.visualizer.plot_venn_diagram(
            [set(results[m].index) for m in methods if m in results],
            labels=methods,
            output_path=de_dir / 'venn_diagram.png'
        )
        
        return results
    
    def _enrichment_analysis(self, de_results: Dict) -> Dict:
        """
        Functional enrichment analysis of DEGs.
        
        Performs:
        - Gene Ontology (GO) enrichment
        - KEGG pathway analysis
        - Reactome pathway analysis
        - Gene Set Enrichment Analysis (GSEA)
        
        Reference:
            Subramanian et al. (2005) PNAS - GSEA
        
        Args:
            de_results: Differential expression results
            
        Returns:
            Dictionary of enrichment results
        """
        enrich_config = self.config['enrichment']
        enrich_dir = self.output_dir / 'functional_enrichment'
        enrich_dir.mkdir(exist_ok=True)
        
        analyzer = EnrichmentAnalyzer(
            organism=enrich_config['organism'],
            fdr_cutoff=enrich_config['fdr_cutoff']
        )
        
        # Get DEG lists
        upregulated = de_results['consensus_degs'][
            de_results['consensus_degs']['log2FoldChange'] > 0
        ].index.tolist()
        
        downregulated = de_results['consensus_degs'][
            de_results['consensus_degs']['log2FoldChange'] < 0
        ].index.tolist()
        
        results = {}
        
        # GO Enrichment
        if 'GO_BP' in enrich_config['databases']:
            self.logger.info("Performing GO enrichment...")
            results['go_up'] = analyzer.go_enrichment(upregulated, ontology='BP')
            results['go_down'] = analyzer.go_enrichment(downregulated, ontology='BP')
            
            results['go_up'].to_csv(enrich_dir / 'go_enrichment_upregulated.csv')
            results['go_down'].to_csv(enrich_dir / 'go_enrichment_downregulated.csv')
        
        # KEGG Pathways
        if 'KEGG' in enrich_config['databases']:
            self.logger.info("Performing KEGG pathway analysis...")
            results['kegg'] = analyzer.kegg_enrichment(
                upregulated + downregulated
            )
            results['kegg'].to_csv(enrich_dir / 'kegg_pathways.csv')
        
        # Reactome Pathways
        if 'REACTOME' in enrich_config['databases']:
            self.logger.info("Performing Reactome pathway analysis...")
            results['reactome'] = analyzer.reactome_enrichment(
                upregulated + downregulated
            )
            results['reactome'].to_csv(enrich_dir / 'reactome_pathways.csv')
        
        # Generate enrichment plots
        if 'go_up' in results:
            self.visualizer.plot_enrichment_bar(
                results['go_up'],
                output_path=enrich_dir / 'go_enrichment_plot.png'
            )
        
        return results
    
    def _survival_analysis(
        self,
        expression_data: pd.DataFrame,
        clinical_data: pd.DataFrame,
        deg_list: pd.DataFrame
    ) -> Dict:
        """
        Comprehensive survival analysis.
        
        Performs:
        - Kaplan-Meier survival curves
        - Cox proportional hazards regression
        - Risk stratification
        - Time-dependent ROC analysis
        
        Reference:
            Therneau & Grambsch (2000) - Cox Models
        
        Args:
            expression_data: Normalized expression matrix
            clinical_data: Clinical data with survival information
            deg_list: List of DEGs to test
            
        Returns:
            Dictionary of survival analysis results
        """
        surv_config = self.config['survival_analysis']
        surv_dir = self.output_dir / 'survival_analysis'
        surv_dir.mkdir(exist_ok=True)
        
        analyzer = SurvivalAnalyzer(
            expression_data=expression_data,
            clinical_data=clinical_data,
            time_col=surv_config['time_column'],
            event_col=surv_config['event_column']
        )
        
        results = {}
        
        # Single gene survival analysis for top DEGs
        self.logger.info("Performing single-gene survival analysis...")
        top_degs = deg_list.head(50).index.tolist()
        
        gene_survival_results = []
        for gene in top_degs:
            result = analyzer.analyze_gene(gene)
            if result['pvalue'] < 0.05:
                gene_survival_results.append(result)
                
                # Generate KM curve for significant genes
                analyzer.plot_km_curve(
                    gene,
                    output_path=surv_dir / 'kaplan_meier' / f'{gene}_km_curve.png'
                )
        
        results['gene_survival'] = pd.DataFrame(gene_survival_results)
        results['gene_survival'].to_csv(surv_dir / 'gene_survival_results.csv')
        
        # Multi-gene signature
        self.logger.info("Building multi-gene risk signature...")
        significant_genes = results['gene_survival']['gene'].tolist()
        
        if len(significant_genes) >= 3:
            risk_model = analyzer.build_risk_model(significant_genes[:10])
            results['risk_model'] = risk_model
            
            # Stratify patients
            high_risk, low_risk = analyzer.stratify_patients(
                risk_model,
                method=surv_config['stratification_method']
            )
            
            # Plot risk stratification
            analyzer.plot_risk_stratification(
                high_risk, low_risk,
                output_path=surv_dir / 'risk_stratification.png'
            )
            
            # Save risk scores
            risk_scores = pd.DataFrame({
                'patient': expression_data.columns,
                'risk_score': risk_model['scores'],
                'risk_group': ['High' if s > risk_model['threshold'] else 'Low' 
                              for s in risk_model['scores']]
            })
            risk_scores.to_csv(surv_dir / 'risk_scores.csv', index=False)
        
        return results
    
    def _machine_learning(
        self,
        expression_data: pd.DataFrame,
        clinical_data: pd.DataFrame,
        deg_list: pd.DataFrame
    ) -> Dict:
        """
        Machine learning-based risk stratification.
        
        Implements:
        - LASSO Cox regression for feature selection
        - Random Forest classification
        - Support Vector Machines
        - Cross-validation and performance evaluation
        
        Args:
            expression_data: Normalized expression matrix
            clinical_data: Clinical metadata
            deg_list: Differentially expressed genes
            
        Returns:
            Dictionary of ML results and model performance
        """
        ml_config = self.config['machine_learning']
        ml_dir = self.output_dir / 'machine_learning'
        ml_dir.mkdir(exist_ok=True)
        
        stratifier = RiskStratifier(
            cv_folds=ml_config['cv_folds'],
            test_size=ml_config['test_size'],
            random_seed=ml_config['random_seed']
        )
        
        results = {}
        
        # Prepare features (use top DEGs)
        feature_genes = deg_list.head(100).index.tolist()
        X = expression_data.loc[feature_genes].T
        y = clinical_data[self.config['survival_analysis']['event_column']]
        
        # LASSO Cox Regression
        if 'lasso' in ml_config['algorithms']:
            self.logger.info("Training LASSO Cox regression model...")
            lasso_results = stratifier.fit_lasso_cox(X, y, clinical_data)
            results['lasso'] = lasso_results
            
            # Save selected features
            selected_features = lasso_results['selected_features']
            pd.DataFrame(selected_features, columns=['gene', 'coefficient']).to_csv(
                ml_dir / 'lasso_selected_features.csv', index=False
            )
        
        # Random Forest
        if 'random_forest' in ml_config['algorithms']:
            self.logger.info("Training Random Forest classifier...")
            rf_results = stratifier.fit_random_forest(X, y)
            results['random_forest'] = rf_results
            
            # Feature importance
            self.visualizer.plot_feature_importance(
                rf_results['feature_importance'],
                output_path=ml_dir / 'rf_feature_importance.png'
            )
        
        # Model performance summary
        performance_summary = pd.DataFrame({
            'Model': list(results.keys()),
            'C-Index': [results[m]['c_index'] for m in results],
            'AUC': [results[m].get('auc', np.nan) for m in results]
        })
        performance_summary.to_csv(ml_dir / 'model_performance.csv', index=False)
        
        # Plot ROC curves
        self.visualizer.plot_roc_curves(
            results,
            output_path=ml_dir / 'roc_curves.png'
        )
        
        return results
    
    def _generate_visualizations(
        self,
        expression_data: pd.DataFrame,
        de_results: Dict,
        survival_results: Dict,
        enrichment_results: Dict
    ):
        """Generate comprehensive visualizations."""
        viz_dir = self.output_dir / 'visualizations'
        viz_dir.mkdir(exist_ok=True)
        
        # Volcano plot
        self.logger.info("Generating volcano plot...")
        self.visualizer.plot_volcano(
            de_results['consensus_degs'],
            output_path=viz_dir / 'volcano_plot.png'
        )
        
        # Heatmap of top DEGs
        self.logger.info("Generating expression heatmap...")
        top_genes = de_results['consensus_degs'].head(50).index
        self.visualizer.plot_heatmap(
            expression_data.loc[top_genes],
            output_path=viz_dir / 'deg_heatmap.png'
        )
        
        # Pathway network
        if 'kegg' in enrichment_results:
            self.logger.info("Generating pathway network...")
            self.visualizer.plot_pathway_network(
                enrichment_results['kegg'],
                output_path=viz_dir / 'pathway_network.png'
            )
    
    def _generate_clinical_report(
        self,
        de_results: Dict,
        survival_results: Dict,
        enrichment_results: Dict,
        ml_results: Dict
    ):
        """Generate comprehensive clinical report."""
        report_generator = ClinicalReportGenerator(
            output_dir=self.output_dir,
            project_name=self.config['project']['name']
        )
        
        report_generator.generate_pdf_report(
            de_results=de_results,
            survival_results=survival_results,
            enrichment_results=enrichment_results,
            ml_results=ml_results,
            output_filename='clinical_report.pdf'
        )
        
        self.logger.info("Clinical report generated successfully")


def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description="Cancer Transcriptomics Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze TCGA breast cancer data
  python pipeline.py --config configs/brca_analysis.yaml
  
  # Custom analysis with 16 threads
  python pipeline.py --config configs/custom_analysis.yaml --threads 16
  
  # Dry run to validate configuration
  python pipeline.py --config configs/test.yaml --dry-run
        """
    )
    
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to YAML configuration file'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads to use (default: 4)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Validate configuration without running analysis'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='Cancer Transcriptomics Pipeline v1.0.0'
    )
    
    args = parser.parse_args()
    
    # Set number of threads
    os.environ['OMP_NUM_THREADS'] = str(args.threads)
    
    try:
        # Initialize pipeline
        pipeline = CancerTranscriptomicsPipeline(args.config)
        
        if args.dry_run:
            print("Configuration validated successfully!")
            print(f"Project: {pipeline.config['project']['name']}")
            print(f"Output directory: {pipeline.output_dir}")
            sys.exit(0)
        
        # Run analysis
        pipeline.run()
        
    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
