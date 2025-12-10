"""
Survival Analysis Module for Cancer Transcriptomics
Implements Kaplan-Meier curves, Cox regression, and risk stratification

References:
    - Therneau, T.M., Grambsch, P.M. (2000). Modeling Survival Data
    - Davidson-Pilon, C. (2019). lifelines: survival analysis in Python

Author: Dorra Rjaibi
Date: 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler
from scipy import stats

import warnings
warnings.filterwarnings('ignore')


class SurvivalAnalyzer:
    """
    Comprehensive survival analysis for cancer patients.
    
    Implements:
    - Kaplan-Meier survival curves
    - Log-rank test for survival differences
    - Cox proportional hazards regression
    - Risk score calculation and stratification
    - Time-dependent ROC analysis
    - Prognostic gene identification
    
    Attributes:
        expression_data (pd.DataFrame): Gene expression matrix (genes x samples)
        clinical_data (pd.DataFrame): Clinical metadata with survival info
        time_col (str): Column name for survival time
        event_col (str): Column name for event occurrence
    """
    
    def __init__(
        self,
        expression_data: pd.DataFrame,
        clinical_data: pd.DataFrame,
        time_col: str = 'overall_survival_days',
        event_col: str = 'deceased'
    ):
        """
        Initialize the SurvivalAnalyzer.
        
        Args:
            expression_data: Normalized gene expression matrix
            clinical_data: Clinical data with survival information
            time_col: Name of column containing survival time
            event_col: Name of column containing event status (0/1 or True/False)
        """
        self.expression_data = expression_data
        self.clinical_data = clinical_data
        self.time_col = time_col
        self.event_col = event_col
        
        # Ensure common samples
        common_samples = self.expression_data.columns.intersection(
            self.clinical_data.index
        )
        self.expression_data = self.expression_data[common_samples]
        self.clinical_data = self.clinical_data.loc[common_samples]
        
        # Convert event column to numeric if needed
        if self.clinical_data[event_col].dtype == bool:
            self.clinical_data[event_col] = self.clinical_data[event_col].astype(int)
        
        print(f"Initialized SurvivalAnalyzer with {len(common_samples)} samples")
        print(f"Events: {self.clinical_data[event_col].sum()}")
        print(f"Censored: {len(common_samples) - self.clinical_data[event_col].sum()}")
    
    def analyze_gene(
        self,
        gene: str,
        stratification_method: str = 'median',
        output_path: Optional[str] = None
    ) -> Dict:
        """
        Perform survival analysis for a single gene.
        
        Stratifies patients into high and low expression groups and compares
        survival curves using log-rank test.
        
        Args:
            gene: Gene symbol to analyze
            stratification_method: Method to split patients ('median', 'quartile', 'optimal')
            output_path: Path to save Kaplan-Meier plot (optional)
        
        Returns:
            Dictionary containing:
                - pvalue: Log-rank test p-value
                - hazard_ratio: Hazard ratio (high vs low)
                - median_survival_high: Median survival for high expression group
                - median_survival_low: Median survival for low expression group
                - n_high: Number of patients in high expression group
                - n_low: Number of patients in low expression group
        """
        if gene not in self.expression_data.index:
            raise ValueError(f"Gene {gene} not found in expression data")
        
        # Get gene expression values
        gene_expr = self.expression_data.loc[gene]
        
        # Stratify patients
        if stratification_method == 'median':
            threshold = gene_expr.median()
            high_expr = gene_expr > threshold
        elif stratification_method == 'quartile':
            q25 = gene_expr.quantile(0.25)
            q75 = gene_expr.quantile(0.75)
            # Use top and bottom quartiles only
            high_expr = gene_expr >= q75
            low_expr = gene_expr <= q25
            # Exclude middle 50%
            mask = high_expr | low_expr
            gene_expr = gene_expr[mask]
            high_expr = high_expr[mask]
            self.clinical_data = self.clinical_data.loc[mask]
        elif stratification_method == 'optimal':
            # Find optimal cutpoint using log-rank test
            threshold = self._find_optimal_cutpoint(gene_expr)
            high_expr = gene_expr > threshold
        else:
            raise ValueError(f"Unknown stratification method: {stratification_method}")
        
        # Create survival data for each group
        high_group = self.clinical_data[high_expr].copy()
        high_group['expression_group'] = 'High'
        
        low_group = self.clinical_data[~high_expr].copy()
        low_group['expression_group'] = 'Low'
        
        # Perform log-rank test
        results = logrank_test(
            durations_A=high_group[self.time_col],
            durations_B=low_group[self.time_col],
            event_observed_A=high_group[self.event_col],
            event_observed_B=low_group[self.event_col]
        )
        
        # Fit Kaplan-Meier curves
        kmf_high = KaplanMeierFitter()
        kmf_high.fit(
            durations=high_group[self.time_col],
            event_observed=high_group[self.event_col],
            label='High Expression'
        )
        
        kmf_low = KaplanMeierFitter()
        kmf_low.fit(
            durations=low_group[self.time_col],
            event_observed=low_group[self.event_col],
            label='Low Expression'
        )
        
        # Calculate hazard ratio using Cox regression
        cox_data = pd.concat([high_group, low_group])
        cox_data['expression_binary'] = (
            cox_data['expression_group'] == 'High'
        ).astype(int)
        
        cph = CoxPHFitter()
        cph.fit(
            cox_data[[self.time_col, self.event_col, 'expression_binary']],
            duration_col=self.time_col,
            event_col=self.event_col
        )
        
        hazard_ratio = np.exp(cph.params_['expression_binary'])
        
        # Plot Kaplan-Meier curves if output path provided
        if output_path:
            self.plot_km_curve_comparison(
                kmf_high, kmf_low,
                gene=gene,
                pvalue=results.p_value,
                output_path=output_path
            )
        
        return {
            'gene': gene,
            'pvalue': results.p_value,
            'test_statistic': results.test_statistic,
            'hazard_ratio': hazard_ratio,
            'median_survival_high': kmf_high.median_survival_time_,
            'median_survival_low': kmf_low.median_survival_time_,
            'n_high': len(high_group),
            'n_low': len(low_group),
            'threshold': threshold if stratification_method != 'quartile' else None
        }
    
    def _find_optimal_cutpoint(self, gene_expr: pd.Series) -> float:
        """
        Find optimal expression cutpoint using log-rank test.
        
        Tests multiple cutpoints and selects the one with minimum p-value.
        
        Args:
            gene_expr: Gene expression values
        
        Returns:
            Optimal threshold value
        """
        # Test cutpoints from 25th to 75th percentile
        cutpoints = np.percentile(gene_expr, np.arange(25, 76, 5))
        
        best_pval = 1.0
        best_cutpoint = gene_expr.median()
        
        for cutpoint in cutpoints:
            high_expr = gene_expr > cutpoint
            
            # Need reasonable sample sizes in each group
            if high_expr.sum() < 10 or (~high_expr).sum() < 10:
                continue
            
            high_surv = self.clinical_data[high_expr]
            low_surv = self.clinical_data[~high_expr]
            
            try:
                results = logrank_test(
                    durations_A=high_surv[self.time_col],
                    durations_B=low_surv[self.time_col],
                    event_observed_A=high_surv[self.event_col],
                    event_observed_B=low_surv[self.event_col]
                )
                
                if results.p_value < best_pval:
                    best_pval = results.p_value
                    best_cutpoint = cutpoint
            except:
                continue
        
        return best_cutpoint
    
    def build_risk_model(
        self,
        genes: List[str],
        method: str = 'cox'
    ) -> Dict:
        """
        Build multi-gene risk prediction model.
        
        Uses Cox regression to combine multiple genes into a single risk score.
        
        Args:
            genes: List of gene symbols to include in model
            method: Model type ('cox' or 'lasso_cox')
        
        Returns:
            Dictionary containing:
                - model: Fitted Cox model
                - scores: Risk scores for each patient
                - coefficients: Gene coefficients
                - c_index: Concordance index (model performance)
                - genes: List of genes used
        """
        # Prepare feature matrix
        X = self.expression_data.loc[genes].T
        
        # Standardize features
        scaler = StandardScaler()
        X_scaled = pd.DataFrame(
            scaler.fit_transform(X),
            index=X.index,
            columns=X.columns
        )
        
        # Add survival data
        model_data = X_scaled.copy()
        model_data[self.time_col] = self.clinical_data[self.time_col]
        model_data[self.event_col] = self.clinical_data[self.event_col]
        
        # Fit Cox proportional hazards model
        cph = CoxPHFitter(penalizer=0.1 if method == 'lasso_cox' else 0.0)
        cph.fit(
            model_data,
            duration_col=self.time_col,
            event_col=self.event_col,
            show_progress=False
        )
        
        # Calculate risk scores
        risk_scores = cph.predict_partial_hazard(X_scaled)
        
        # Calculate concordance index
        c_index = concordance_index(
            self.clinical_data[self.time_col],
            -risk_scores,
            self.clinical_data[self.event_col]
        )
        
        print(f"\nRisk Model Performance:")
        print(f"  C-index: {c_index:.3f}")
        print(f"  Number of genes: {len(genes)}")
        
        return {
            'model': cph,
            'scores': risk_scores.values,
            'coefficients': dict(zip(genes, cph.params_[:len(genes)])),
            'c_index': c_index,
            'genes': genes,
            'scaler': scaler
        }
    
    def stratify_patients(
        self,
        risk_model: Dict,
        method: str = 'median',
        n_groups: int = 2
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Stratify patients into risk groups based on risk scores.
        
        Args:
            risk_model: Risk model from build_risk_model()
            method: Stratification method ('median', 'tertile', 'quartile', 'optimal')
            n_groups: Number of risk groups (2 or 3)
        
        Returns:
            Tuple of (high_risk_patients, low_risk_patients)
        """
        risk_scores = risk_model['scores']
        
        if method == 'median':
            threshold = np.median(risk_scores)
            high_risk = risk_scores > threshold
        elif method == 'tertile' or (method == 'quantile' and n_groups == 3):
            thresholds = np.percentile(risk_scores, [33.33, 66.67])
            high_risk = risk_scores > thresholds[1]
            low_risk = risk_scores < thresholds[0]
            # Return three groups
            mid_risk = ~(high_risk | low_risk)
            mid_group = self.clinical_data[mid_risk].copy()
            mid_group['risk_group'] = 'Medium'
            mid_group['risk_score'] = risk_scores[mid_risk]
        elif method == 'quartile' or (method == 'quantile' and n_groups == 4):
            thresholds = np.percentile(risk_scores, [25, 50, 75])
            high_risk = risk_scores > thresholds[2]
            low_risk = risk_scores < thresholds[0]
        else:
            raise ValueError(f"Unknown stratification method: {method}")
        
        high_risk_patients = self.clinical_data[high_risk].copy()
        high_risk_patients['risk_group'] = 'High'
        high_risk_patients['risk_score'] = risk_scores[high_risk]
        
        low_risk_patients = self.clinical_data[low_risk].copy()
        low_risk_patients['risk_group'] = 'Low'
        low_risk_patients['risk_score'] = risk_scores[low_risk]
        
        print(f"\nRisk Stratification ({method}):")
        print(f"  High risk: {len(high_risk_patients)} patients")
        print(f"  Low risk: {len(low_risk_patients)} patients")
        
        if method in ['tertile', 'quartile']:
            print(f"  Medium risk: {len(mid_group)} patients")
            return high_risk_patients, mid_group, low_risk_patients
        
        return high_risk_patients, low_risk_patients
    
    def plot_km_curve_comparison(
        self,
        kmf_high: KaplanMeierFitter,
        kmf_low: KaplanMeierFitter,
        gene: str,
        pvalue: float,
        output_path: str
    ):
        """
        Plot Kaplan-Meier survival curves for two groups.
        
        Args:
            kmf_high: KM fitter for high expression group
            kmf_low: KM fitter for low expression group
            gene: Gene name
            pvalue: Log-rank test p-value
            output_path: Path to save figure
        """
        plt.figure(figsize=(10, 7))
        
        # Plot survival curves
        kmf_high.plot_survival_function(ci_show=True, color='red', linewidth=2)
        kmf_low.plot_survival_function(ci_show=True, color='blue', linewidth=2)
        
        plt.xlabel('Time (days)', fontsize=12, fontweight='bold')
        plt.ylabel('Survival Probability', fontsize=12, fontweight='bold')
        plt.title(
            f'Kaplan-Meier Curve: {gene}\nLog-rank p = {pvalue:.2e}',
            fontsize=14,
            fontweight='bold'
        )
        plt.legend(loc='best', fontsize=11)
        plt.grid(alpha=0.3)
        
        # Add at-risk table
        add_at_risk_table(plt.gca(), kmf_high, kmf_low)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_risk_stratification(
        self,
        high_risk: pd.DataFrame,
        low_risk: pd.DataFrame,
        output_path: str
    ):
        """
        Plot Kaplan-Meier curves for risk groups.
        
        Args:
            high_risk: High-risk patient data
            low_risk: Low-risk patient data
            output_path: Path to save figure
        """
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Kaplan-Meier curves
        ax1 = axes[0]
        plt.sca(ax1)
        
        kmf_high = KaplanMeierFitter()
        kmf_high.fit(
            durations=high_risk[self.time_col],
            event_observed=high_risk[self.event_col],
            label='High Risk'
        )
        
        kmf_low = KaplanMeierFitter()
        kmf_low.fit(
            durations=low_risk[self.time_col],
            event_observed=low_risk[self.event_col],
            label='Low Risk'
        )
        
        kmf_high.plot_survival_function(ci_show=True, color='red', linewidth=2, ax=ax1)
        kmf_low.plot_survival_function(ci_show=True, color='blue', linewidth=2, ax=ax1)
        
        # Log-rank test
        results = logrank_test(
            durations_A=high_risk[self.time_col],
            durations_B=low_risk[self.time_col],
            event_observed_A=high_risk[self.event_col],
            event_observed_B=low_risk[self.event_col]
        )
        
        ax1.set_xlabel('Time (days)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Survival Probability', fontsize=12, fontweight='bold')
        ax1.set_title(
            f'Risk Stratification\nLog-rank p = {results.p_value:.2e}',
            fontsize=14,
            fontweight='bold'
        )
        ax1.legend(loc='best', fontsize=11)
        ax1.grid(alpha=0.3)
        
        # Risk score distribution
        ax2 = axes[1]
        
        all_scores = pd.concat([
            high_risk[['risk_score']].assign(group='High Risk'),
            low_risk[['risk_score']].assign(group='Low Risk')
        ])
        
        sns.violinplot(
            data=all_scores,
            x='group',
            y='risk_score',
            palette=['red', 'blue'],
            ax=ax2
        )
        
        ax2.set_xlabel('Risk Group', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Risk Score', fontsize=12, fontweight='bold')
        ax2.set_title('Risk Score Distribution', fontsize=14, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    
    def generate_survival_report(
        self,
        gene_results: List[Dict],
        risk_model: Optional[Dict] = None,
        output_path: str = 'survival_report.html'
    ):
        """
        Generate comprehensive HTML report of survival analysis.
        
        Args:
            gene_results: List of single-gene analysis results
            risk_model: Multi-gene risk model (optional)
            output_path: Path to save HTML report
        """
        # Convert results to DataFrame
        results_df = pd.DataFrame(gene_results)
        results_df = results_df.sort_values('pvalue')
        
        # Create HTML report
        html = f"""
        <html>
        <head>
            <title>Survival Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                h1 {{ color: #2c3e50; }}
                h2 {{ color: #34495e; }}
                table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
                th {{ background-color: #3498db; color: white; }}
                tr:nth-child(even) {{ background-color: #f2f2f2; }}
                .significant {{ background-color: #d5f4e6; }}
            </style>
        </head>
        <body>
            <h1>Survival Analysis Report</h1>
            <p>Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            
            <h2>Summary Statistics</h2>
            <ul>
                <li>Total genes analyzed: {len(results_df)}</li>
                <li>Significant genes (p < 0.05): {(results_df['pvalue'] < 0.05).sum()}</li>
                <li>Median survival difference > 100 days: {
                    (abs(results_df['median_survival_high'] - results_df['median_survival_low']) > 100).sum()
                }</li>
            </ul>
            
            <h2>Top Prognostic Genes</h2>
            {results_df.head(20).to_html(index=False, classes='dataframe')}
        """
        
        if risk_model is not None:
            html += f"""
            <h2>Multi-Gene Risk Model</h2>
            <ul>
                <li>C-index: {risk_model['c_index']:.3f}</li>
                <li>Number of genes: {len(risk_model['genes'])}</li>
            </ul>
            
            <h3>Gene Coefficients</h3>
            <table>
                <tr>
                    <th>Gene</th>
                    <th>Coefficient</th>
                    <th>Effect</th>
                </tr>
            """
            
            for gene, coef in sorted(
                risk_model['coefficients'].items(),
                key=lambda x: abs(x[1]),
                reverse=True
            ):
                effect = "Protective" if coef < 0 else "Risk"
                html += f"""
                <tr>
                    <td>{gene}</td>
                    <td>{coef:.4f}</td>
                    <td>{effect}</td>
                </tr>
                """
            
            html += "</table>"
        
        html += """
        </body>
        </html>
        """
        
        with open(output_path, 'w') as f:
            f.write(html)
        
        print(f"Survival report saved to: {output_path}")


def add_at_risk_table(ax, kmf_high, kmf_low):
    """Add at-risk counts below survival curve."""
    from matplotlib.patches import Rectangle
    
    # This is a simplified version - full implementation would be more complex
    times = [0, 365, 730, 1095, 1460]  # 0, 1, 2, 3, 4 years
    
    # Calculate at-risk counts at specific timepoints
    # (Implementation details omitted for brevity)
    pass


# Example usage
if __name__ == '__main__':
    # Example with simulated data
    np.random.seed(42)
    
    # Simulated expression data (100 genes x 200 samples)
    genes = [f'GENE{i}' for i in range(1, 101)]
    samples = [f'SAMPLE{i}' for i in range(1, 201)]
    expression = pd.DataFrame(
        np.random.randn(100, 200),
        index=genes,
        columns=samples
    )
    
    # Simulated clinical data
    clinical = pd.DataFrame({
        'overall_survival_days': np.random.exponential(500, 200),
        'deceased': np.random.binomial(1, 0.4, 200)
    }, index=samples)
    
    # Initialize analyzer
    analyzer = SurvivalAnalyzer(expression, clinical)
    
    # Analyze single gene
    result = analyzer.analyze_gene('GENE1', output_path='km_curve.png')
    print(f"\nGene analysis result: p-value = {result['pvalue']:.4f}")
    
    # Build multi-gene risk model
    top_genes = [f'GENE{i}' for i in range(1, 11)]
    risk_model = analyzer.build_risk_model(top_genes)
    
    # Stratify patients
    high_risk, low_risk = analyzer.stratify_patients(risk_model)
    
    print("\nSurvival analysis complete!")
