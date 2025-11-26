#!/usr/bin/env python3
"""
Visualize confidence score distributions and component correlations.

Creates:
1. Histogram of confidence scores with tier annotations
2. Component score distributions
3. Correlation matrix heatmap
4. Scatter plots: components vs final score
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys


def plot_confidence_histogram(df, output_prefix):
    """Plot histogram of confidence scores with tier annotations."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot histogram
    n, bins, patches = ax.hist(df['confidence_score'], bins=30, edgecolor='black', alpha=0.7, color='skyblue')

    # Color bars by tier
    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i+1]) / 2
        if bin_center >= 0.8:
            patch.set_facecolor('#27ae60')  # Green for high confidence
        elif bin_center >= 0.5:
            patch.set_facecolor('#f39c12')  # Orange for medium
        else:
            patch.set_facecolor('#e74c3c')  # Red for low

    # Add tier boundaries
    ax.axvline(0.8, color='darkgreen', linestyle='--', linewidth=2, label='High confidence (≥0.8)')
    ax.axvline(0.5, color='darkorange', linestyle='--', linewidth=2, label='Medium confidence (≥0.5)')

    # Add statistics
    mean_score = df['confidence_score'].mean()
    median_score = df['confidence_score'].median()
    ax.axvline(mean_score, color='blue', linestyle='-', linewidth=2, label=f'Mean ({mean_score:.3f})')
    ax.axvline(median_score, color='purple', linestyle='-', linewidth=2, label=f'Median ({median_score:.3f})')

    # Labels and formatting
    ax.set_xlabel('Confidence Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_title('Crossover Confidence Score Distribution', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    # Add tier counts as text
    high_count = (df['confidence_score'] >= 0.8).sum()
    med_count = ((df['confidence_score'] >= 0.5) & (df['confidence_score'] < 0.8)).sum()
    low_count = (df['confidence_score'] < 0.5).sum()

    textstr = f'High: {high_count} ({high_count/len(df)*100:.1f}%)\nMedium: {med_count} ({med_count/len(df)*100:.1f}%)\nLow: {low_count} ({low_count/len(df)*100:.1f}%)'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.98, 0.97, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_confidence_histogram.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Created: {output_prefix}_confidence_histogram.png")


def plot_component_distributions(df, output_prefix):
    """Plot distributions of all component scores."""
    components = ['mapping_consensus', 'segment_quality', 'mapping_quality',
                 'parent_discrimination', 'crossover_resolution']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, component in enumerate(components):
        ax = axes[i]
        ax.hist(df[component], bins=20, edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xlabel('Score', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(component.replace('_', ' ').title(), fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        # Add mean line
        mean_val = df[component].mean()
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.3f}')
        ax.legend(fontsize=8)

    # Remove empty subplot
    fig.delaxes(axes[5])

    plt.suptitle('Component Score Distributions', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_component_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Created: {output_prefix}_component_distributions.png")


def plot_correlation_matrix(df, output_prefix):
    """Plot correlation matrix of all scores."""
    score_cols = ['confidence_score', 'mapping_consensus', 'segment_quality',
                 'mapping_quality', 'parent_discrimination', 'crossover_resolution']

    corr_matrix = df[score_cols].corr()

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='coolwarm', center=0,
                square=True, linewidths=1, cbar_kws={"shrink": 0.8}, ax=ax)

    # Rotate labels
    ax.set_xticklabels([col.replace('_', '\n') for col in score_cols], rotation=45, ha='right')
    ax.set_yticklabels([col.replace('_', ' ').title() for col in score_cols], rotation=0)

    ax.set_title('Confidence Score Component Correlations', fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_correlation_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Created: {output_prefix}_correlation_matrix.png")


def plot_component_vs_confidence(df, output_prefix):
    """Scatter plots: each component vs final confidence score."""
    components = ['mapping_consensus', 'segment_quality', 'mapping_quality',
                 'parent_discrimination', 'crossover_resolution']
    weights = [0.30, 0.25, 0.20, 0.15, 0.10]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, (component, weight) in enumerate(zip(components, weights)):
        ax = axes[i]

        # Scatter plot with color by confidence tier
        colors = df['confidence_score'].apply(lambda x: '#27ae60' if x >= 0.8 else ('#f39c12' if x >= 0.5 else '#e74c3c'))
        ax.scatter(df[component], df['confidence_score'], c=colors, alpha=0.6, s=30, edgecolor='black', linewidth=0.5)

        # Add regression line
        z = np.polyfit(df[component], df['confidence_score'], 1)
        p = np.poly1d(z)
        x_line = np.linspace(df[component].min(), df[component].max(), 100)
        ax.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2)

        # Calculate correlation
        corr = df[component].corr(df['confidence_score'])

        ax.set_xlabel(component.replace('_', ' ').title(), fontsize=10, fontweight='bold')
        ax.set_ylabel('Confidence Score', fontsize=10, fontweight='bold')
        ax.set_title(f'{component.replace("_", " ").title()}\n(Weight: {weight:.0%}, Corr: {corr:.3f})',
                    fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)

    # Remove empty subplot and add legend
    fig.delaxes(axes[5])
    legend_ax = fig.add_subplot(2, 3, 6)
    legend_ax.axis('off')
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#27ae60', edgecolor='black', label='High (≥0.8)'),
        Patch(facecolor='#f39c12', edgecolor='black', label='Medium (0.5-0.8)'),
        Patch(facecolor='#e74c3c', edgecolor='black', label='Low (<0.5)')
    ]
    legend_ax.legend(handles=legend_elements, loc='center', fontsize=12, title='Confidence Tier', title_fontsize=13)

    plt.suptitle('Component Scores vs Final Confidence', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_component_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Created: {output_prefix}_component_scatter.png")


def plot_region_comparison(df, output_prefix):
    """Compare ARMS vs CEN confidence scores."""
    if 'region_type' not in df.columns:
        print("Skipping region comparison (no region_type column)")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Histogram comparison
    ax1 = axes[0]
    arms_scores = df[df['region_type'] == 'ARMS']['confidence_score']
    cen_scores = df[df['region_type'] == 'CEN']['confidence_score']

    ax1.hist(arms_scores, bins=20, alpha=0.6, label=f'ARMS (n={len(arms_scores)})', color='#3498db', edgecolor='black')
    ax1.hist(cen_scores, bins=20, alpha=0.6, label=f'CEN (n={len(cen_scores)})', color='#e74c3c', edgecolor='black')
    ax1.axvline(arms_scores.mean(), color='#2471a3', linestyle='--', linewidth=2, label=f'ARMS mean: {arms_scores.mean():.3f}')
    ax1.axvline(cen_scores.mean(), color='#c0392b', linestyle='--', linewidth=2, label=f'CEN mean: {cen_scores.mean():.3f}')

    ax1.set_xlabel('Confidence Score', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax1.set_title('Confidence Score by Region Type', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, axis='y')

    # Box plot comparison
    ax2 = axes[1]
    region_data = [arms_scores, cen_scores]
    bp = ax2.boxplot(region_data, labels=['ARMS', 'CEN'], patch_artist=True, widths=0.6)

    colors = ['#3498db', '#e74c3c']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax2.set_ylabel('Confidence Score', fontsize=12, fontweight='bold')
    ax2.set_title('Confidence Score Distribution\nby Region Type', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')

    # Add statistical test result
    from scipy import stats
    if len(arms_scores) > 0 and len(cen_scores) > 0:
        stat, pval = stats.mannwhitneyu(arms_scores, cen_scores, alternative='two-sided')
        ax2.text(0.5, 0.05, f'Mann-Whitney U test\np-value: {pval:.4f}',
                transform=ax2.transAxes, ha='center', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_region_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Created: {output_prefix}_region_comparison.png")


def main():
    parser = argparse.ArgumentParser(description='Visualize confidence score distributions')
    parser.add_argument('--input', required=True, help='Input TSV file with confidence scores')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')

    args = parser.parse_args()

    # Load data
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error loading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Check required columns
    required_cols = ['confidence_score', 'mapping_consensus', 'segment_quality',
                    'mapping_quality', 'parent_discrimination', 'crossover_resolution']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)

    if len(df) == 0:
        print("Warning: No data to plot", file=sys.stderr)
        sys.exit(0)

    print(f"Loaded {len(df)} confidence scores")
    print("Generating visualizations...")

    # Create plots
    sns.set_style("white")
    plot_confidence_histogram(df, args.output_prefix)
    plot_component_distributions(df, args.output_prefix)
    plot_correlation_matrix(df, args.output_prefix)
    plot_component_vs_confidence(df, args.output_prefix)
    plot_region_comparison(df, args.output_prefix)

    print("\nVisualization complete!")
    print(f"Created 5 plots with prefix: {args.output_prefix}")


if __name__ == '__main__':
    main()
