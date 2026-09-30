#!/usr/bin/python3

import argparse
import re
import os
import pandas as pd
import matplotlib.pyplot as plt


def parse_busco_summary(file_path):
    """Parses the short_summary.txt file from BUSCO output (supports v4 and v5)."""
    with open(file_path, 'r') as f:
        content = f.read()

    # Search for the statistics line
    match = re.search(r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)', content)

    if not match:
        raise ValueError(f"Failed to parse file {file_path}.")

    c, s, d, f, m, n = match.groups()
    # Extract name from file for default use
    assembly_name = os.path.basename(file_path).replace('short_summary.', '').replace('.txt', '')

    return {
        'Assembly': assembly_name,
        'Single': float(s),
        'Duplicated': float(d),
        'Fragmented': float(f),
        'Missing': float(m),
        'Total BUSCOs': int(n)
    }


def plot_busco(data_df, output_file, y_labels):
    """Builds and saves a horizontal stacked bar chart."""
    # Colors similar to official BUSCO plots
    colors = {
        'Single': '#1E90FF',  # Blue
        'Duplicated': '#87CEFA',  # Light blue
        'Fragmented': '#FF8C00',  # Orange
        'Missing': '#CD5C5C'  # Red-brown
    }

    fig, ax = plt.subplots(figsize=(10, max(2, len(data_df) * 1.2)))
    bar_width = 0.6
    y_pos = range(len(data_df))

    left = [0] * len(data_df)
    for status in ['Single', 'Duplicated', 'Fragmented', 'Missing']:
        values = data_df[status].values
        ax.barh(y_pos, values, left=left, color=colors[status], label=status, edgecolor='white', height=bar_width)

        for i, v in enumerate(values):
            if v > 2:
                ax.text(left[i] + v / 2, i, f"{v:.1f}%", va='center', ha='center', color='black', fontsize=9,
                        fontweight='bold')
        left = [left[j] + values[j] for j in range(len(values))]

    # Set Y-axis labels
    ax.set_yticks(y_pos)

    # Use custom labels if provided and their count matches
    if y_labels and len(y_labels) == len(data_df):
        ax.set_yticklabels(y_labels, fontsize=11)
    else:
        if y_labels and len(y_labels) != len(data_df):
            print("Warning: The number of titles in --title does not match the number of files. Using filenames instead.")
        ax.set_yticklabels(data_df['Assembly'], fontsize=11)

    ax.set_xlim(0, 100)
    ax.set_xlabel('BUSCO (%)', fontsize=12)

    # Remove spines
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)

    # Legend in the upper right, but placed outside the plot area to avoid overlapping data
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), ncol=1, frameon=False)

    # tight_layout accounts for the external legend when saving
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot successfully saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Visualize BUSCO results from short_summary.txt")
    parser.add_argument('files', nargs='+', help="Paths to short_summary.txt files (multiple allowed)")
    parser.add_argument('-o', '--output', default='busco_plot.png',
                        help="Output filename (e.g., busco.png, busco.pdf)")
    parser.add_argument('--title', nargs='+', default=None,
                        help="Custom labels for the Y-axis in the order of passed files (e.g., --title \"Assembly A\" \"Assembly B\")")

    args = parser.parse_args()

    data = []
    for file_path in args.files:
        try:
            data.append(parse_busco_summary(file_path))
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    if not data:
        print("No data to plot. Exiting.")
        return

    df = pd.DataFrame(data)
    plot_busco(df, args.output, args.title)


if __name__ == '__main__':
    main()