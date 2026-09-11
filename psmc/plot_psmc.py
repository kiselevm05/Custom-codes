#!/usr/bin/python3
# Script for multiple PSMC results visualisation

import argparse
import os
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def parse_psmc(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    blocks = content.split('//')
    results = []
    theta_0 = 0
    n_iter = 25
    dt = 0
    skip_val = 100

    # Extract theta_0 and number of iterations from the header
    for line in content.split('\n'):
        if line.startswith('TR') and len(line.split()) >= 2:
            try:
                theta_0 = float(line.split()[1])
            except Exception:
                pass
        if line.startswith('DT') and len(line.split()) >= 2:
            try:
                dt = float(line.split()[1])
            except Exception:
                    pass
        if line.startswith('MM') and 'n_iterations:' in line:
            try:
                n_iter_str = line.split('n_iterations:')[1].split(',')[0].strip()
                n_iter = int(n_iter_str)
            except Exception:
                pass
        if line.startswith('MM') and 'skip:' in line:
            try:
                skip_str = line.split('skip:')[1].split(',')[0].strip()
                file_skip = int(skip_str)
                skip_val = file_skip * 100
            except Exception:
                pass

    # Parse blocks with final iterations
    for block in blocks:
        if not block.strip():
            continue

        rd_val = -1
        for line in block.strip().split('\n'):
            if line.startswith('RD'):
                try:
                    rd_val = int(line.split()[1])
                except Exception:
                    pass
                break

        if rd_val == n_iter:
            t_ks = []
            lambda_ks = []
            for line in block.strip().split('\n'):
                if line.startswith('RS'):
                    parts = line.split()
                    if len(parts) >= 4:
                        t_k = float(parts[2])
                        lambda_k = float(parts[3])
                        if t_k > 0:
                            t_ks.append(t_k)
                            lambda_ks.append(lambda_k)
            if t_ks:
                results.append((t_ks, lambda_ks))

    return theta_0, dt, skip_val, results


def main():
    parser = argparse.ArgumentParser(
        description="The script plots PSMC results from multiple files considering bootstrap replicates. "
                    "It also allows calculating bootstrap intervals (Pseudo-CI) and saving them to a CSV file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-f', '--files', nargs='+', required=True,
                        help='List of paths to input PSMC files (space-separated).')
    parser.add_argument('-o', '--output', default='psmc_plot.jpg',
                        help='Output path and filename for the plot (supports .jpg, .pdf, .png).')

    parser.add_argument('-c', '--conf-int', nargs='?', const='conf-int.csv', default=None,
                        help='Calculate 95%% bootstrap intervals (Pseudo-CI) and save to the specified CSV file. '
                             'If no filename is provided, defaults to conf-int.csv. '
                             'If output path is specified, uses a directory for the output file.')

    parser.add_argument('--legend', nargs='+', default=None,
                        help='Space-separated list of legend names corresponding to the input files. '
                             'Legend is displayed only if there are 2 or more input files.')
    parser.add_argument('--legend-title', default='PSMC-files',
                        help='Legend title')

    # psmc_plot.pl calibration options
    parser.add_argument('-u', '--mut', type=float, default=2.5e-8,
                        help='Mutation rate per nucleotide per generation.')
    parser.add_argument('-g', '--gen', type=float, default=25,
                        help='Generation time in years.')

    # Axis limit options
    parser.add_argument('-x', '--min-x', type=float, default=1e4, help='Minimum X-axis value in years.')
    parser.add_argument('-X', '--max-x', type=float, default=1e7, help='Maximum X-axis value in years.')
    parser.add_argument('-y', '--min-y', type=float, default=0,
                        help='Minimum Y-axis value (x10^4 of individuals).')
    parser.add_argument('-Y', '--max-y', type=float, default=None,
                        help='Maximum Y-axis value (x10^4 of individuals).')

    parser.add_argument('--add-generations', action='store_true',
                        help='Display the number of generations on the X-axis in parentheses after the year value.')
    parser.add_argument('--format-numbers', nargs='?', const='full', default=None,
                        help='Format X-axis numbers. Use "k" for thousands of years, or no value for full comma-separated numbers.')

    args = parser.parse_args()

    plt.figure(figsize=(10, 6))
    colors = plt.cm.tab10.colors

    csv_file = None
    csv_writer = None
    csv_filepath = None
    if args.conf_int:
        if args.output:
            csv_out_dir = os.path.dirname(os.path.abspath(args.output))
            csv_base = os.path.basename(args.conf_int)
            csv_filepath = os.path.join(csv_out_dir, csv_base)
            csv_file = open(csv_filepath, 'w', newline='')
        else:
            csv_file = open(args.conf_int, 'w', newline='')
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(
            ['Source_File', 't_k_raw', 'Time_Years', 'Ne_Median(x10^4)', 'Ne_CI_2.5(x10^4)', 'Ne_CI_97.5(x10^4)'])

    for i, filepath in enumerate(args.files):
        color = colors[i % len(colors)]

        if args.legend and i < len(args.legend):
            label = args.legend[i]
        else:
            label = os.path.basename(filepath)

        theta_0, dt, skip_val, results = parse_psmc(filepath)

        if not results:
            print(f"Warning: No final results found in file {filepath}.")
            continue

        # Calibration formulas
        n0 = (theta_0 / skip_val) / (4.0 * args.mut)
        scale_t = 2.0 * n0 * args.gen

        main_t = []
        main_ne = []

        for j, (t_ks, lambda_ks) in enumerate(results):
            scaled_t = [(t + dt) * scale_t for t in t_ks]
            scaled_ne = [l * n0 / 10000 for l in lambda_ks]

            if j == 0:
                plt.plot(scaled_t, scaled_ne, color=color, alpha=1.0, linewidth=1.5, label=label)
                main_t = scaled_t
                main_ne = scaled_ne
            else:
                plt.plot(scaled_t, scaled_ne, color=color, alpha=0.1, linewidth=1.0)

        # Calculate Pseudo-CI for CSV
        if args.conf_int and len(results) > 1:
            bs_matrix = []
            for t_ks, lambda_ks in results[1:]:
                scaled_ne_bs = [l * n0 / 10000 for l in lambda_ks]
                if len(scaled_ne_bs) == len(main_t):
                    bs_matrix.append(scaled_ne_bs)

            if bs_matrix:
                bs_array = np.array(bs_matrix)
                ci_lower = np.percentile(bs_array, 2.5, axis=0)
                ci_upper = np.percentile(bs_array, 97.5, axis=0)

                for idx in range(len(main_t)):
                    csv_writer.writerow(
                        [label, results[0][0][idx], main_t[idx], main_ne[idx], ci_lower[idx], ci_upper[idx]])

    if csv_file:
        csv_file.close()
        print(f"Confidence intervals saved to: {args.conf_int}")

    # X-axis configuration
    plt.xscale('log')

    # Generation display
    if args.add_generations:
        plt.xlabel('Years (generations) before present')

        # Create a formatter function for the X-axis tick labels
        def x_formatter(x, pos):
            if x <= 0:
                return ""

            year_str = ""
            if args.format_numbers == 'k':
                val = x / 1000.0
                year_str = f"{val:.1f}" if val < 10 else f"{int(val)}"
            elif args.format_numbers == 'full':
                year_str = f"{x:.1f}" if x < 10 else f"{int(x):,}"
            else:  # None
                if not args.add_generations:
                    return ""  # Let matplotlib handle default log format
                year_str = f"{x:.1f}" if x < 10 else f"{int(x)}"

            if args.add_generations:
                gen_val = x / args.gen
                gen_str = f"{gen_val:.1f}" if gen_val < 10 else f"{int(gen_val):,}"
                return f"{year_str} ({gen_str})"

            return year_str

        # Determine X-axis label
        if args.format_numbers == 'k':
            if args.add_generations:
                x_label = 'Thousand years (generations) before present'
            else:
                x_label = 'Thousand years before present'
        else:
            if args.add_generations:
                x_label = 'Years (generations) before present'
            else:
                x_label = 'Years before present'

            # Add "log scale" to label if default formatting is used
            if args.format_numbers is None:
                x_label += ', log scale'

        plt.xlabel(x_label)

        # Apply formatter if not using default log scale
        if args.format_numbers is not None or args.add_generations:
            plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(x_formatter))
            plt.gca().xaxis.set_minor_formatter(ticker.FuncFormatter(x_formatter))


    plt.ylabel(r'Effective population size ($\times 10^4$)')

    # Apply axis limits
    plt.xlim(left=args.min_x, right=args.max_x)
    if args.min_y is not None or args.max_y is not None:
        plt.ylim(bottom=args.min_y, top=args.max_y)

    if len(args.files) >= 2:
        plt.legend(title=args.legend_title)

    plt.grid(True, which="both", ls="--", alpha=0.5)

    # Increase bottom margin so long labels with generations are not cut off
    plt.gcf().subplots_adjust(bottom=0.15)

    plt.savefig(args.output, dpi=300)
    print(f"Plot successfully saved to file: {args.output}")


if __name__ == '__main__':
    main()