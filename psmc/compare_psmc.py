#!/usr/bin/python3
# The script uses Mann-Whitney U-test to measure the probability
# of one group of observations (Ne values)
# lie significantly greater than another.
# Suited for comparing PSMC results from different samples/assemblies etc.

import argparse
import numpy as np
from scipy import stats
import os


def parse_psmc(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    blocks = content.split('//')
    results = []
    theta_0 = 0.0
    n_iter = 25
    skip_val = 100

    # Extract theta_0, number of iterations, and skip from the header
    for line in content.split('\n'):
        if line.startswith('TR') and len(line.split()) >= 2:
            try:
                theta_0 = float(line.split()[1])
            except Exception:
                pass
        if line.startswith('MM') and 'n_iterations:' in line:
            try:
                n_iter = int(line.split('n_iterations:')[1].split(',')[0].strip())
            except Exception:
                pass
        if line.startswith('MM') and 'skip:' in line:
            try:
                file_skip = int(line.split('skip:')[1].split(',')[0].strip())
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

    return theta_0, skip_val, results


def main():
    parser = argparse.ArgumentParser(
        description="Compare Ne values of two PSMC assemblies at a specific time point using bootstrap replicates and the Mann-Whitney U test.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-f1', '--file1', required=True, help='First .psmc file')
    parser.add_argument('-f2', '--file2', required=True, help='Second .psmc file')
    parser.add_argument('-t', '--time', type=int, required=True, help='Target time in years before present, integer')
    parser.add_argument('-u', '--mut', type=float, default=2.5e-8, help='Mutation rate per site per generation (-u in psmc_plot)')
    parser.add_argument('-g', '--gen', type=float, default=25, help='Generation time (-g in psmc_plot)')

    args = parser.parse_args()

    theta1, skip1, res1 = parse_psmc(args.file1)
    theta2, skip2, res2 = parse_psmc(args.file2)

    if len(res1) < 2 or len(res2) < 2:
        print("Error: No bootstrap replicates found in the files. Comparison is impossible.")
        return

    n0_1 = (theta1 / skip1) / (4.0 * args.mut)
    n0_2 = (theta2 / skip2) / (4.0 * args.mut)

    def get_ne_at_time(results, n0, gen, target_year):
        scale_t = 2.0 * n0 * gen
        main_ne = None
        bs_nes = []

        for i, (t_ks, lambda_ks) in enumerate(results):
            # Scale axes to absolute values
            times = [t * scale_t for t in t_ks]
            nes = [l * n0 for l in lambda_ks]

            # Linear interpolation of Ne value for the given time
            ne_val = np.interp(target_year, times, nes)

            if i == 0:
                main_ne = ne_val
            else:
                bs_nes.append(ne_val)

        return main_ne, bs_nes

    # Get values for the specified time
    main_val1, bs_vals1 = get_ne_at_time(res1, n0_1, args.gen, args.time)
    main_val2, bs_vals2 = get_ne_at_time(res2, n0_2, args.gen, args.time)

    # Calculate 95% confidence intervals (CI)
    ci1_lower = np.percentile(bs_vals1, 2.5)
    ci1_upper = np.percentile(bs_vals1, 97.5)

    ci2_lower = np.percentile(bs_vals2, 2.5)
    ci2_upper = np.percentile(bs_vals2, 97.5)

    # Perform Mann-Whitney U test
    stat, p_value = stats.mannwhitneyu(bs_vals1, bs_vals2, alternative='two-sided')

    print("\n" + "=" * 60)
    print(f"COMPARISON OF PSMC VALUES AT t = {args.time} YEARS BEFORE PRESENT")
    print("=" * 60)
    print(f"File 1: {os.path.basename(args.file1)}")
    print(f"  Ne (main run):       {main_val1:,.0f} individuals")
    print(f"  95% CI (bootstrap):  [{ci1_lower:,.0f} - {ci1_upper:,.0f}] individuals")
    print("-" * 60)
    print(f"File 2: {os.path.basename(args.file2)}")
    print(f"  Ne (main run):       {main_val2:,.0f} individuals")
    print(f"  95% CI (bootstrap):  [{ci2_lower:,.0f} - {ci2_upper:,.0f}] individuals")
    print("=" * 60)

    print("\nSTATISTICAL CONCLUSION:")
    print(f"Mann-Whitney U-test p-value = {p_value:.5e}")

    if p_value < 0.05:
        print("The difference in Ne values is STATISTICALLY SIGNIFICANT (p < 0.05).")
    else:
        print("The difference in Ne values is NOT statistically significant (p >= 0.05).")

    if ci1_lower > ci2_upper:
        print("CI Conclusion: Ne value of File 1 is significantly HIGHER (CIs do not overlap).")
    elif ci2_lower > ci1_upper:
        print("CI Conclusion: Ne value of File 2 is significantly HIGHER (CIs do not overlap).")
    else:
        print("CI Conclusion: Confidence intervals overlap, the difference may not be significant.")


if __name__ == '__main__':
    main()