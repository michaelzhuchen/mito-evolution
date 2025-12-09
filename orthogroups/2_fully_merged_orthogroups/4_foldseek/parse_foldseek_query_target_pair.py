import os
import pandas as pd
import argparse
import numpy as np
from scipy.stats import gamma

## Compute the mean bitscore and probability for all vs all search between a pair of OGs.

# Function to count files larger than 200 bytes in a directory. Files larger than 200 bytes are bona fide PDB structures, files smaller than 200 bytes are failed downloads.
def count_large_files(directory):
    return sum(
        os.path.getsize(os.path.join(directory, file)) > 200
        for file in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, file))
    )

# Main function
def compute_result(input_file, output_file, pdb_dir):
    # Get the basename of the file
    tsv_filename = os.path.basename(input_file)

    # Extract directory names from the basename
    try:
        dir1_name, rest = tsv_filename.split('_vs_')
        dir2_name = rest.split('_exhaustive_expect10.tsv')[0]  # Remove the trailing part
    except ValueError:
        print("Error: The filename format is incorrect. Expected format: dir1_vs_dir2_exhaustive_expect10.tsv")
        return

    dir1 = os.path.join(pdb_dir, 'OG_structs', dir1_name, "structs_raw")
    dir2 = os.path.join(pdb_dir, 'OG_structs', dir2_name, "structs_raw")

    # Check if directories exist
    if not (os.path.isdir(dir1) and os.path.isdir(dir2)):
        print(f"Error: One or both directories {dir1}, {dir2} do not exist.")
        return

    # Load the TSV file
    try:
        tsv_data = pd.read_csv(input_file, sep='\t', header=None)
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        return

    # Compute the sum of the fourth column
    try:
        sum_fourth_column = tsv_data.iloc[:, 3].sum()
    except IndexError:
        print("Error: The TSV file does not have a fourth column.")
        return

    # Count files larger than 200 bytes
    count_dir1 = count_large_files(dir1)
    count_dir2 = count_large_files(dir2)
    
    # Compute the product of the counts
    product_counts = count_dir1 * count_dir2

    # Compute the mean score
    if product_counts > 0:  # Avoid division by zero
        score = sum_fourth_column / product_counts
        score = round(score, 4)
    else:
        print("Product of file counts is zero, cannot divide.")
        return

    # Calculate the probability based on the fitted score distributions for TPs and FPs based on SCOPe40 (https://github.com/steineggerlab/foldseek/src/commons/CalcProbTP.h)
    d_tp = (0.8279 * gamma.pdf(score, a=1.8123, scale=46.0042) +
            0.1721 * gamma.pdf(score, a=1.0057, scale=563.5014)) * 0.1023
    d_fp = (0.34 * gamma.pdf(score, a=4.9259, scale=4.745) +
            0.66 * gamma.pdf(score, a=9.4834, scale=1.3136)) * 0.8977
    prob = 1 / (1 + d_fp / d_tp)
    prob = prob * 100
    prob = round(prob, 4)

    # Write the score to the output file
    try:
        result_data = pd.DataFrame({
            "dir1": [dir1_name],
            "dir2": [dir2_name],
            "score": [score],
            "prob": [prob],
        })
        result_data.to_csv(output_file, sep='\t', index=False, header=False)
    except Exception as e:
        print(f"Error writing to output file: {e}")

# Entry point
if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Compute result from a TSV file and directory information.")
    parser.add_argument("input_file", help="Path to the input TSV file.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    parser.add_argument("pdb_dir", help="Path to the PDB directory.")

    args = parser.parse_args()

    # Call the function with arguments
    compute_result(args.input_file, args.output_file, args.pdb_dir)
