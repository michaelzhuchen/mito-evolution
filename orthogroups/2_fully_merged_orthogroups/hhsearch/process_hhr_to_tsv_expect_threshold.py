#!/bin/python

import re
import argparse

def extract_low_evalue_lines(input_file, output_file):
    prob_threshold = 30
    blank_line_count = 0
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Check if the line is blank
            if line.strip() == '':
                blank_line_count += 1
                # Terminate after the second blank line
                if blank_line_count == 2:
                    break
                continue
            
            # Ensure there's a space before '(' if it's missing
            line = re.sub(r'(\S)\(', r'\1 (', line)

            # Split the line into columns
            columns = line.split()
            if len(columns) >= 4:
                try:
                    ## For hmm+foldseek: Use probability threshold
                    prob = float(columns[2])
                    # Check if E-value is below the threshold
                    if prob >= prob_threshold:
                        outfile.write('\t'.join(columns) + '\n')
                        
                except ValueError:
                    # Continue if conversion fails (e.g., header lines)
                    continue

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract lines with E-value below a threshold.")
    parser.add_argument("input_file", help="Input file containing data to process.")
    parser.add_argument("output_file", help="Output file to save the filtered results.")
    args = parser.parse_args()

    # Call the function with the provided arguments
    extract_low_evalue_lines(args.input_file, args.output_file)

if __name__ == "__main__":
    main()