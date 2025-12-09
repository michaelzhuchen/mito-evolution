#!/usr/bin/env python
"""
Module Docstring
Remove sequences from the trimmed alignment that are only composed of gaps.
"""

def read_fasta(file_path):
    """
    Read a FASTA file and return a dictionary mapping sequence IDs to sequences.
    """
    sequences = {}
    current_sequence_id = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence_id = line[1:]
                sequences[current_sequence_id] = ""
            else:
                sequences[current_sequence_id] += line
    return sequences

def write_fasta(filtered_sequences, output_file_path):
    """
    Write filtered sequences to a new FASTA file.
    """
    with open(output_file_path, 'w') as file:
        for sequence_id, sequence in filtered_sequences.items():
            file.write(f'>{sequence_id}\n{sequence}\n')

def remove_dash_sequences(sequences):
    """
    Remove sequences composed only of "-".
    """
    filtered_sequences = {sequence_id: sequence for sequence_id, sequence in sequences.items() if set(sequence) != {'-'}}
    return filtered_sequences

def main(input_file, output_file):
    """
    Main function to remove dash-only sequences from a FASTA file.
    """
    sequences = read_fasta(input_file)
    filtered_sequences = remove_dash_sequences(sequences)
    write_fasta(filtered_sequences, output_file)

if __name__ == "__main__":
    import sys

    input_file_path = sys.argv[1]
    output_file_path = input_file_path # Overwrite the initial file

    main(input_file_path, output_file_path)

