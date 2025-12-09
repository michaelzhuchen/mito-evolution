#!/usr/bin/env python3
"""
Script: generate_orthogroups_and_clusters.py

Description:
    Efficiently generates OrthoFinder Orthogroups.txt and
    clusters_OrthoFinder_<runID>.txt_id_pairs.txt in MCL matrix format
    (dimensions <num_sequences>x<num_ogs>), using original FASTA headers.

Usage:
    python generate_orthogroups_and_clusters.py \
        --og_fasta /path/to/OG_fasta \
        --sequence_ids /path/to/SequenceIDs.txt \
        --run_id I1.5 \
        --output_dir /path/to/WorkingDirectory

Arguments:
    --og_fasta:     Directory containing one FASTA per orthogroup
    --sequence_ids: Path to SequenceIDs.txt mapping compositeID: header
    --run_id:       Identifier for naming the clusters output
    --output_dir:   Directory to write Orthogroups.txt and clusters file

This version reads SequenceIDs once into a header->compID dict, streams
through each OG FASTA exactly once, and writes output lines directly,
avoiding intermediate storage and O(n^2) lookups.
"""
import os
import sys
import argparse
from pathlib import Path


def load_header_map(seqids_path):
    """Load SequenceIDs.txt into a dict header->compID."""
    header_to_comp = {}
    with open(seqids_path) as f:
        for line in f:
            line = line.strip()
            if not line or ':' not in line:
                continue
            comp, hdr = line.split(':', 1)
            header_to_comp[hdr.strip()] = comp.strip()
    return header_to_comp


def parse_fasta_headers(fasta_path):
    """Yield each header (without '>') from a FASTA file."""
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                yield line[1:].strip()


def main():
    parser = argparse.ArgumentParser(
        description='Generate Orthogroups and clusters efficiently')
    parser.add_argument('--og_fasta', required=True,
                        help='Directory with per-OG FASTA files')
    parser.add_argument('--sequence_ids', required=True,
                        help='Path to SequenceIDs.txt')
    parser.add_argument('--run_id', required=True,
                        help='Run identifier for clusters file name')
    parser.add_argument('--output_dir', required=True,
                        help='Directory to write output files')
    args = parser.parse_args()

    og_dir = Path(args.og_fasta)
    seqids_path = Path(args.sequence_ids)
    out_dir = Path(args.output_dir)
    run_id = args.run_id

    # Validate inputs
    if not og_dir.is_dir():
        sys.exit(f"Error: OG_fasta dir '{og_dir}' not found")
    if not seqids_path.is_file():
        sys.exit(f"Error: SequenceIDs.txt '{seqids_path}' not found")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load mapping and compute number of sequences
    header_to_comp = load_header_map(seqids_path)
    num_sequences = len(header_to_comp)

    # Collect OG FASTA files
    og_files = sorted(og_dir.glob('*.fa*'))
    num_ogs = len(og_files)
    if num_ogs == 0:
        sys.exit(f"Error: no FASTA files found in {og_dir}")

    orthogroups_path = out_dir / 'Orthogroups.txt'
    clusters_path = out_dir / f"clusters_OrthoFinder_{run_id}.txt_id_pairs.txt"

    # Open outputs
    with orthogroups_path.open('w') as ogf, clusters_path.open('w') as cpf:
        # Write MCL header
        cpf.write('(mclheader\n')
        cpf.write('mcltype matrix\n')
        cpf.write(f'dimensions {num_sequences}x{num_ogs}\n')
        cpf.write(')\n')
        cpf.write('(mclmatrix\n')
        cpf.write('begin\n')

        # Process each OG
        for idx, fasta_file in enumerate(og_files):
            og_name = fasta_file.stem

            # Read headers and map to comp IDs
            headers = list(parse_fasta_headers(fasta_file))
            comps = []
            for hdr in headers:
                comp = header_to_comp.get(hdr)
                if comp is None:
                    sys.exit(f"Error: header '{hdr}' not in SequenceIDs.txt")
                comps.append(comp)

            # Write Orthogroups line
            ogf.write(f"{og_name}: {' '.join(headers)}\n")

            # Write clusters row
            cpf.write(f"{idx} {' '.join(comps)} $\n")

        # Close matrix
        cpf.write(')\n')

    print(f"Generated Orthogroups.txt at {orthogroups_path}")
    print(f"Generated clusters matrix at {clusters_path}")

if __name__ == '__main__':
    main()
