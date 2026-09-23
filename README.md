# Protein K-mer Spectrum Comparator

A Python CLI tool to compare protein sequences based on k-mer spectrum similarity.

## Features
- Parses standard FASTA files.
- Computes k-mer frequency profiles for query sequences.
- Ranks target sequences based on shared k-mer overlap ratios.
- Outputs tab-separated (TSV) result files.

## Usage

This tool relies solely on the Python Standard Library (no external dependencies required).

## Command-line Arguments

    -q, --query : Path to the query FASTA file (required).

    -t, --targets : Path to target FASTA file or directory containing FASTA files (required).

    -k : K-mer length parameter (default: 6).

    -o, --output-dir : Directory to store output TSV files (default: current directory .).

## Output Files

    query_kmer_counts.tsv : K-mer frequency table for the query (sorted by occurrence).

    similarity_results.tsv : Similarity metrics for target sequences (sorted by shared k-mers).

```bash
python kmer_comparator.py -q example_data/query.fasta -t example_data/targets/ -k 6 -o results/
