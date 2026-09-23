#!/usr/bin/env python3
"""
Protein K-mer Spectrum Comparator
----------------------------------
A tool for analyzing protein sequence similarity using k-mer spectrum comparisons.
"""

import argparse
import os
from typing import Dict, List, Set, Tuple


def read_fasta(file_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file and return a dictionary mapping sequence IDs to sequences.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        Dict[str, str]: Dictionary mapping sequence IDs to uppercase sequence strings.
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(file_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq).upper()
                # 1. Extract first word after '>' as the sequence identifier
                raw_id = line[1:].split()[0]
                # 2. If it's a UniProt format (ex: sp|Q62559|IFT52_MOUSE), extract the third part
                if "|" in raw_id:
                    parts = raw_id.split("|")
                    if len(parts) >= 3:
                        current_id = parts[2]
                    else:
                        current_id = raw_id
                else:
                    current_id = raw_id
                
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = "".join(current_seq).upper()

    return sequences


def count_kmers(sequence: str, k: int) -> Dict[str, int]:
    """
    Extract k-mers from a sequence and compute their frequencies.

    Args:
        sequence: Amino acid sequence.
        k: Length of the k-mers.

    Returns:
        Dict[str, int]: Dictionary mapping k-mer strings to their counts.
    """
    kmer_counts = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts


def get_unique_kmers(kmer_counts: Dict[str, int]) -> Set[str]:
    """Return the set of k-mers that appear exactly once in a sequence."""
    return {kmer for kmer, count in kmer_counts.items() if count == 1}


def compare_kmer_sets(
    query_kmers: Set[str], target_kmers: Set[str]
) -> Tuple[int, float]:
    """
    Compare two sets of k-mers and calculate shared count and proportion.

    Args:
        query_kmers: K-mer set from the query sequence.
        target_kmers: K-mer set from the target sequence.

    Returns:
        Tuple[int, float]: Number of shared k-mers and overlap proportion.
    """
    shared = query_kmers.intersection(target_kmers)
    shared_count = len(shared)
    min_kmer_count = min(len(query_kmers), len(target_kmers))

    ratio = shared_count / min_kmer_count if min_kmer_count > 0 else 0.0
    return shared_count, ratio


def write_query_kmer_counts(kmer_counts: Dict[str, int], output_path: str) -> None:
    """Export query k-mer frequencies to a TSV file sorted by occurrence."""
    sorted_counts = sorted(kmer_counts.items(), key=lambda item: item[1], reverse=True)
    with open(output_path, "w") as out_file:
        out_file.write("kmer\toccurrence\n")
        for kmer, count in sorted_counts:
            out_file.write(f"{kmer}\t{count}\n")


def write_similarity_results(results: List[Tuple], output_path: str) -> None:
    """Export target similarity metrics to a TSV file."""
    headers = [
        "target_id",
        "shared_kmers",
        "shared_kmers_ratio",
        "shared_unique_kmers",
        "shared_unique_kmers_ratio",
    ]
    with open(output_path, "w") as out_file:
        out_file.write("\t".join(headers) + "\n")
        for row in results:
            out_file.write(f"{row[0]}\t{row[1]}\t{row[2]:.4f}\t{row[3]}\t{row[4]:.4f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Compare protein sequences based on k-mer spectrum similarity."
    )
    parser.add_argument(
        "-q", "--query", required=True, help="Path to the query FASTA file"
    )
    parser.add_argument(
        "-t",
        "--targets",
        required=True,
        help="Path to target FASTA file or directory containing FASTA files",
    )
    parser.add_argument(
        "-k", type=int, default=6, help="K-mer length parameter (default: 6)"
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Output directory for generated files (default: current directory)",
    )

    args = parser.parse_args()

    # Load query sequence
    query_dict = read_fasta(args.query)
    if not query_dict:
        raise ValueError(f"No valid FASTA sequence found in: {args.query}")
    query_id, query_seq = list(query_dict.items())[0]

    # Load target sequences (from a single file or a directory)
    target_sequences = {}
    if os.path.isdir(args.targets):
        for file_name in os.listdir(args.targets):
            if file_name.endswith((".fasta", ".fa", ".faa")):
                file_path = os.path.join(args.targets, file_name)
                target_sequences.update(read_fasta(file_path))
    elif os.path.isfile(args.targets):
        target_sequences = read_fasta(args.targets)
    else:
        raise FileNotFoundError(f"Target path not found: {args.targets}")

    # Exclude query from targets if present in the target pool
    target_sequences.pop(query_id, None)

    # Calculate query k-mers
    query_counts = count_kmers(query_seq, args.k)
    query_all_kmers = set(query_counts.keys())
    query_unique_kmers = get_unique_kmers(query_counts)

    # Export query k-mer profile
    query_out_path = os.path.join(args.output_dir, "query_kmer_counts.tsv")
    write_query_kmer_counts(query_counts, query_out_path)

    # Process target sequence comparisons
    results = []
    best_target_id = None
    max_shared_kmers = -1

    for target_id, target_seq in target_sequences.items():
        target_counts = count_kmers(target_seq, args.k)
        target_all_kmers = set(target_counts.keys())
        target_unique_kmers = get_unique_kmers(target_counts)

        # Compute metrics
        shared_kmers, prop_kmers = compare_kmer_sets(query_all_kmers, target_all_kmers)
        shared_uniq, prop_uniq = compare_kmer_sets(query_unique_kmers, target_unique_kmers)

        results.append((target_id, shared_kmers, prop_kmers, shared_uniq, prop_uniq))

        if shared_kmers > max_shared_kmers:
            max_shared_kmers = shared_kmers
            best_target_id = target_id

    # Sort results by shared k-mers count (descending)
    results.sort(key=lambda x: x[1], reverse=True)

    # Export results table
    results_out_path = os.path.join(args.output_dir, "similarity_results.tsv")
    write_similarity_results(results, results_out_path)

    # Terminal summary output
    print(
        f"Cette protéine ({query_id}) est composée de {len(query_unique_kmers)} {args.k}-mers uniques. "
        f"La séquence la plus proche possède {max_shared_kmers} {args.k}-mers en commun, "
        f"avec l’identifiant {best_target_id}."
    )


if __name__ == "__main__":
    main()