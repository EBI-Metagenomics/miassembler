#!/usr/bin/env python

import csv
from typing import Tuple
import json
import argparse
import gzip

def get_assembled_base_pairs_and_length(jgi_summarize_coverage_file_gz: str) -> Tuple[int, int]:
    """
    Computes the total assembled pairs and assembly length from a JGI summarize coverage file.

    :param jgi_summarize_coverage_file: Path to the JGI summarize coverage file.
    :type jgi_summarize_coverage_file: str

    :return: A tuple containing the total number of assembled pairs and the total assembly length.
    :rtype: Tuple[int, int]

    :raises ValueError: If 'contigLen' or 'totalAvgDepth' columns contain non-numeric values.
    """
    assembled_base_pairs = 0
    assembly_length = 0
    with gzip.open(jgi_summarize_coverage_file_gz, "rt") as file_handle:
        csv_reader = csv.DictReader(file_handle, delimiter="\t")
        for row in csv_reader:

            contig_length = float(row["contigLen"])

            if not contig_length.is_integer():
                raise ValueError(f"The column 'contigLen' has an invalid value: {contig_length}")

            if int(contig_length) == 0:
                raise ValueError("The column 'contigLen' cannot have a contig of len 0")

            total_avg_depth_str = row["totalAvgDepth"]
            # If total total_avg_depth_str is not a float, a ValueError should be raised
            total_avg_depth = float(total_avg_depth_str)

            assembled_base_pairs += contig_length * total_avg_depth
            assembly_length += contig_length

    return assembled_base_pairs, assembly_length


def get_total_bases_before_filtering(fastp_json: str) -> int:
    """
    Retrieves the total number of bases before filtering from a Fastp JSON file.

    :param fastp_json: Path to the Fastp JSON file.
    :type fastp_json: str

    :return: The total number of bases before filtering.
    :rtype: int

    :raises ValueError: If the 'total_bases' value in the Fastp JSON is not numeric.
    """
    with open(fastp_json, "r") as json_handle:
        fastp_json_data = json.load(json_handle)
        # Get the total bases before filtering (the raw bp count)
        total_bases_before_filtering = fastp_json_data.get("summary").get("before_filtering").get("total_bases")

        if not isinstance(total_bases_before_filtering, int):
            raise ValueError(
                f"The value of 'total_bases_before_filtering' in the fastp json is not numeric: {total_bases_before_filtering}"
            )

    return total_bases_before_filtering


def main():
    parser = argparse.ArgumentParser(
        description="Calculate coverage and coverage depth from the jgi_summarize_bam_contig_depths and fastp data."
    )
    parser.add_argument("-j", "--jgi", type=str, required=True, help="Path to the JGI summarize coverage compressed (gz) file.")
    parser.add_argument("-f", "--fastp", type=str, required=True, help="Path to the qc fastp json report file.")
    parser.add_argument("-o", "--outfile", type=str, default="coverage_report.json", help="Name of the output file.")

    args = parser.parse_args()

    assembled_base_pairs, assembly_length = get_assembled_base_pairs_and_length(args.jgi)
    total_bases_before_filtering = get_total_bases_before_filtering(args.fastp)

    coverage = assembled_base_pairs / total_bases_before_filtering
    coverage_depth = assembled_base_pairs / assembly_length if assembly_length > 0 else 0

    with open(args.outfile, "w") as outfile_file:
        json.dump({"coverage": coverage, "coverage_depth": coverage_depth}, outfile_file)


if __name__ == "__main__":
    main()
