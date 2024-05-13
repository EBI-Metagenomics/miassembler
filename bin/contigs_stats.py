#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018-2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import os
from Bio import SeqIO
import gzip
import json

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Get contigs stats")
    parser.add_argument('-f', '--fasta', required=True, help="Contigs file")
    parser.add_argument('-c', '--coverage', required=True, help="Coverage file")
    parser.add_argument('-p', '--fastp-json', required=True, help="Fastp json report")
    parser.add_argument('-o', '--output', type=str, required=False, help="Name of output file",
                        default="assembly_stats.json")
    return parser.parse_args()


def open_file(filename):
    if filename.endswith('.gz'):
        f = gzip.open(filename, "rt")
    else:
        f = open(filename, "rt")
    return f


def get_filtered_stats(contigs, limit):
    contigs = list(filter(lambda x: x >= limit, contigs))
    return len(contigs), sum(contigs)


def get_contigs_lenghts(contigs):
    contigs_lengths = []
    handle = open_file(contigs)
    for record in SeqIO.parse(handle, "fasta"):
        contigs_lengths.append(len(record.seq))
    handle.close()
    return contigs_lengths


def get_n50_l50(contigs_unsorted):
    contigs = sorted(contigs_unsorted, reverse=True)
    if not len(contigs):
        return 0, 0
    half_contigs = sum(contigs) / 2
    l50 = 0
    total_top = contigs[l50]
    while total_top < half_contigs:
        l50 += 1
        total_top += contigs[l50]
    return contigs[l50], l50 + 1


def calc_assembled_pairs(coverage_file):
    assembled_pairs, assembly_length = 0, 0
    if os.stat(coverage_file).st_size:
        f = open_file(coverage_file)
        next(f)
        for line in f:
            line = line.split()
            length = float(line[1])
            read_depth = float(line[2])
            assembled_pairs += length * read_depth
            assembly_length += length
        f.close()
    return assembled_pairs, assembly_length


def parse_fastp_report(fastp):
    with open(fastp, 'r', encoding='utf-8') as f:
        data = json.load(f)
        return data["summary"]["after_filtering"]["total_reads"], data["summary"]["after_filtering"]["total_bases"]

def main():
    args = parse_args()
    assembly_stats = {}
    contigs_lengths = get_contigs_lenghts(args.fasta)
    n50, l50 = get_n50_l50(contigs_lengths)
    assembled_pairs, assembly_length = calc_assembled_pairs(args.coverage)
    read_count, input_base_count = parse_fastp_report(args.fastp_json)
    coverage = assembled_pairs / input_base_count
    if not sum(contigs_lengths):
        coverage_depth = float(0)
    else:
        coverage_depth = assembled_pairs / assembly_length

    assembly_stats["input_read_count"] = read_count
    assembly_stats["limited_1000"] = get_filtered_stats(contigs_lengths, 1000)
    assembly_stats["limited_10000"] = get_filtered_stats(contigs_lengths, 10000)
    assembly_stats["limited_50000"] = get_filtered_stats(contigs_lengths, 50000)
    assembly_stats["num_contigs"] = len(contigs_lengths)
    assembly_stats["assembly_length"] = sum(contigs_lengths)
    assembly_stats["largest_contig"] = max(contigs_lengths)
    assembly_stats["n50"] = n50
    assembly_stats["l50"] = l50
    assembly_stats["coverage"] = coverage
    assembly_stats["coverage_depth"] = coverage_depth

    with open(args.output, 'w', encoding='utf-8') as file_out:
        json.dump(assembly_stats, file_out, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    main()
