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
import re
from Bio import SeqIO
import gzip

STATS_MEGAHIT_FASTA_HEADER_REGEX = r"k\d+_\d+ flag=(\d) multi=([\d.]+) len=(\d+)"

def open_fasta_file(filename):
    if filename.endswith('.gz'):
        f = gzip.open(filename, "rt")
    else:
        f = open(filename, "rt")
    return f


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Rename megahit assembly headers to spades style")
    parser.add_argument('-i', '--input', required=True, help="Megahit contigs file")
    parser.add_argument('-o', '--output', type=str, required=False, help="Name of output file",
                        default="renamed_megahit.fasta")
    return parser.parse_args()


def rename_megahit_headers(input_file, output_name):
    counter = 1
    with open(output_name, 'w') as file_out:
        handle = open_fasta_file(input_file)
        for record in SeqIO.parse(handle, "fasta"):
            match = re.search(STATS_MEGAHIT_FASTA_HEADER_REGEX, record.description)
            if match:
                multi_value = match.group(2)
                length_value = match.group(3)
                line = f"NODE_{counter}_length_{length_value}_cov_{multi_value}\n"
                counter += 1
            else:
                line = f">NODE_{counter}"
            record.id = line
            record.description = ""
            SeqIO.write(record, file_out, "fasta")
        handle.close()


def main():
    args = parse_args()
    rename_megahit_headers(args.input, args.output)


if __name__ == "__main__":
    main()
