#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018-2023 EMBL - European Bioinformatics Institute
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

import sys
import re
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-1',
                    required=True,
                    help='first bin folder name')

version_re = re.compile("SPAdes version: (.+)")
mode_re = re.compile(r"\s+(.+) mode")

megahit_version_re = re.compile(r"MEGAHIT v(\d+\.\d+\.\d+)")

def get_assembler_and_version(log_file):
    confirmed = None
    version = None
    mode = None
    with open(log_file, "r") as f:
        line = next(f, None)
        while line:
            if "System information" in line:
                match = version_re.findall(next(f))
                if len(match) > 0:
                    version = match[0]
            if "Dataset parameters:" in line:
                match = mode_re.findall(next(f))
                if len(match) > 0:
                    mode = match[0]
            if "Assembling finished." in line and version and mode:
                confirmed = (mode, version)
            if megahit_version_re.search(line):
                version = megahit_version_re.findall(line)[0]
                return "megahit", version
            line = next(f, "")
    if confirmed:
        mode, version = confirmed
    if mode == "Metagenomic":
        name = "metaspades"
    else:
        name = "spades"
    return name, version

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("file", help='assembler log file')

    argv=sys.argv[1:]
    args = parser.parse_args(argv)
    
    name, version = get_assembler_and_version(args.file)
    
    print(f'{{"assembler_name": "{name}", ')
    print(f'"assembler_version": "{version}", ')
