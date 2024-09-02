#!/usr/bin/env python3
 
import json
import argparse

parser = argparse.ArgumentParser(description="Evaluate run quality from fastp output")
parser.add_argument('--json','-j',help='Fastp json output',required=True)

argv = parser.parse_args()

fastp_out = argv.json
data = json.load(open(fastp_out))

q20_bases = float(data['read1_before_filtering']['q20_bases'])
total_bases = float(data['read1_before_filtering']['total_bases'])
q20_percentage = q20_bases/total_bases*100

quality = "low"
if q20_percentage >= 80:
    quality = "high"

print(quality)