#!/usr/bin/env python3

from __future__ import print_function
import argparse
import json
import gzip

parser = argparse.ArgumentParser(description='Convert primers/amplicons from json to tsv format')
parser.add_argument('-j', '--json', metavar='js.gz', required=True, dest='json', help='json file (required)')
parser.add_argument('-m', '--mode', metavar='primer', required=True, dest='mode', help='mode [primer|amplicon] (required)')
args = parser.parse_args()


primerKeys = ["Name", "Tm", "Chrom", "Pos", "End", "Ori", "MatchTm", "Seq", "Genome"]
ampliconKeys = ["Id", "Length", "Penalty", "Chrom", "ForPos", "ForEnd", "ForTm", "ForName", "ForSeq", "Chrom", "RevPos", "RevEnd", "RevTm", "RevName", "RevSeq", "Seq"]

if args.json:
    with gzip.open(args.json, 'rt') as f:
        df = json.load(f)
        if "errors" in df:
            for err in df['errors']:
                print(err.get('title', ''))
        meta = df.get('meta')
        if meta and meta.get('subcommand') == "search":
            data = df.get('data', {})
            if args.mode == "primer":
                print("Primer", end="")
                for k in primerKeys:
                    print("\t", k, sep="", end="")
                print()
                for hit in data.get('primers', []):
                    print("Primer", end="")
                    for k in primerKeys:
                        print("\t", hit.get(k, ""), sep="", end="")
                    print()
            else:
                print("Amplicon", end="")
                for k in ampliconKeys:
                    print("\t", k, sep="", end="")
                print()
                for hit in data.get('amplicons', []):
                    print("Amplicon", end="")
                    for k in ampliconKeys:
                        print("\t", hit.get(k, ""), sep="", end="")
                    print()
