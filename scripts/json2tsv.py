#!/usr/bin/env python

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
    with gzip.open(args.json, 'r') as f:
        df = json.load(f)
        if "errors" in df.keys():
            for err in df['errors']:
                print(err['title'])
        if "meta" in df.keys():
            if df['meta']['subcommand'] == "search":
                if "data" in df.keys():
                    if args.mode == "primer":
                        print("Primer", end="")
                        for k in primerKeys:
                            print("\t", k, sep="", end="")
                        print()
                        if "primers" in df['data'].keys():
                            for hit in df['data']['primers']:
                                print("Primer", end="")
                                for k in primerKeys:
                                    print("\t", hit[k], sep="", end="")
                                print()
                    else:
                        print("Amplicon", end="")
                        for k in ampliconKeys:
                            print("\t", k, sep="", end="")
                        print()
                        if "amplicons" in df['data'].keys():
                            for hit in df['data']['amplicons']:
                                print("Amplicon", end="")
                                for k in ampliconKeys:
                                    print("\t", hit[k], sep="", end="")
                                print()
