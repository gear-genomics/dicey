#!/usr/bin/env python

from __future__ import print_function
import fileinput
import json


for line in fileinput.input():
    df = json.loads(line)
    if "errors" in df.keys():
        for err in df['errors']:
            print(err['title'])
        print()
    if "meta" in df.keys():
        # dicey hunt
        if df['meta']['subcommand'] == "hunt":
            if "data" in df.keys():
                for hit in df['data']:
                    print(">", hit['chr'], ":", hit['start'], "-", hit['end'], " (Strand: ", hit['strand'], ", Distance: ", hit['distance'], ")", sep="")
                    print(hit['queryalign'])
                    print(hit['refalign'])
                    print()

        # dicey search
        if df['meta']['subcommand'] == "search":
            if "data" in df.keys():
                if "primers" in df['data'].keys():
                    for hit in df['data']['primers']:
                        print("Primer_", hit["Id"], "_", "Tm", "=", hit["Tm"], sep="")
                        print("Primer_", hit["Id"], "_", "Pos", "=", hit["Chrom"], ":", hit["Pos"], "-", hit["End"], sep="")
                        print("Primer_", hit["Id"], "_", "Ori", "=", hit["Ori"], sep="")
                        print("Primer_", hit["Id"], "_", "Name", "=", hit["Name"], sep="")
                        print("Primer_", hit["Id"], "_", "MatchTm", "=", hit["MatchTm"], sep="")
                        print("Primer_", hit["Id"], "_", "Seq", "=", hit["Seq"], sep="")
                if "amplicons" in df['data'].keys():
                    for hit in df['data']['amplicons']:
                        print("Amplicon_", hit["Id"], "_", "Length", "=", hit["Length"], sep="")
                        print("Amplicon_", hit["Id"], "_", "Penalty", "=", hit["Penalty"], sep="")
                        print("Amplicon_", hit["Id"], "_", "ForPos", "=", hit["Chrom"], ":", hit["ForPos"], "-", hit["ForEnd"], sep="")
                        print("Amplicon_", hit["Id"], "_", "ForTm", "=", hit["ForTm"], sep="")
                        print("Amplicon_", hit["Id"], "_", "ForName", "=", hit["ForName"], sep="")
                        print("Amplicon_", hit["Id"], "_", "ForSeq", "=", hit["ForSeq"], sep="")
                        print("Amplicon_", hit["Id"], "_", "RevPos", "=", hit["Chrom"], ":", hit["RevPos"], "-", hit["RevEnd"], sep="")
                        print("Amplicon_", hit["Id"], "_", "RevTm", "=", hit["RevTm"], sep="")
                        print("Amplicon_", hit["Id"], "_", "RevName", "=", hit["RevName"], sep="")
                        print("Amplicon_", hit["Id"], "_", "RevSeq", "=", hit["RevSeq"], sep="")
