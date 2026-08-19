#!/usr/bin/env python3

from __future__ import print_function
import fileinput
import json


for line in fileinput.input():
    line = line.strip()
    if not line:
        continue
    df = json.loads(line)
    if "errors" in df:
        for err in df['errors']:
            print(err.get('title', ''))
        print()
    meta = df.get('meta')
    if meta:
        subcommand = meta.get('subcommand')
        # dicey hunt
        if subcommand == "hunt":
            for hit in df.get('data', []):
                print(">", hit.get('chr'), ":", hit.get('start'), "-", hit.get('end'), " (Strand: ", hit.get('strand'), ", Distance: ", hit.get('distance'), ")", sep="")
                print(hit.get('queryalign', ''))
                print(hit.get('refalign', ''))
                print()

        # dicey search
        if subcommand == "search":
            data = df.get('data', {})
            for hit in data.get('primers', []):
                print("Primer_", hit.get("Id"), "_", "Tm", "=", hit.get("Tm"), sep="")
                print("Primer_", hit.get("Id"), "_", "Pos", "=", hit.get("Chrom"), ":", hit.get("Pos"), "-", hit.get("End"), sep="")
                print("Primer_", hit.get("Id"), "_", "Ori", "=", hit.get("Ori"), sep="")
                print("Primer_", hit.get("Id"), "_", "Name", "=", hit.get("Name"), sep="")
                print("Primer_", hit.get("Id"), "_", "MatchTm", "=", hit.get("MatchTm"), sep="")
                print("Primer_", hit.get("Id"), "_", "Seq", "=", hit.get("Seq"), sep="")
                print("Primer_", hit.get("Id"), "_", "Genome", "=", hit.get("Genome"), sep="")
            for hit in data.get('amplicons', []):
                print("Amplicon_", hit.get("Id"), "_", "Length", "=", hit.get("Length"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "Penalty", "=", hit.get("Penalty"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "ForPos", "=", hit.get("Chrom"), ":", hit.get("ForPos"), "-", hit.get("ForEnd"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "ForTm", "=", hit.get("ForTm"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "ForName", "=", hit.get("ForName"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "ForSeq", "=", hit.get("ForSeq"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "RevPos", "=", hit.get("Chrom"), ":", hit.get("RevPos"), "-", hit.get("RevEnd"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "RevTm", "=", hit.get("RevTm"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "RevName", "=", hit.get("RevName"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "RevSeq", "=", hit.get("RevSeq"), sep="")
                print("Amplicon_", hit.get("Id"), "_", "Seq", "=", hit.get("Seq"), sep="")
