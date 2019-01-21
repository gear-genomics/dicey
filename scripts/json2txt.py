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
    if "data" in df.keys():
        for hit in df['data']:
            print(">", hit['chr'], ":", hit['start'], "-", hit['end'], " (Strand: ", hit['strand'], ", Distance: ", hit['distance'], ")", sep="")
            print(hit['queryalign'])
            print(hit['refalign'])
            print()
