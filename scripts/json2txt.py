#!/usr/bin/env python

from __future__ import print_function
import fileinput
import json


for line in fileinput.input():
    df = json.loads(line)
    for hit in df['data']:
        print(">", hit['chr'], ":", hit['start'], "-", hit['end'], " - Score:", hit['score'], sep="")
        print(hit['queryalign'])
        print(hit['refalign'])
        print()
