#!/usr/bin/env python

from __future__ import print_function
import fileinput
import json


for line in fileinput.input():
    df = json.loads(line)
    for hit in df['data']:
        print(">", hit['chr'], ":", hit['start'], "-", hit['end'], " - Score:", hit['score'], sep="")
        {u'queryalign': u'-ACTGATGTTAGTAATG--', u'end': 18800, u'refalign': u'TAGTGATGTTAGTAATGCA', u'start': 18785, u'chr': u'KN150570.1', u'score': 14, u'strand': u'-'}
        print(hit['queryalign'])
        print(hit['refalign'])
        print()
