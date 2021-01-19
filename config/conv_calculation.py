
import numpy as np
import sys
from copy import deepcopy
from pprint import pprint

cells = ['pyramidalcell', 'axoaxoniccell', 'bistratifiedcell', 'pvbasketcell', 'cckcell', 
        'ngfcell', 'ivycell', 'olmcell', 'scacell', 'ca3cell', 'eccell']
cell_numbers = [311500, 1470, 2210, 5530, 3600, 3580, 8810, 1640, 400, 204700, 250000]

info = {}


def read_and_parse(filename):
    f = open(filename, 'r')
    lcount = 0
    for line in f.readlines():
        if lcount != 0:
            line = line.strip('\n')
            pre, post, wgt, ncons, nsyns = line.split(' ')
            
            if pre not in info: 
                info[pre] = {}
            if post not in info:
                info[post] = {}
            val = float(ncons) / cell_numbers[cells.index(post)] * float(nsyns)
            info[post][pre] = val
        lcount += 1
    f.close()

    total_incoming = {c: 0 for c in cells}
    fraction = deepcopy(info)
    for cell in cells:
        for pre in list(info[cell].keys()):
            total_incoming[cell] += info[cell][pre]
        for pre in list(info[cell].keys()):
            fraction[cell][pre] /= float(total_incoming[cell])
    pprint(total_incoming)
    pprint(fraction)


filename = sys.argv[1]
read_and_parse(filename)

