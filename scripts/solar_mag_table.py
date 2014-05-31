#!/usr/bin/env python
# encoding: utf-8
"""
Make a table of solar magnitudes in various bands by bootstrapping
(python)-FSPS.
"""

import os
import numpy as np
from astropy.table import Table, join
from fsps.fsps import FILTERS


def main():
    msun_table = read_msun_table()
    ftable = load_filters()
    tbl = join(msun_table, ftable, keys='index')
    print tbl['index', 'name', 'ABMAG']
    dirname = "androphotsys/data"
    if os.path.exists(dirname):
        out_path = os.path.join(dirname, "msun.txt")
    else:
        out_path = "msun.txt"
    tbl.write(out_path, format="ascii.commented_header")


def read_msun_table():
    path = os.path.expandvars(os.path.join("$SPS_HOME", "data", "magsun.dat"))
    # msun = Table.read(path, format="ascii.commented_header")
    dt = [('index', int), ('ABMAG', float), ('VEGAMAG', float)]
    data = np.loadtxt(path, dtype=np.dtype(dt))
    msun = Table(data)
    return msun


def load_filters():
    keys, indices, description = [], [], []
    for key, f in FILTERS.iteritems():
        keys.append(key)
        indices.append(f.index + 1)
        description.append(f.fullname)
    ftable = Table([indices, keys, description],
                   names=['index', 'name', 'fullname'])
    return ftable


if __name__ == '__main__':
    main()
