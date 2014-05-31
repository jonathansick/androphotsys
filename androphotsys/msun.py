#!/usr/bin/env python
# encoding: utf-8
"""
Solar magnitude database
"""

from pkg_resources import resource_stream, resource_exists

import numpy as np
from astropy.table import Table

SYSTEMS = ('ABMAG', 'VEGAMAG')


def solar_mag(band, zp='ABMAG'):
    """Get the solar absolute magnitude for the given bandpass.

    Zeropoint is ABMAG by default, can be switched to VEGAMAG.
    """
    global SYSTEMS
    assert zp in SYSTEMS
    msun_tbl = _load_msun_table()
    i = np.where(msun_tbl['name'] == band.lower())[0][0]
    return msun_tbl[i][zp]


MSUN_TABLE = None


def _load_msun_table():
    global MSUN_TABLE
    if MSUN_TABLE is None:
        path = "data/msun.txt"
        assert resource_exists(__name__, path)
        f = resource_stream(__name__, path)
        MSUN_TABLE = Table.read(f, format='ascii.commented_header')
    return MSUN_TABLE
