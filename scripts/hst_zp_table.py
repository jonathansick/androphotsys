#!/usr/bin/env python
# encoding: utf-8
"""
Write a latex table of HST zeropoints.
"""

from androphotsys.hst import WFC3_ZP, ACS_ZP
from astropy.table import Table


def main():
    ins_col = []
    band_col = []
    vega_col = []
    st_col = []
    ab_col = []

    wfc_bands = ['f275w', 'f336w', 'f110w', 'f160w']
    acs_bands = ['f475w', 'f606w', 'f814w']

    for band in wfc_bands:
        ins_col.append("WFC3")
        band_col.append(band.upper())
        vega_col.append(WFC3_ZP['VEGAMAG'][band])
        st_col.append(WFC3_ZP['STMAG'][band])
        ab_col.append(WFC3_ZP['ABMAG'][band])

    for band in acs_bands:
        ins_col.append("ACS")
        band_col.append(band.upper())
        vega_col.append(ACS_ZP['VEGAMAG'][band])
        st_col.append(ACS_ZP['STMAG'][band])
        ab_col.append(ACS_ZP['ABMAG'][band])

    table = Table(
        [ins_col, band_col, vega_col, st_col, ab_col],
        names=['Instrument', 'Filter', 'VEGAMAG', 'STMAG', 'ABMAG'])
    print table
    table.write("hst_zp.tex", format="ascii.latex")


if __name__ == '__main__':
    main()
