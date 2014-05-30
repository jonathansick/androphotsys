#!/usr/bin/env python
# encoding: utf-8
"""
Photometric conversions from HST bandpasses to the Cousins system.

See: Sirianni 2005, WFC3 handbook and the FSPS colour transformation modelling
scripts in the scripts/ directory.
"""

from pkg_resources import resource_stream, resource_exists

import numpy as np
from astropy.table import Table, vstack


WFC3_ZP = {"STMAG": {'f110w': 28.4401,
                     'f160w': 28.1875,
                     'f275w': 22.6034,
                     'f336w': 23.6044},
           "VEGAMAG": {'f110w': 26.0628,
                       'f160w': 24.6949,
                       'f275w': 22.6322,
                       'f336w': 23.4836},
           "ABMAG": {'f110w': 26.8223,
                     'f160w': 25.9463,
                     'f275w': 24.1305,
                     'f336w': 24.6682},
           "OBMAG": {'f110w': 0.,
                     'f160w': 0.,
                     'f275w': 0.,
                     'f336w': 0.}}

ACS_ZP = {"STMAG": {'f475w': 25.757,
                    'f606w': 26.655,
                    'f814w': 26.776},
          "VEGAMAG": {'f475w': 26.168,
                      'f606w': 26.398,
                      'f814w': 25.501},
          "ABMAG": {'f475w': 25.673,
                    'f606w': 26.486,
                    'f814w': 25.937},
          "OBMAG": {'f475w': 0.,
                    'f606w': 0.,
                    'f814w': 0.}}


def transform_acs_zp(mags, band, zp_in, zp_out):
    """Convert HST photometry from one ZP system to another.
    
    Uses values from Sirianni 2005 Table 10.
    """
    zp_systems = ('ABMAG', "STMAG", "VEGAMAG", 'OBMAG')
    assert band in ('f475w', 'f606w', 'f814w')
    assert zp_in in zp_systems
    assert zp_out in zp_systems
    return mags - ACS_ZP[zp_in][band] + ACS_ZP[zp_out][band]


def transform_wfc3_zp(mags, band, zp_in, zp_out):
    """Convert WFC3 photometry from one ZP system to another.
    
    Uses values from the WFC3 instrument handbook
    (http://www.stsci.edu/hst/wfc3/phot_zp_lbn).
    """
    zp_systems = ('ABMAG', "STMAG", "VEGAMAG", "OBMAG")
    assert band in ('f110w', 'f160w', 'f275w', 'f336w')
    assert zp_in in zp_systems
    assert zp_out in zp_systems
    return mags - WFC3_ZP[zp_in][band] + WFC3_ZP[zp_out][band]


def WFC3_275_336_to_u(m275_input, e275, m336_input, e336, zp="VEGAMAG"):
    """Convert WFC3 photometry to u-band using the FSPS modelling
    in scripts/model_wfc3_u.py
    """
    m275 = transform_wfc3_zp(m275_input, 'f275w', zp, 'ABMAG')
    m336 = transform_wfc3_zp(m336_input, 'f336w', zp, 'ABMAG')
    x = [-0.08689923, -0.35246386, 0.14181587]
    x = query_colour_transformation("u", "F336W", ("F275W", "F336W"))
    C = m275 - m336
    return m336 + x[0] + x[1] * C + x[2] * C ** 2


def WFC3_110_160_to_JK(m110_input, e110, m160_input, e160, zp="VEGAMAG"):
    """Convert WFC3 photometry to JK-band using the FSPS modelling in
    scripts/model_wfc3_wircam.py
    """
    m110 = transform_wfc3_zp(m110_input, 'f110w', zp, 'ABMAG')
    m160 = transform_wfc3_zp(m160_input, 'f160w', zp, 'ABMAG')
    jx = query_colour_transformation("J", "F160W", ("F110W", "F160W"))
    kx = query_colour_transformation("Ks", "F160W", ("F110W", "F160W"))
    C = m110 - m160
    J = m160 + jx[0] + jx[1] * C + jx[2] * C ** 2
    K = m160 + kx[0] + kx[1] * C + kx[2] * C ** 2
    return J, K


def ACS_475_814_to_gri(m475_input, e475, m814_input, e814, zp="VEGAMAG"):
    """Convert PHAT ACS photometry to gri-band using the FSPS modelling in
    scripts/model_phat_acs_gri.py
    """
    m475 = transform_acs_zp(m475_input, 'f475w', zp, 'ABMAG')
    m814 = transform_acs_zp(m814_input, 'f814w', zp, 'ABMAG')
    gx = query_colour_transformation("g", "F475W", ("F475W", "F814W"))
    rx = query_colour_transformation("r", "F814W", ("F475W", "F814W"))
    ix = query_colour_transformation("i", "F814W", ("F475W", "F814W"))
    C = m475 - m814
    g = m475 + gx[0] + gx[1] * C + gx[2] * C ** 2
    r = m814 + rx[0] + rx[1] * C + rx[2] * C ** 2
    i = m814 + ix[0] + ix[1] * C + ix[2] * C ** 2
    return g, r, i


def ACS_606_814_to_gri(m606_input, e606, m814_input, e814, zp="STMAG"):
    """Convert Brown+ ACS photometry to gri-band using the FSPS modelling in
    scripts/model_acs_brown_gri.py
    """
    m606 = transform_acs_zp(m606_input, 'f606w', zp, 'ABMAG')
    m814 = transform_acs_zp(m814_input, 'f814w', zp, 'ABMAG')
    gx = query_colour_transformation("g", "F606W", ("F606W", "F814W"))
    rx = query_colour_transformation("r", "F814W", ("F606W", "F814W"))
    ix = query_colour_transformation("i", "F814W", ("F606W", "F814W"))
    C = m606 - m814
    g = m606 + gx[0] + gx[1] * C + gx[2] * C ** 2
    r = m814 + rx[0] + rx[1] * C + rx[2] * C ** 2
    i = m814 + ix[0] + ix[1] * C + ix[2] * C ** 2
    return g, r, i


def query_colour_transformation(tmag, smag, colour):
    """Get the colour transformation for the target bandpass `tmag`
    given a source bandpass `smag` and colour `colour`. Colour transformation
    coefficients are s.t.

    tmag = smag + x[0] + x[1] * colour + x[2] * colour^2
    """
    tbl = load_colour_transformation_table()
    i = np.where((tbl['TMAG'] == tmag) &
                 (tbl['SMAG'] == smag) &
                 (tbl['COL'] == "-".join(colour)))[0][0]
    return (tbl['x0'][i], tbl['x1'][i], tbl['x2'][i])


CT_TABLE = None


def load_colour_transformation_table():
    global CT_TABLE
    if CT_TABLE is None:
        table_names = ['wfc3_u', 'phat_acs', 'acs_brown', 'wfc3_wircam']
        tables = [Table.read(_read_colour_transformation_table(n),
                             format='ascii.commented_header')
                  for n in table_names]
        CT_TABLE = vstack(tables)
    return CT_TABLE


def _read_colour_transformation_table(name):
    path = "data/{0}.txt".format(name)
    assert resource_exists(__name__, path)
    return resource_stream(__name__, path)
