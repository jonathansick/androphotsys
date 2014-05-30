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
    
    Uses values from Sirianni 2005 Table 10."""
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
    # jx = [0.00547237, 0.66899301, -0.04243658]
    # kx = [0.31956416, -0.79017206, -0.50508958]
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
    # gx = [0.00011405, 0.01957309, 0.00617221]
    # rx = [0.01195939, 0.4330676, -0.07053753]
    # ix = [0.00137411, 0.10460093, -0.01623761]
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
    # gx = [-0.01747885, 0.92625851, 0.48930022]
    # rx = [0.00288843, 0.79799345, -0.02999439]
    # ix = [-0.00116191, 0.19501024, -0.00117196]
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


def ACS_606_814_to_VRI(m606_input, e606, m814_input, e814, zp="STMAG"):
    """Uses the transformations presented by Sirianni 2005 Table 22.
    Transforms ACS magnitudes from either the STMAG or VEGAMAG systems as well
    based on Sirianni 2005 Table 11.
    
    Here's the Sage notebook::

        var('m606 e606 m814 e814 V R I c0VVI c1VVI c0IVI c1IVI c0IRI c1IRI
            c0VVIe c1VVIe c0IVIe c1IVIe c0IRIe c1IRIe')
        # TMAG = SMAG + c0 + c1 * TCOL
        eq1 = V == m606 + c0VVI + c1VVI * (V - I)
        eq2 = I == m814 + c0IVI + c1IVI * (V - I)
        eq3 = I == m814 + c0IRI + c1IRI * (R - I)
        s = solve([eq1, eq2, eq3], V, R, I, solution_dict=True)[0]
        print "V = %s" % s[V]
        print "R = %s" % s[R]
        print "I = %s" % s[I]
        Ve = sqrt((s[V].diff(m606)*e606)^2
            + (s[V].diff(m814)*e814)^2
            + (s[V].diff(c0VVI)*c0VVIe)^2
            + (s[V].diff(c1VVI)*c1VVIe)^2
            + (s[V].diff(c0IVI)*c0IVIe)^2
            + (s[V].diff(c1IVI)*c1IVIe)^2
            + (s[V].diff(c0IRI)*c0IRIe)^2
            + (s[V].diff(c1IRI)*c1IRIe)^2)
        Re = sqrt((s[R].diff(m606)*e606)^2
            + (s[R].diff(m814)*e814)^2
            + (s[R].diff(c0VVI)*c0VVIe)^2
            + (s[R].diff(c1VVI)*c1VVIe)^2
            + (s[R].diff(c0IVI)*c0IVIe)^2
            + (s[R].diff(c1IVI)*c1IVIe)^2
            + (s[R].diff(c0IRI)*c0IRIe)^2
            + (s[R].diff(c1IRI)*c1IRIe)^2)
        Ie = sqrt((s[I].diff(m606)*e606)^2
            + (s[I].diff(m814)*e814)^2
            + (s[I].diff(c0VVI)*c0VVIe)^2
            + (s[I].diff(c1VVI)*c1VVIe)^2
            + (s[I].diff(c0IVI)*c0IVIe)^2
            + (s[I].diff(c1IVI)*c1IVIe)^2
            + (s[I].diff(c0IRI)*c0IRIe)^2
            + (s[I].diff(c1IRI)*c1IRIe)^2)
        print "Ve = %s" % Ve
        print "Re = %s" % Re
        print "Ie = %s" % Ie
    """
    # Coefficients from Table 22
    c0VVI = 26.325  # pm 0.057
    c0VVIe = 0.057
    c1VVI = 0.236  # pm 0.058
    c1VVIe = 0.058
    c0IVI = 25.495  # pm 0.015
    c0IVIe = 0.015
    c1IVI = -0.002  # pm 0.017
    c1IVIe = 0.017
    c0IRI = 25.492  # pm 0.013
    c0IRIe = 0.013
    c1IRI = 0.002  # pm 0.003
    c1IRIe = 0.003

    m606 = transform_acs_zp(m606_input, 'f606w', zp, 'OBMAG')
    m814 = transform_acs_zp(m814_input, 'f814w', zp, 'OBMAG')

    V = (c0VVI * (c1IVI + 1) - (c0IVI + m814) * c1VVI + (c1IVI + 1) * m606) \
            / (c1IVI - c1VVI + 1)
    R = (c0VVI * (c1IRI + 1) * c1IVI + (c1IRI + 1) * c1IVI * m606
            + c0IVI * (c1IRI + 1) - (c0IRI + m814) * c1IVI
            - (c0IVI * (c1IRI + 1) + c1IRI * m814 - c0IRI) * c1VVI
            + c1IRI * m814 - c0IRI) / (c1IRI * c1IVI - c1IRI * c1VVI + c1IRI)
    I = (c0VVI * c1IVI - (c0IVI + m814) * c1VVI + c1IVI * m606 + c0IVI
            + m814) / (c1IVI - c1VVI + 1)

    Ve = np.sqrt(c1VVIe ** 2 * ((c0IVI + m814) / (c1IVI - c1VVI + 1) - (c0VVI * (c1IVI +
        1) - (c0IVI + m814) * c1VVI + (c1IVI + 1) * m606) / (c1IVI - c1VVI + 1) ** 2) ** 2 +
        c1IVIe ** 2 * ((c0VVI + m606) / (c1IVI - c1VVI + 1) - (c0VVI * (c1IVI + 1) -
        (c0IVI + m814) * c1VVI + (c1IVI + 1) * m606) / (c1IVI - c1VVI + 1) ** 2) ** 2 +
        c0VVIe ** 2 * (c1IVI + 1) ** 2 / (c1IVI - c1VVI + 1) ** 2 + c0IVIe ** 2 * c1VVI ** 2 / (c1IVI -
        c1VVI + 1) ** 2 + (c1IVI + 1) ** 2 * e606 ** 2 / (c1IVI - c1VVI + 1) ** 2 +
        c1VVI ** 2 * e814 ** 2 / (c1IVI - c1VVI + 1) ** 2)
    Re = np.sqrt(c0VVIe ** 2 * (c1IRI + 1) ** 2 * c1IVI ** 2 / (c1IRI * c1IVI - c1IRI * c1VVI +
        c1IRI) ** 2 + (c1IRI + 1) ** 2 * c1IVI ** 2 * e606 ** 2 / (c1IRI * c1IVI - c1IRI * c1VVI +
        c1IRI) ** 2 + c1VVIe ** 2 * ((c0IVI * (c1IRI + 1) + c1IRI * m814 -
        c0IRI) / (c1IRI * c1IVI - c1IRI * c1VVI + c1IRI) - (c0VVI * (c1IRI + 1) * c1IVI +
        (c1IRI + 1) * c1IVI * m606 + c0IVI * (c1IRI + 1) - (c0IRI + m814) * c1IVI -
        (c0IVI * (c1IRI + 1) + c1IRI * m814 - c0IRI) * c1VVI + c1IRI * m814 -
        c0IRI) * c1IRI / (c1IRI * c1IVI - c1IRI * c1VVI + c1IRI) ** 2) ** 2 +
        c1IVIe ** 2 * ((c0VVI * (c1IRI + 1) + (c1IRI + 1) * m606 - c0IRI -
        m814) / (c1IRI * c1IVI - c1IRI * c1VVI + c1IRI) - (c0VVI * (c1IRI + 1) * c1IVI +
        (c1IRI + 1) * c1IVI * m606 + c0IVI * (c1IRI + 1) - (c0IRI + m814) * c1IVI -
        (c0IVI * (c1IRI + 1) + c1IRI * m814 - c0IRI) * c1VVI + c1IRI * m814 -
        c0IRI) * c1IRI / (c1IRI * c1IVI - c1IRI * c1VVI + c1IRI) ** 2) ** 2 +
        c1IRIe ** 2 * ((c0VVI * c1IVI - (c0IVI + m814) * c1VVI + c1IVI * m606 + c0IVI +
        m814) / (c1IRI * c1IVI - c1IRI * c1VVI + c1IRI) - (c0VVI * (c1IRI + 1) * c1IVI +
        (c1IRI + 1) * c1IVI * m606 + c0IVI * (c1IRI + 1) - (c0IRI + m814) * c1IVI -
        (c0IVI * (c1IRI + 1) + c1IRI * m814 - c0IRI) * c1VVI + c1IRI * m814 -
        c0IRI) * (c1IVI - c1VVI + 1) / (c1IRI * c1IVI - c1IRI * c1VVI + c1IRI) ** 2) ** 2 +
        ((c1IRI + 1) * c1VVI - c1IRI - 1) ** 2 * c0IVIe ** 2 / (c1IRI * c1IVI - c1IRI * c1VVI +
        c1IRI) ** 2 + c0IRIe ** 2 * (c1IVI - c1VVI + 1) ** 2 / (c1IRI * c1IVI - c1IRI * c1VVI +
        c1IRI) ** 2 + (c1IRI * c1VVI - c1IRI + c1IVI) ** 2 * e814 ** 2 / (c1IRI * c1IVI -
        c1IRI * c1VVI + c1IRI) ** 2)
    Ie = np.sqrt(c1VVIe ** 2 * ((c0IVI + m814) / (c1IVI - c1VVI + 1) - (c0VVI * c1IVI -
        (c0IVI + m814) * c1VVI + c1IVI * m606 + c0IVI + m814) / (c1IVI - c1VVI +
        1) ** 2) ** 2 + c1IVIe ** 2 * ((c0VVI + m606) / (c1IVI - c1VVI + 1) - (c0VVI * c1IVI -
        (c0IVI + m814) * c1VVI + c1IVI * m606 + c0IVI + m814) / (c1IVI - c1VVI +
        1) ** 2) ** 2 + c0VVIe ** 2 * c1IVI ** 2 / (c1IVI - c1VVI + 1) ** 2 + c0IVIe ** 2 * (c1VVI -
        1) ** 2 / (c1IVI - c1VVI + 1) ** 2 + c1IVI ** 2 * e606 ** 2 / (c1IVI - c1VVI + 1) ** 2 +
        (c1VVI - 1) ** 2 * e814 ** 2 / (c1IVI - c1VVI + 1) ** 2)
    return V, Ve, R, Re, I, Ie


def ACS_475_814_to_BVRI(m475_input, e475, m814_input, e814, zp="STMAG"):
    """Transform F475W, F814W to BVRI."""
    m475 = transform_acs_zp(m475_input, 'f475w', zp, 'OBMAG')
    m814 = transform_acs_zp(m814_input, 'f814w', zp, 'OBMAG')

    c0BBV = 26.146
    c0BBVe = 0.002
    c1BBV = 0.389
    c1BBVe = 0.004
    c2BBV = 0.032
    c2BBVe = 0.001

    c0BBR = 26.150
    c0BBRe = 0.004
    c1BBR = 0.291
    c1BBRe = 0.008
    c2BBR = -0.110
    c2BBRe = 0.028

    c0IVI = 25.495  # pm 0.015
    c0IVIe = 0.015
    c1IVI = -0.002  # pm 0.017
    c1IVIe = 0.017

    c0IRI = 25.492  # pm 0.013
    c0IRIe = 0.013
    c1IRI = 0.002  # pm 0.003
    c1IRIe = 0.003

    B = (c0BBV*(c1IRI + 1)*c1IVI*c2BBR + (c1IRI + 1)*c1IVI*c2BBR*m475 -
        ((c1IRI*c1IVI + c1IRI)*c0BBR - (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI
        + c1IRI*m814 - c0IRI)*c2BBR + (c1IRI*c1IVI + c1IRI)*m475)*c1BBV)/((c1IRI
        + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV)
    V = ((c1IRI*c1IVI + c1IRI)*c2BBR*m475 + (c1IRI*c1IVI + c1IRI)*c0BBR -
        (c1IRI*c1IVI - (c1IRI*c1IVI + c1IRI)*c2BBR + c1IRI)*c0BBV -
        ((c1IRI*c1IVI + c1IRI)*c0BBR - (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI
        + c1IRI*m814 - c0IRI)*c2BBR + (c1IRI*c1IVI + c1IRI)*m475)*c1BBV -
        (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 -
        c0IRI)*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV)
    R = ((c1IRI + 1)*c1IVI*c2BBR*m475 + c0BBR*(c1IRI + 1)*c1IVI + ((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI + 1)*c1IVI)*c0BBV - (c0BBR*(c1IRI + 1)*c1IVI +
        (c1IRI + 1)*c1IVI*m475 + c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI -
        (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR +
        c1IRI*m814 - c0IRI)*c1BBV)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI -
        (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV)
    I = (c1IRI*c1IVI*c2BBR*m475 + c0BBR*c1IRI*c1IVI + (c0IRI +
        m814)*c1IVI*c2BBR + (c1IRI*c1IVI*c2BBR - c1IRI*c1IVI)*c0BBV -
        (c0BBR*c1IRI*c1IVI + c1IRI*c1IVI*m475 + c0IVI*c1IRI - (c0IVI*c1IRI -
        (c0IRI + m814)*c1IVI + c1IRI*m814)*c2BBR + c1IRI*m814)*c1BBV)/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV)

    Be = np.sqrt(c0IVIe ** 2*c1BBV ** 2*(c1IRI + 1) ** 2*c2BBR ** 2/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2
        + c0IRIe ** 2*c1BBV ** 2*(c1IVI + 1) ** 2*c2BBR ** 2/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + c0BBVe ** 2*(c1IRI
        + 1) ** 2*c1IVI ** 2*c2BBR ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI
        - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + c1BBV ** 2*(c1IRI -
        c1IVI) ** 2*c2BBR ** 2*e814 ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI
        - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (c1IRI*c1IVI +
        c1IRI) ** 2*c0BBRe ** 2*c1BBV ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI -
        (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + c2BBRe ** 2*((c0BBV*(c1IRI +
        1)*c1IVI + (c1IRI + 1)*c1IVI*m475 + (c0IVI*(c1IRI + 1) - (c0IRI +
        m814)*c1IVI + c1IRI*m814 - c0IRI)*c1BBV)/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) - (c0BBV*(c1IRI +
        1)*c1IVI*c2BBR + (c1IRI + 1)*c1IVI*c2BBR*m475 - ((c1IRI*c1IVI +
        c1IRI)*c0BBR - (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 -
        c0IRI)*c2BBR + (c1IRI*c1IVI + c1IRI)*m475)*c1BBV)*(c1BBV*(c1IRI - c1IVI)
        + (c1IRI + 1)*c1IVI)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2) ** 2 + c1IVIe ** 2*((c0BBV*(c1IRI + 1)*c2BBR +
        (c1IRI + 1)*c2BBR*m475 - (c0BBR*c1IRI + (c0IRI + m814)*c2BBR +
        c1IRI*m475)*c1BBV)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) + (c0BBV*(c1IRI + 1)*c1IVI*c2BBR + (c1IRI +
        1)*c1IVI*c2BBR*m475 - ((c1IRI*c1IVI + c1IRI)*c0BBR - (c0IVI*(c1IRI + 1)
        - (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR + (c1IRI*c1IVI +
        c1IRI)*m475)*c1BBV)*(c1BBV*(c1IRI + c2BBR) - (c1IRI + 1)*c2BBR)/((c1IRI
        + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)*c1BBV) ** 2) ** 2 + c1IRIe ** 2*((c0BBV*c1IVI*c2BBR + c1IVI*c2BBR*m475 -
        (c0BBR*(c1IVI + 1) - (c0IVI + m814)*c2BBR + (c1IVI +
        1)*m475)*c1BBV)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) + (c0BBV*(c1IRI + 1)*c1IVI*c2BBR + (c1IRI +
        1)*c1IVI*c2BBR*m475 - ((c1IRI*c1IVI + c1IRI)*c0BBR - (c0IVI*(c1IRI + 1)
        - (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR + (c1IRI*c1IVI +
        c1IRI)*m475)*c1BBV)*(c1BBV*(c1IVI - c2BBR + 1) - c1IVI*c2BBR)/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)*c1BBV) ** 2) ** 2 + c1BBVe ** 2*(((c1IRI*c1IVI + c1IRI)*c0BBR -
        (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR +
        (c1IRI*c1IVI + c1IRI)*m475)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI -
        (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) - (c0BBV*(c1IRI + 1)*c1IVI*c2BBR +
        (c1IRI + 1)*c1IVI*c2BBR*m475 - ((c1IRI*c1IVI + c1IRI)*c0BBR -
        (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR +
        (c1IRI*c1IVI + c1IRI)*m475)*c1BBV)*(c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR
        + c1IRI)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR
        + c1IRI)*c1BBV) ** 2) ** 2 + ((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI +
        c1IRI)*c1BBV) ** 2*e475 ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI
        - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2)

    Ve = np.sqrt(c1IVIe ** 2*((c1IRI*c2BBR*m475 + (c1IRI*c2BBR - c1IRI)*c0BBV -
        (c0BBR*c1IRI + (c0IRI + m814)*c2BBR + c1IRI*m475)*c1BBV + c0BBR*c1IRI +
        (c0IRI + m814)*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) + ((c1IRI*c1IVI + c1IRI)*c2BBR*m475 +
        (c1IRI*c1IVI + c1IRI)*c0BBR - (c1IRI*c1IVI - (c1IRI*c1IVI + c1IRI)*c2BBR
        + c1IRI)*c0BBV - ((c1IRI*c1IVI + c1IRI)*c0BBR - (c0IVI*(c1IRI + 1) -
        (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR + (c1IRI*c1IVI +
        c1IRI)*m475)*c1BBV - (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI +
        c1IRI*m814 - c0IRI)*c2BBR)*(c1BBV*(c1IRI + c2BBR) - (c1IRI +
        1)*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2) ** 2 + c1IRIe ** 2*(((c1IVI + 1)*c2BBR*m475 +
        ((c1IVI + 1)*c2BBR - c1IVI - 1)*c0BBV - (c0BBR*(c1IVI + 1) - (c0IVI +
        m814)*c2BBR + (c1IVI + 1)*m475)*c1BBV + c0BBR*(c1IVI + 1) - (c0IVI +
        m814)*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) + ((c1IRI*c1IVI + c1IRI)*c2BBR*m475 +
        (c1IRI*c1IVI + c1IRI)*c0BBR - (c1IRI*c1IVI - (c1IRI*c1IVI + c1IRI)*c2BBR
        + c1IRI)*c0BBV - ((c1IRI*c1IVI + c1IRI)*c0BBR - (c0IVI*(c1IRI + 1) -
        (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR + (c1IRI*c1IVI +
        c1IRI)*m475)*c1BBV - (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI +
        c1IRI*m814 - c0IRI)*c2BBR)*(c1BBV*(c1IVI - c2BBR + 1) -
        c1IVI*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2) ** 2 + c1BBVe ** 2*(((c1IRI*c1IVI +
        c1IRI)*c0BBR - (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 -
        c0IRI)*c2BBR + (c1IRI*c1IVI + c1IRI)*m475)/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) - ((c1IRI*c1IVI +
        c1IRI)*c2BBR*m475 + (c1IRI*c1IVI + c1IRI)*c0BBR - (c1IRI*c1IVI -
        (c1IRI*c1IVI + c1IRI)*c2BBR + c1IRI)*c0BBV - ((c1IRI*c1IVI +
        c1IRI)*c0BBR - (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 -
        c0IRI)*c2BBR + (c1IRI*c1IVI + c1IRI)*m475)*c1BBV - (c0IVI*(c1IRI + 1) -
        (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR)*(c1IRI*c1IVI - (c1IRI
        - c1IVI)*c2BBR + c1IRI)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI
        - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2) ** 2 + c2BBRe ** 2*(((c1IRI*c1IVI +
        c1IRI)*c0BBV + (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 -
        c0IRI)*c1BBV - c0IVI*(c1IRI + 1) + (c0IRI + m814)*c1IVI + (c1IRI*c1IVI +
        c1IRI)*m475 - c1IRI*m814 + c0IRI)/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) - ((c1IRI*c1IVI +
        c1IRI)*c2BBR*m475 + (c1IRI*c1IVI + c1IRI)*c0BBR - (c1IRI*c1IVI -
        (c1IRI*c1IVI + c1IRI)*c2BBR + c1IRI)*c0BBV - ((c1IRI*c1IVI +
        c1IRI)*c0BBR - (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 -
        c0IRI)*c2BBR + (c1IRI*c1IVI + c1IRI)*m475)*c1BBV - (c0IVI*(c1IRI + 1) -
        (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR)*(c1BBV*(c1IRI - c1IVI)
        + (c1IRI + 1)*c1IVI)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2) ** 2 + ((c1IRI*c1IVI + c1IRI)*c1BBV -
        c1IRI*c1IVI - c1IRI) ** 2*c0BBRe ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI
        - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (c1IRI*c1IVI - (c1IRI*c1IVI
        + c1IRI)*c2BBR + c1IRI) ** 2*c0BBVe ** 2/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (c1BBV*(c1IVI +
        1)*c2BBR - (c1IVI + 1)*c2BBR) ** 2*c0IRIe ** 2/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (c1BBV*(c1IRI +
        1)*c2BBR - (c1IRI + 1)*c2BBR) ** 2*c0IVIe ** 2/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + ((c1IRI*c1IVI +
        c1IRI)*c1BBV - (c1IRI*c1IVI + c1IRI)*c2BBR) ** 2*e475 ** 2/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2
        + (c1BBV*(c1IRI - c1IVI)*c2BBR - (c1IRI - c1IVI)*c2BBR) ** 2*e814 ** 2/((c1IRI
        + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)*c1BBV) ** 2)

    Re = np.sqrt(((c1IVI + 1)*c2BBR - c1IVI - 1) ** 2*c0IRIe ** 2*c1BBV ** 2/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2
        + ((c1IRI + 1)*c2BBR - c1IRI - 1) ** 2*c0IVIe ** 2*c1BBV ** 2/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2
        + ((c1IRI - c1IVI)*c2BBR - c1IRI + c1IVI) ** 2*c1BBV ** 2*e814 ** 2/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2
        + c1BBVe ** 2*((c0BBR*(c1IRI + 1)*c1IVI + (c1IRI + 1)*c1IVI*m475 +
        c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI - (c0IVI*(c1IRI + 1) - (c0IRI +
        m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR + c1IRI*m814 - c0IRI)/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) -
        ((c1IRI + 1)*c1IVI*c2BBR*m475 + c0BBR*(c1IRI + 1)*c1IVI + ((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI + 1)*c1IVI)*c0BBV - (c0BBR*(c1IRI + 1)*c1IVI +
        (c1IRI + 1)*c1IVI*m475 + c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI -
        (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR +
        c1IRI*m814 - c0IRI)*c1BBV)*(c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)*c1BBV) ** 2) ** 2 + c2BBRe ** 2*((c0BBV*(c1IRI + 1)*c1IVI + (c1IRI +
        1)*c1IVI*m475 + (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 -
        c0IRI)*c1BBV)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) - ((c1IRI + 1)*c1IVI*c2BBR*m475 +
        c0BBR*(c1IRI + 1)*c1IVI + ((c1IRI + 1)*c1IVI*c2BBR - (c1IRI +
        1)*c1IVI)*c0BBV - (c0BBR*(c1IRI + 1)*c1IVI + (c1IRI + 1)*c1IVI*m475 +
        c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI - (c0IVI*(c1IRI + 1) - (c0IRI +
        m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR + c1IRI*m814 -
        c0IRI)*c1BBV)*(c1BBV*(c1IRI - c1IVI) + (c1IRI + 1)*c1IVI)/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)*c1BBV) ** 2) ** 2 + c1IVIe ** 2*(((c1IRI + 1)*c2BBR*m475 + ((c1IRI +
        1)*c2BBR - c1IRI - 1)*c0BBV - (c0BBR*(c1IRI + 1) + (c0IRI + m814)*c2BBR
        + (c1IRI + 1)*m475 - c0IRI - m814)*c1BBV + c0BBR*(c1IRI + 1))/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) +
        ((c1IRI + 1)*c1IVI*c2BBR*m475 + c0BBR*(c1IRI + 1)*c1IVI + ((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI + 1)*c1IVI)*c0BBV - (c0BBR*(c1IRI + 1)*c1IVI +
        (c1IRI + 1)*c1IVI*m475 + c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI -
        (c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR +
        c1IRI*m814 - c0IRI)*c1BBV)*(c1BBV*(c1IRI + c2BBR) - (c1IRI +
        1)*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2) ** 2 + c1IRIe ** 2*((c1IVI*c2BBR*m475 +
        (c1IVI*c2BBR - c1IVI)*c0BBV - (c0BBR*c1IVI - (c0IVI + m814)*c2BBR +
        c1IVI*m475 + c0IVI + m814)*c1BBV + c0BBR*c1IVI)/((c1IRI + 1)*c1IVI*c2BBR
        - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) + ((c1IRI +
        1)*c1IVI*c2BBR*m475 + c0BBR*(c1IRI + 1)*c1IVI + ((c1IRI + 1)*c1IVI*c2BBR
        - (c1IRI + 1)*c1IVI)*c0BBV - (c0BBR*(c1IRI + 1)*c1IVI + (c1IRI +
        1)*c1IVI*m475 + c0IVI*(c1IRI + 1) - (c0IRI + m814)*c1IVI - (c0IVI*(c1IRI
        + 1) - (c0IRI + m814)*c1IVI + c1IRI*m814 - c0IRI)*c2BBR + c1IRI*m814 -
        c0IRI)*c1BBV)*(c1BBV*(c1IVI - c2BBR + 1) - c1IVI*c2BBR)/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)*c1BBV) ** 2) ** 2 + (c1BBV*(c1IRI + 1)*c1IVI - (c1IRI +
        1)*c1IVI) ** 2*c0BBRe ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + ((c1IRI + 1)*c1IVI*c2BBR - (c1IRI +
        1)*c1IVI) ** 2*c0BBVe ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (c1BBV*(c1IRI + 1)*c1IVI - (c1IRI +
        1)*c1IVI*c2BBR) ** 2*e475 ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI -
        (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2)

    Ie = np.sqrt((c1IRI*c2BBR - c1IRI) ** 2*c0IVIe ** 2*c1BBV ** 2/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2
        + c1BBVe ** 2*((c0BBR*c1IRI*c1IVI + c1IRI*c1IVI*m475 + c0IVI*c1IRI -
        (c0IVI*c1IRI - (c0IRI + m814)*c1IVI + c1IRI*m814)*c2BBR +
        c1IRI*m814)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) - (c1IRI*c1IVI*c2BBR*m475 +
        c0BBR*c1IRI*c1IVI + (c0IRI + m814)*c1IVI*c2BBR + (c1IRI*c1IVI*c2BBR -
        c1IRI*c1IVI)*c0BBV - (c0BBR*c1IRI*c1IVI + c1IRI*c1IVI*m475 + c0IVI*c1IRI
        - (c0IVI*c1IRI - (c0IRI + m814)*c1IVI + c1IRI*m814)*c2BBR +
        c1IRI*m814)*c1BBV)*(c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)/((c1IRI
        + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)*c1BBV) ** 2) ** 2 + c2BBRe ** 2*((c0BBV*c1IRI*c1IVI + c1IRI*c1IVI*m475 +
        (c0IVI*c1IRI - (c0IRI + m814)*c1IVI + c1IRI*m814)*c1BBV + (c0IRI +
        m814)*c1IVI)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) - (c1IRI*c1IVI*c2BBR*m475 +
        c0BBR*c1IRI*c1IVI + (c0IRI + m814)*c1IVI*c2BBR + (c1IRI*c1IVI*c2BBR -
        c1IRI*c1IVI)*c0BBV - (c0BBR*c1IRI*c1IVI + c1IRI*c1IVI*m475 + c0IVI*c1IRI
        - (c0IVI*c1IRI - (c0IRI + m814)*c1IVI + c1IRI*m814)*c2BBR +
        c1IRI*m814)*c1BBV)*(c1BBV*(c1IRI - c1IVI) + (c1IRI + 1)*c1IVI)/((c1IRI +
        1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR +
        c1IRI)*c1BBV) ** 2) ** 2 + c1IVIe ** 2*((c1IRI*c2BBR*m475 + (c1IRI*c2BBR -
        c1IRI)*c0BBV - (c0BBR*c1IRI + (c0IRI + m814)*c2BBR + c1IRI*m475)*c1BBV +
        c0BBR*c1IRI + (c0IRI + m814)*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) +
        (c1IRI*c1IVI*c2BBR*m475 + c0BBR*c1IRI*c1IVI + (c0IRI + m814)*c1IVI*c2BBR
        + (c1IRI*c1IVI*c2BBR - c1IRI*c1IVI)*c0BBV - (c0BBR*c1IRI*c1IVI +
        c1IRI*c1IVI*m475 + c0IVI*c1IRI - (c0IVI*c1IRI - (c0IRI + m814)*c1IVI +
        c1IRI*m814)*c2BBR + c1IRI*m814)*c1BBV)*(c1BBV*(c1IRI + c2BBR) - (c1IRI +
        1)*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2) ** 2 + c1IRIe ** 2*((c1IVI*c2BBR*m475 +
        (c1IVI*c2BBR - c1IVI)*c0BBV - (c0BBR*c1IVI - (c0IVI + m814)*c2BBR +
        c1IVI*m475 + c0IVI + m814)*c1BBV + c0BBR*c1IVI)/((c1IRI + 1)*c1IVI*c2BBR
        - (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) +
        (c1IRI*c1IVI*c2BBR*m475 + c0BBR*c1IRI*c1IVI + (c0IRI + m814)*c1IVI*c2BBR
        + (c1IRI*c1IVI*c2BBR - c1IRI*c1IVI)*c0BBV - (c0BBR*c1IRI*c1IVI +
        c1IRI*c1IVI*m475 + c0IVI*c1IRI - (c0IVI*c1IRI - (c0IRI + m814)*c1IVI +
        c1IRI*m814)*c2BBR + c1IRI*m814)*c1BBV)*(c1BBV*(c1IVI - c2BBR + 1) -
        c1IVI*c2BBR)/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI -
        c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2) ** 2 + (c1BBV*c1IRI*c1IVI -
        c1IRI*c1IVI) ** 2*c0BBRe ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI
        - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (c1IRI*c1IVI*c2BBR -
        c1IRI*c1IVI) ** 2*c0BBVe ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI
        - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (c1BBV*c1IVI*c2BBR -
        c1IVI*c2BBR) ** 2*c0IRIe ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI - (c1IRI
        - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (c1BBV*c1IRI*c1IVI -
        c1IRI*c1IVI*c2BBR) ** 2*e475 ** 2/((c1IRI + 1)*c1IVI*c2BBR - (c1IRI*c1IVI -
        (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2 + (((c1IRI - c1IVI)*c2BBR -
        c1IRI)*c1BBV + c1IVI*c2BBR) ** 2*e814 ** 2/((c1IRI + 1)*c1IVI*c2BBR -
        (c1IRI*c1IVI - (c1IRI - c1IVI)*c2BBR + c1IRI)*c1BBV) ** 2)
    return V, Be, V, Ve, R, Re, I, Ie
