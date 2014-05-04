#!/usr/bin/env python
# encoding: utf-8
"""
Photometric conversions from HST bandpasses to the Cousins system.

See: Sirianni 2005.
"""

import numpy as np


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

    assert zp in ('STMAG', 'VEGAMAG', 'OBMAG')
    if zp == 'STMAG':
        # Subtract the STMAG zeropoints to get OBMAG (Sirianni 2005 Table 11)
        m606 = m606_input - 26.655
        m814 = m814_input - 26.776
    elif zp == 'VEGAMAG':
        # Subtract the VEGAMAG zeropoints to get OBMAG (Sirianni 2005 Table 11)
        m606 = m606_input - 26.398
        m814 = m814_input - 25.501
    elif zp == 'ABMAG':
        # Subtract the VEGAMAG zeropoints to get OBMAG (Sirianni 2005 Table 11)
        m606 = m606_input - 26.486
        m814 = m814_input - 25.937
    elif zp == 'OBMAG':
        m606 = m606_input
        m814 = m814_input

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
    assert zp in ('STMAG', 'VEGAMAG', 'OBMAG', 'ABMAG')
    if zp == 'STMAG':
        # Subtract the STMAG zeropoints to get OBMAG (Sirianni 2005 Table 11)
        m475 = m475_input - 25.757
        m814 = m814_input - 26.776
    elif zp == 'VEGAMAG':
        # Subtract the VEGAMAG zeropoints to get OBMAG (Sirianni 2005 Table 11)
        m475 = m475_input - 26.168
        m814 = m814_input - 25.501
    elif zp == 'ABMAG':
        # Subtract the VEGAMAG zeropoints to get OBMAG (Sirianni 2005 Table 11)
        m475 = m475_input - 26.068
        m814 = m814_input - 25.937
    elif zp == 'OBMAG':
        m475 = m475_input
        m814 = m814_input

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
