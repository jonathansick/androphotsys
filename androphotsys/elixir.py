#!/usr/bin/env python
# encoding: utf-8
"""
Photometric system conversions for MegaCam/Elixir-LSB data.

2013-11-21 - Created by Jonathan Sick
"""


def megacam_to_sdss(u=None, g=None, r=None, i=None,
        u_err=None, g_err=None, r_err=None, i_err=None,
        XU=0.241, XG=0.153, XR=0.024, XI=0.003):
    """Convert MegaCam magnitudes to SDSS magnitudes (both in AB system)."""

    if g is not None:
        # depends on g, r
        g1 = g * (1. - XR) - r * XG
        g2 = (1. - XG) * (1. - XR) - XR * XG
        g_sdss = g1 / g2
    else:
        g_sdss = None

    if u is not None:
        # depends on u, g, r
        u_sdss = (u - XU * g_sdss) / (1. - XU)
    else:
        u_sdss = None

    if r is not None:
        # depends on g, r
        r_sdss = (r + XR * g_sdss) / (1. + XR)
    else:
        r_sdss = None

    if i is not None:
        # depends on g, r, i
        i_sdss = i / (1. + XI) + XI * r / (1. + XR) \
            + XI * XR / (1. + XR) \
            * (g * (1. + XR) - XG * r) / ((1. - XG) * (1. - XR) - XG * XR)
    else:
        i_sdss = None

    return {"u": u_sdss, "g": g_sdss, "r": r_sdss, "i": i_sdss}
