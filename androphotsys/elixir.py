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
    if u is not None:
        # depends on u, g, r
        u_sdss = u / (1. - XU) \
            - XU / (1. - XU) * ((g * (1. + XR) - XG * r) \
            / ((1. - XG) * (1. - XR) + XG * XR))
    else:
        u_sdss = None

    if g is not None:
        # depends on g, r
        g_sdss = (g * (1. + XR) - XG * r) / ((1. - XG) * (1. - XR) + XG * XR)
    else:
        g_sdss = None

    if r is not None:
        # depends on g, r
        r_sdss = r / (1. + XR) + XG * g / ((1. - XG) * (1. - XR) + XG * XR) \
                - XG ** 2. * r / (1. + XR) / ((1. - XG) * (1. - XR) + XG * XR)
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
