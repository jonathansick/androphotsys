#!/usr/bin/env python
# encoding: utf-8
"""
Convert WIRCam photometry.
"""

# Constants added to vega mags to obtain AB mags
# From Blanton et al 2005 AJ 129 2565 eq 5.
WIRCAM_VEGA_TO_AB = {
    "J": 0.91,
    "H": 1.39,
    "Ks": 1.85}


def wircam_vega_to_ab(vegamag, band):
    """Returns the AB mag of a WIRCam AB magnitude.

    Uses zeropoint transformations form Blanton et al 2005 AJ 129 2565 eq 5.
    """
    return vegamag + WIRCAM_VEGA_TO_AB[band]
