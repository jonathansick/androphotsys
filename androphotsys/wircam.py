#!/usr/bin/env python
# encoding: utf-8
"""
Convert WIRCam photometry.
"""

# Constants added to vega mags to obtain AB mags
# From Blanton et al 2005 AJ 129 2565
WIRCAM_VEGA_TO_AB = {
    "J": 0.91,
    "H": 1.39,
    "Ks": 1.85}
