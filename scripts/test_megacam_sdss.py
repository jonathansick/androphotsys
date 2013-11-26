#!/usr/bin/env python
# encoding: utf-8
"""
Test colour transformation between megacam and SDSS bands.

2013-11-25 - Created by Jonathan Sick
"""

import numpy as np

from androphotsys.elixir import megacam_to_sdss


def main():
    # These are SDSS magnitudes
    # input_sdss_u = 4 * np.random.randn(4)
    # input_sdss_g = 4 * np.random.randn(4)
    # input_sdss_r = 4 * np.random.randn(4)
    # input_sdss_i = 4 * np.random.randn(4)
    input_sdss_u = 18.
    input_sdss_g = 19.
    input_sdss_r = 20.
    input_sdss_i = 21.

    # These are MegaCam magnitudes computed by the forward transformation
    XU = 0.241
    XG = 0.153
    XR = 0.024
    XI = 0.003
    input_mega_u = input_sdss_u - XU * (input_sdss_u - input_sdss_g)
    input_mega_g = input_sdss_g - XG * (input_sdss_g - input_sdss_r)
    input_mega_r = input_sdss_r - XR * (input_sdss_g - input_sdss_r)
    input_mega_i = input_sdss_i - XI * (input_sdss_r - input_sdss_i)

    direct_g_sdss = (input_mega_g - XG * input_sdss_r) / (1. - XG)
    direct_r_sdss = (input_mega_r + XR * input_sdss_g) / (1. + XR)

    trans = megacam_to_sdss(u=input_mega_u, g=input_mega_g,
            r=input_mega_r, i=input_mega_i,
            XU=XU, XG=XG, XR=XR, XI=XI)
    print "Residuals"
    print "u"
    print trans['u'] - input_sdss_u
    print "g"
    print trans['g'] - input_sdss_g
    print input_mega_g - direct_g_sdss
    print input_mega_g - input_sdss_g
    print "r"
    print trans['r'] - input_sdss_r
    print input_mega_r - direct_r_sdss
    print input_mega_r - input_sdss_r
    print "i"
    print trans['i'] - input_sdss_i


if __name__ == '__main__':
    main()
