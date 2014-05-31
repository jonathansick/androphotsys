#!/usr/bin/env python
# encoding: utf-8
"""
Photometric transformations between Johnson Cousins and SDSS magnitudes
"""

import numpy as np


def FIG96_UBVRc_to_ugr(U, B, V,
                       Ue=None, Be=None, Ve=None):
    """
    Conversions to Johnson-Cousins VEGAMAG magnitudes to SDSS ABMAG.

    Fukugita, M., Ichikawa, T., Gunn, J. E. et al
    1996AJ....111.1748F
    """
    g = V + 0.56 * (B - V) - 0.12
    r = V - 0.49 * (B - V) + 0.11
    # Alt: r = V - 0.84 * (V - R) + 0.13
    ug = 1.38 * (U - B) + 1.14
    u = ug + g
    return u, g, r


def KBT05_UBVRc_to_ugr(U, B, V,
                       Ue=None, Be=None, Ve=None):
    """Karaali, Bilir & Tunçel (2005)

    **Not implemented**

    These transformations appeared in Karaali, Bilir & Tunçel (2005).
    They are based on Landolt (1992) UBV data for 224 stars in the color
    range 0.3 < B-V < 1.1 with SDSS ugr photometry from the CASU INT Wide
    Field Survey. An improvement over previous SDSS - UBVRcIc transformations
    is the use of two colors in each equation, which is particularly helpful
    for the u-g transformation.

    UBVRcIc -> ugriz

    Stars with  0.3 < B-V < 1.1
        u-g    =    0.779*(U-B) + 0.755*(B-V)  + 0.801
        g-r    =    1.023*(B-V) + 0.016*(U-B)  - 0.187

    I *think* that that input magnitudes are meant to be in Vega mag, while
    output magnitudes are meant to be AB.
    
    See:
    http://www.sdss3.org/dr10/algorithms/sdssUBVRITransform.php#Bilir2005
    """
    ug = 0.779 * (U - B) + 0.755 * (B - V) + 0.801
    gr = 1.023 * (B - V) + 0.016 * (U - B) - 0.187
    return None


def VRI_to_gri(V, Ve, R, Re, I, Ie):
    """Transformation of Johnson Cousins VRI magnitudes to SDSS AB mags.
    
    Based on Lupton 2005
    http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php#Lupton2005
    
    Solution of::
        V = g - 0.5784*(g - r) - 0.0038
        R = r - 0.2936*(r - i) - 0.1439
        I = r - 1.2444*(r - i) - 0.3820
    """
    g = 530682. / 1252679 * I - 132309. / 73687. * R + 1250. / 527. * V \
        - 1096562327. / 12526790000.
    r = -734. / 2377. * I + 3111. / 2377. * R + 1672849. / 23770000.
    i = 1766. / 2377. * I + 611. / 2377. * R + 7625349. / 23770000.
    g_e = np.sqrt(281623385124. / 1569204677041 * Ie ** 2
            + 17505671481. / 5429773969. * Re ** 2
            + 1562500. / 277729. * Ve ** 2)
    r_e = np.sqrt(538756. / 5650129. * Ie ** 2 + 9678321. / 5650129. * Re ** 2)
    i_e = np.sqrt(3118756. / 5650129. * Ie ** 2 + 373321. / 5650129. * Re ** 2)
    return g, g_e, r, r_e, i, i_e


def BVRI_to_ugri(B, Be, V, Ve, R, Re, I, Ie):
    """Transformation of Johnson Cousins BVRI magnitudes to SDSS AB mags.
    
    Based on Lupton 2005
    http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php#Lupton2005
    
    Solution of::

        B = u - 0.8116*(u - g) + 0.1313
        V = g - 0.5784*(g - r) - 0.0038
        R = r - 0.2936*(r - i) - 0.1439
        I = r - 1.2444*(r - i) - 0.3820
    """
    u = 2500. / 471. * B - 358917926. / 196670603. * I \
        + 89484987. / 11568859. * R \
        - 2536250. / 248217. * V - 1886993856017. / 5900118090000.
    g = 530682. / 1252679. * I - 132309. / 73687. * R + 1250. / 527 * V \
        - 1096562327. / 12526790000.
    r = -734. / 2377. * I + 3111. / 2377. * R + 1672849. / 23770000.
    i = 1766. / 2377. * I + 611. / 2377. * R + 7625349. / 23770000.

    u_e = np.sqrt(6432564062500. / 61611679089. * Be ** 2 +
        128822077604141476. / 38679326084383609. * Ie ** 2 +
        8007562898390169. / 133838498561881. * Re ** 2 +
        6432564062500. / 61611679089. * Ve ** 2)
    g_e = np.sqrt(1562500. / 277729. * Be ** 2
        + 281623385124. / 1569204677041. * Ie ** 2
        + 17505671481. / 5429773969. * Re ** 2 + 1562500. / 277729. * Ve ** 2)
    r_e = np.sqrt(538756. / 5650129. * Ie ** 2 + 9678321. / 5650129. * Re ** 2)
    i_e = np.sqrt(3118756. / 5650129. * Ie ** 2 + 373321. / 5650129. * Re ** 2)
    return u, u_e, g, g_e, r, r_e, i, i_e
