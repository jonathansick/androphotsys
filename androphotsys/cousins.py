#!/usr/bin/env python
# encoding: utf-8
"""
Photometric transformations between johnson cousins and other magnitude
systems.
"""

import numpy as np


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
