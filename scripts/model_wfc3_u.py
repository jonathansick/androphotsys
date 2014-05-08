#!/usr/bin/env python
# encoding: utf-8
"""
Transform the PHAT f275w-f336w photometry to SDSS u-band.
"""

import numpy as np
import scipy.optimize

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec

import fsps


def main():
    bands = ['sdss_u', 'wfc3_uvis_f275w', 'wfc3_uvis_f336w']
    sp = fsps.StellarPopulation(compute_vega_mags=False, sfh=0)
    data = [sp.get_mags(zmet=i, bands=bands) for i in xrange(1, 23)]
    data = np.vstack(data)

    u = data[:, 0]
    F275W = data[:, 1]
    F336W = data[:, 2]

    u_result = scipy.optimize.leastsq(cfunc_u,
            [0., 0.1, 0.1], args=(u, F275W, F336W))
    print "u coefficients", u_result[0]

    fig = Figure(figsize=(4, 5))
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, 1, left=0.2, right=0.95, bottom=0.15, top=0.95,
        wspace=None, hspace=0.3, width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    axJ = fig.add_subplot(gs[1])

    ax.scatter(u, F275W - F336W, s=2, edgecolors='None', facecolors='k')
    ax.set_xlabel(r"$u$")
    ax.set_ylabel("F275W-F336W")

    uc = u_result[0]
    axJ.scatter(u, cfunc_u(uc, u, F275W, F336W),
            s=2, edgecolors='None', facecolors='k')
    axJ.set_xlabel(r"$u$")
    axJ.set_ylabel(r"$\Delta_J$")

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure("uv_cal.pdf", format="pdf")


def cfunc_u(x, u, F275W, F336W):
    C = F275W - F336W
    return F336W + x[0] + x[1] * C + x[2] * C * C - u


if __name__ == '__main__':
    main()
