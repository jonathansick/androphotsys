#!/usr/bin/env python
# encoding: utf-8
"""
Transform the Brown/ACS f606w-f814w photometry to SDSS gri-bands.
"""

import os
import numpy as np
import scipy.optimize
from astropy.table import Table

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec

import fsps


def main():
    bands = ['sdss_g', 'sdss_r', 'sdss_i', 'wfc_acs_f606w', 'wfc_acs_f814w']
    sp = fsps.StellarPopulation(compute_vega_mags=False, sfh=0)
    data = [sp.get_mags(zmet=i, bands=bands) for i in xrange(1, 23)]
    data = np.vstack(data)

    g = data[:, 0]
    r = data[:, 1]
    i = data[:, 2]
    F606W = data[:, 3]
    F814W = data[:, 4]

    g_result = scipy.optimize.leastsq(
        cfunc_g,
        [0., 0.1, 0.1], args=(g, F606W, F814W))
    r_result = scipy.optimize.leastsq(
        cfunc_r,
        [0., 0.1, 0.1], args=(r, F606W, F814W))
    i_result = scipy.optimize.leastsq(
        cfunc_i,
        [0., 0.1, 0.1], args=(i, F606W, F814W))
    print "g coefficients", g_result[0]
    print "r coefficients", r_result[0]
    print "i coefficients", i_result[0]

    colnames = ['TMAG', 'SMAG', 'COL', 'x0', 'x1', 'x2']
    tmags = ['g', 'r', 'i']
    smags = ['F606W', 'F814W', 'F814W']
    colors = ['F606W-F814W', 'F606W-F814W', 'F606W-F814W']
    x0 = [g_result[0][0], r_result[0][0], i_result[0][0]]
    x1 = [g_result[0][1], r_result[0][1], i_result[0][1]]
    x2 = [g_result[0][2], r_result[0][2], i_result[0][2]]
    tbl = Table([tmags, smags, colors, x0, x1, x2], names=colnames)
    data_dir = os.path.join("androphotsys/data/")
    name = "acs_brown.txt"
    if os.path.exists(data_dir):
        path = os.path.join(data_dir, name)
    else:
        path = name
    tbl.write(path, format='ascii.commented_header')

    fig = Figure(figsize=(4, 5))
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(
        3, 1, left=0.2, right=0.95, bottom=0.15, top=0.95,
        wspace=None, hspace=0.3, width_ratios=None, height_ratios=None)
    axg = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1])
    axi = fig.add_subplot(gs[2])

    gc = g_result[0]
    axg.scatter(
        g, cfunc_g(gc, g, F606W, F814W),
        s=2, edgecolors='None', facecolors='k')
    axg.set_xlabel(r"$g$")
    axg.set_ylabel(r"$\Delta_g$")

    rc = r_result[0]
    axr.scatter(
        r, cfunc_r(rc, r, F606W, F814W),
        s=2, edgecolors='None', facecolors='k')
    axr.set_xlabel(r"$r$")
    axr.set_ylabel(r"$\Delta_r$")

    ic = i_result[0]
    axi.scatter(
        i, cfunc_i(ic, i, F606W, F814W),
        s=2, edgecolors='None', facecolors='k')
    axi.set_xlabel(r"$i$")
    axi.set_ylabel(r"$\Delta_i$")

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure("acs_brown_cal.pdf", format="pdf")


def cfunc_g(x, g, F606W, F814W):
    C = F606W - F814W
    return F606W + x[0] + x[1] * C + x[2] * C * C - g


def cfunc_r(x, r, F606W, F814W):
    C = F606W - F814W
    return F814W + x[0] + x[1] * C + x[2] * C * C - r


def cfunc_i(x, i, F606W, F814W):
    C = F606W - F814W
    return F814W + x[0] + x[1] * C + x[2] * C * C - i


if __name__ == '__main__':
    main()
