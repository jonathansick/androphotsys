#!/usr/bin/env python
# encoding: utf-8
"""
Model the photometric transformation between F110W-F160W and JKs.
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
    bands = ['2mass_j', '2mass_ks', 'wfc3_ir_f110w', 'wfc3_ir_f160w']
    sp = fsps.StellarPopulation(compute_vega_mags=False, sfh=0)
    data = [sp.get_mags(zmet=i, bands=bands) for i in xrange(1, 23)]
    data = np.vstack(data)

    J = data[:, 0]
    K = data[:, 1]
    F110W = data[:, 2]
    F160W = data[:, 3]

    j_result = scipy.optimize.leastsq(
        cfunc_J,
        [0., 0.1, 0.1], args=(J, F110W, F160W))
    k_result = scipy.optimize.leastsq(
        cfunc_K,
        [0., 0.1, 0.1], args=(K, F110W, F160W))
    print "J coefficients", j_result[0]
    print "Ks coefficients", k_result[0]

    colnames = ['TMAG', 'SMAG', 'COL', 'x0', 'x1', 'x2']
    tmags = ['J', 'Ks']
    smags = ['F160W', 'F160W']
    colors = ['F110W-F160W', 'F110W-F160W']
    x0 = [j_result[0][0], k_result[0][0]]
    x1 = [j_result[0][1], k_result[0][1]]
    x2 = [j_result[0][2], k_result[0][2]]
    tbl = Table([tmags, smags, colors, x0, x1, x2], names=colnames)
    data_dir = os.path.join("androphotsys/data/")
    name = "wfc3_wircam.txt"
    if os.path.exists(data_dir):
        path = os.path.join(data_dir, name)
    else:
        path = name
    tbl.write(path, format='ascii.commented_header')

    fig = Figure(figsize=(4, 7))
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(
        3, 1, left=0.2, right=0.95, bottom=0.15, top=0.95,
        wspace=None, hspace=0.3, width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0], aspect='equal')
    axJ = fig.add_subplot(gs[1])
    axK = fig.add_subplot(gs[2])

    ax.scatter(J - K, F110W - F160W, s=2, edgecolors='None', facecolors='k')
    ax.set_xlabel(r"$J-K_s$")
    ax.set_ylabel("F110W-F160W")

    jc = j_result[0]
    axJ.scatter(
        J - K, cfunc_J(jc, J, F110W, F160W),
        s=2, edgecolors='None', facecolors='k')
    axJ.set_xlabel(r"$J-K_s$")
    axJ.set_ylabel(r"$\Delta_J$")

    kc = k_result[0]
    axK.scatter(
        J - K, cfunc_K(kc, K, F110W, F160W),
        s=2, edgecolors='None', facecolors='k')
    axK.set_xlabel(r"$J - K_s$")
    axK.set_ylabel(r"$\Delta_{K_s}$")

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure("ir_cal.pdf", format="pdf")


def cfunc_J(x, J, F110W, F160W):
    C = F110W - F160W
    return F160W + x[0] + x[1] * C + x[2] * C * C - J


def cfunc_K(x, K, F110W, F160W):
    """0 = F110W + c0 + c1 (J-K) _ c2 (J-K)^2 - J"""
    C = F110W - F160W
    return F160W + x[0] + x[1] * C + x[2] * C * C - K


if __name__ == '__main__':
    main()
