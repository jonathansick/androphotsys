#!/usr/bin/env python
# encoding: utf-8
"""
Write a latex table of HST zeropoints.
"""

from androphotsys.hst import load_colour_transformation_table


def main():
    table = load_colour_transformation_table()
    names = ['TMAG', 'SMAG', r"CMAG", r"$x_0$", r"$x_1$", r"$x_2$"]
    latexfmt = {'tablealign': 't',
                'col_align': 'lll|rrr',
                'caption': 'HST color transformations modelled using FSPS '
                           '(see text) to yield ANDROIDS bandpasses. Both '
                           'HST and SDSS magnitudes use the AB zeropoint.',
                'preamble': r'\begin{center}',
                'tablefoot': r'\end{center}\label{tab:hst_ct}',
                'headerstart': r'\hline',
                'headerend': r'\hline',
                'dataend': r'\hline'}
    fmt = {'$x_0$': '{:.3f}', '$x_1$': '{:.3f}', '$x_2$': '{:.3f}'}
    table.write("hst_ct.tex",
                format="ascii.latex",
                latexdict=latexfmt,
                names=names,
                formats=fmt)


if __name__ == '__main__':
    main()
