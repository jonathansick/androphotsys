#!/usr/bin/env python
# encoding: utf-8
"""
Fit the ZP offset between WIRCam/2MASS Vega (native) and AB zeropoints.
"""
import fsps


def main():
    bands = ['2mass_j', '2mass_ks']
    vega_sp = fsps.StellarPopulation(compute_vega_mags=True, sfh=0, tage=10.)
    vega_19 = vega_sp.get_mags(zmet=19, bands=bands)

    fsps.fsps.driver.is_setup = False
    ab_sp = fsps.StellarPopulation(compute_vega_mags=False, sfh=0, tage=10.)
    ab_19 = ab_sp.get_mags(zmet=19, bands=bands)

    diff = vega_19 - ab_19
    print "VEGAMAG - ABMAG"
    print bands[0], diff[0][0]
    print bands[1], diff[0][1]


if __name__ == '__main__':
    main()
