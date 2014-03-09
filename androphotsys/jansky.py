#!/usr/bin/env python
# encoding: utf-8
"""
Convert photometry to Janskies.
"""


def ab_to_microjy(ab_mag, err=None):
    """Convert AB magnitudes to µJy.

    See http://en.wikipedia.org/wiki/Jansky
    
    Parameters
    ----------
    ab_mag : array
        Magnitudes in the AB system.
    err : array
        Uncertainties in AB magnitudes (optional), interpreted as 1-sigma
        Gaussian errors.

    Returns
    -------
    mjy : array
        Fluxes, as micro-janskies
    mjy_err : array
        (If errors are provided) the 1-sigma uncertainties in flux.
    """
    mjy = 10. ** (-0.4 * (ab_mag + 48.6) + 6. + 23.)
    if err is not None:
        mjy_err = mjy * 10. ** (-0.4 * 48.6) * 0.4 * err  # TODO verify
        return mjy, mjy_err
    else:
        return mjy


def dn_to_microjy(dn, zp, err=None):
    """Convert pixel values (DN) to µJy.

    Parameters
    ----------
    dn : array
        Pixel values.
    zp : array or scalar
        Zeropoint (or zeropoints) corresponding to the DN measurements
        to bring them into the AB system, where

        .. math::
           m_\mathrm{AB} = -2.5 \log_{10} (\mathrm{DN}) + \mathrm{zp}
    err : array
        Uncertainties in DN (optional), interpreted as 1-sigma
        Gaussian errors.

    Returns
    -------
    mjy : array
        Fluxes, as micro-janskies
    mjy_err : array
        (If errors are provided) the 1-sigma uncertainties in flux.
    """
    f = 10. ** (-0.4 * (zp + 48.6) + 6. + 23.)
    mjy = dn * f
    if err is not None:
        mjy_err = err * f
        return mjy, mjy_err
    else:
        return mjy
