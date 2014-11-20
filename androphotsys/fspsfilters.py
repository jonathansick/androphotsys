#!/usr/bin/env python
# encoding: utf-8
"""
Helpers for FSPS filters.
"""


ANDROIDS_TO_FSPS = {'u': 'megacam_u',
                    'g': 'megacam_g',
                    'r': 'megacam_r',
                    'i': 'megacam_i',
                    'J': '2mass_J',
                    'Ks': '2mass_Ks',
                    '36': 'irac_1',
                    '45': 'irac_2',
                    '58': 'irac_3',
                    '80': 'irac_4',
                    'f475w': 'wfc_acs_f475w',
                    'f814w': 'wfc_acs_f814w',
                    'f110w': 'wfc3_ir_f110w',
                    'f160w': 'wfc3_ir_f160w'}
FSPS_TO_ANDROIDS = {v: k for k, v in ANDROIDS_TO_FSPS.iteritems()}


FILTER_LATEX = {"megacam_u": r'$u^*$',
                "megacam_g": r'$g^\prime$',
                "megacam_r": r'$r^\prime$',
                "megacam_i": r'$i^\prime$',
                "sdss_u": r'$u$',
                "sdss_g": r'$g$',
                "sdss_r": r'$r$',
                "sdss_i": r'$i$',
                "sdss_z": r'$z$',
                "2mass_J": r'$J$',
                "2mass_Ks": r'$K_s$',
                "irac_1": r'$[3.6]$',
                "irac_2": r'$[4.5]$',
                "irac_3": r'$[5.8]$',
                "irac_4": r'$[8.0]$',
                "f275w": r'F275W',
                "f336w": r'F336W',
                "f475w": r'F475W',
                "f814w": r'F814W',
                "f110w": r'F110W',
                "f160w": r'F160W'}


def fsps_name_to_androids(fsps_filter_name):
    """Map a filter key name for FSPS to conventional ANDROIDS naming."""
    if fsps_filter_name in FSPS_TO_ANDROIDS:
        return FSPS_TO_ANDROIDS[fsps_filter_name]
    else:
        return fsps_filter_name


def androids_name_to_fsps(androids_filter_name):
    """Map a filter key name for ANDROIDS to conventional FSPS naming."""
    if androids_filter_name in ANDROIDS_TO_FSPS:
        return ANDROIDS_TO_FSPS[androids_filter_name]
    else:
        return androids_filter_name


def latex_name(filter_name, mathmode=True):
    """LateX string for the named filter.

    Parameters
    ----------
    filter_name : str
        An ANDROIDS or FSPS filter name.
    mathmode : bool
        If `True` (default) then math delimeters (``$``) are used as
        appropriate. If `False` then those demimeters are removed.
    """
    if filter_name in FILTER_LATEX:
        lstr = FILTER_LATEX[filter_name]
    elif androids_name_to_fsps(filter_name) in FILTER_LATEX:
        lstr = FILTER_LATEX[androids_name_to_fsps(filter_name)]
    else:
        assert False
    if not mathmode:
        return lstr.replace("$", "")
    else:
        return lstr
