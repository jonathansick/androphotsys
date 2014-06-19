#!/usr/bin/env python
# encoding: utf-8
"""
Helpers for FSPS filters.
"""


ANDROIDS_TO_FSPS = {'u': 'sdss_u',
                    'g': 'sdss_g',
                    'r': 'sdss_r',
                    'i': 'sdss_i',
                    'J': '2mass_J',
                    'Ks': '2mass_Ks',
                    '36': 'irac_1',
                    '45': 'irac_2',
                    '58': 'irac_3',
                    '80': 'irac_4'}
FSPS_TO_ANDROIDS = {v: k for k, v in ANDROIDS_TO_FSPS.iteritems()}


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
