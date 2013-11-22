#!/usr/bin/env python

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages


setup(
    name = "androphotsys",
    version = "0.1",
    packages = find_packages(),

    # metadata for upload to PyPI
    author = "Jonathan Sick",
    author_email = "jonathansick@mac.com",
    description = "Photometry system for ANDROPHOT/CFHT data.",
    license = "BSD",
    keywords = "astronomy",
    url = "http://www.jonathansick.ca",   # project home page, if any
)
