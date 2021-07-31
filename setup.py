#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup
from glob import glob
from os.path import basename
from os.path import splitext



if __name__ == "__main__":

  # Read README.md
  with open("README.md", "r") as f:
    long_description = f.read()

  author = "Jin BENIYAMA"
  email = "beniyama@ioa.s.u-tokyo.ac.jp"

  # from https://pypi.org/classifiers/
  classifiers = [
    "Development Status :: 1 - Planning",
    "Programming Language :: Python :: 3.7",
    "License :: OSI Approved :: MIT License"
    ]
  
  # Upload scripts in 'src' as 'movphot'
  setup(
    name="tdr",
    version="0.0.1",
    author=author,
    author_email=email,
    maintainer=author,
    maintainer_email=email,
    url="https://bitbucket.org/jin_beniyama/triccs-data-reduction/src/master/",
    packages=["tdr"],
    package_dir={"tdr": "src"},
    # Define condole scripts
    entry_points = {
      "console_scripts": ["makedark = tdr.makedark:main",
                          "makeflat = tdr.makeflat:main",
                          "reduce = tdr.reduce:main",
                          "stackfits = tdr.stackfits:main",
                          "commonIDsearch = tdr.commonIDserach:main",
                          ]
    },
    description="TriCCS Data Reduction",
    long_description = long_description,
    long_description_content_type="text/markdown"
    )
