#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Make a master dark frame from a 3-d fits
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits


def chopped_mean(cube):
  """ Calculate chopped mean along axis-zero

  Parameter
  ---------
    cube : numpy.ndarray
      three-dimensional ndarray

  Return
  ------
    img : numpy.ndarray
      stacked two-dimensional image
  """
  nz, ny, nx = cube.shape
  img = np.sum(cube, axis=0)-np.max(cube, axis=0)
  img = img/(nz-1.0)
  return img


def main(args):
  """This is the main function called by the `makedark` script.

  Parameters
  ----------
  args : argparse.Namespace
    Arguments passed from the command-line as defined below.
  """
  hdu = fits.open(args.fits)
  src = hdu[0]
  cube = src.data
  hdr = src.header
  assert len(cube.shape)==3, "Input should be 3-d fits."

  img = chopped_mean(cube)

  # Create 2-d master dark
  dark = fits.PrimaryHDU(data=img, header=hdr)
  # Add history
  hdr = dark.header
  hdr.add_history(
    f"[makedark] original dark frames : {os.path.basename(args.fits)}")

  if args.out is None:
    out = f"d{os.path.basename(args.fits)}"
  else:
    out = args.out

  # Write new fits 
  dark.writeto(out, overwrite=args.overwrite)


if __name__ == "__main__":
  parser = ap(description="Make a mater dark frame from a 3-d fits")
  parser.add_argument(
    "fits", type=str,
    help="fits file to compile a dark-frame")
  parser.add_argument(
    "--out", type=str, default=None,
    help="output fits file")
  parser.add_argument(
    "-f", dest="overwrite", action="store_true",
    help="overwrite a fits image")
  args = parser.parse_args()
  
  main(args)
