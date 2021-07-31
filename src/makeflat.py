#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Make a master flat frame from a 3-d fits
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


def main(args=None):
  """This is the main function called by the `makeflat` script.

  Parameters
  ----------
  args : argparse.Namespace
    Arguments passed from the command-line as defined below.
  """
  parser = ap(description="Make a mater flat frame from a 3-d fits")
  parser.add_argument(
    "--flat", type=str, required=True,
    help="a raw flat frame")
  parser.add_argument(
    "--dark", type=str, required=True,
    help="a raw dark frame")
  parser.add_argument(
    "--out", type=str, default=None,
    help="output fits file")
  parser.add_argument(
    "-f", dest="overwrite", action="store_true",
    help="overwrite a fits image")
  args = parser.parse_args()

  hdu_flat = fits.open(args.flat)
  src_flat = hdu_flat[0]
  flat = src_flat.data
  hdr_flat = src_flat.header
  assert len(flat.shape)==3, "Flat should be 3-d fits."

  hdu_dark = fits.open(args.dark)
  src_dark = hdu_dark[0]
  dark = src_dark.data
  assert len(dark.shape)==2, "Dark should be 2-d fits."
  
  img = chopped_mean(flat) - dark
  # Normalize by mean
  img = img/np.mean(img)

  # Create 2-d master flat
  flat = fits.PrimaryHDU(data=img, header=hdr_flat)
  # Add history
  hdr = flat.header
  hdr.add_history(
    f"[makeflat] created from {os.path.basename(args.flat)}")
  hdr.add_history(
    f"[makeflat] created from {os.path.basename(args.dark)}")


  if args.out is None:
    out = f"f{os.path.join.basename(args.flat)}"
  else:
    out = args.out

  # Write new fits 
  flat.writeto(out, overwrite=args.overwrite)


if __name__ == "__main__":
  parser = ap(description="Make a mater flat frame from a 3-d fits")
  parser.add_argument(
    "--flat", type=str, required=True,
    help="a raw flat frame")
  parser.add_argument(
    "--dark", type=str, required=True,
    help="a raw dark frame")
  parser.add_argument(
    "--out", type=str, default=None,
    help="output fits file")
  parser.add_argument(
    "-f", dest="overwrite", action="store_true",
    help="overwrite a fits image")
  args = parser.parse_args()
  
  main(args)
