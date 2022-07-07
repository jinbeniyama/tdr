#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Cut edge region of 3-d fits. 

1. Cut non-sensitive resions
from 
Toral pixels     2220 x 1360
to 
Sensitive pixels 2160 x 1280
"""
from argparse import ArgumentParser as ap
import astropy.io.fits as fits
import sys
import numpy as np
import os
import datetime 


def main(args):
  """This is the main function called by the `cutedge_tomoe` script.

  Parameters
  ----------
  args : argparse.Namespace
    Arguments passed from the command-line as defined below.
  """
  
  # Create a directory to save output fits
  outdir = "cutedge"
  if os.path.isdir(outdir):
    print(f"Already exists {outdir}")
  else:
    os.makedirs(outdir)
  
  fitsname = os.path.basename(args.fits)
  basename = fitsname.split(".")[0]
  print(basename)

  # Select sensitive pixels
  xmin, xmax = args.xr
  ymin, ymax = args.yr

  # Open a 3-d fits
  # Header keywords are optimized for TriCCS
  hdu = fits.open(args.fits)
  hdr = hdu[0].header
  nz0, ny0, nx0 = hdu[0].data.shape
  print(f"  Original Data Shape (nx,ny,nz)=({nx0},{ny0},{nz0})")

  # Use sensitive pixels
  #temp = data[0:100, (ymin-1):ymax, (xmin-1):xmax]
  hdu[0].data = hdu[0].data[:, ymin:ymax, xmin:xmax]
  nz1, ny1, nx1 = hdu[0].data.shape
  print(f"  Reduced Data Shape (nx,ny,nz)=({nx1},{ny1},{nz1})")
  print(f"    x : {xmin:4d}-{xmax:4d}")
  print(f"    y : {ymin:4d}-{ymax:4d}")

  # Add header keyword
  hdu[0].header.add_history(
    f"[cutedge_cube] original fits: {fitsname}")
  hdu[0].header.add_history(
    f"[cutedge_cube] dim: ({nx0},{ny0},{nz0}) to ({nx1},{ny1},{nz1})")

  # Save the fits
  out = f"{basename}_nx{nx1}ny{ny1}nz{nz1}.fits"
  hdu.writeto(os.path.join(outdir, out), overwrite=True)


if __name__ == "__main__":
  parser = ap(
    description="Cut edge of 3d fits cube")
  parser.add_argument(
    "fits", type=str, 
    help="a reduced 3-d fits")
  parser.add_argument(
    "--xr", nargs=2, type=int, default=[0, 2000], 
    help="x range")
  parser.add_argument(
    "--yr", nargs=2, type=int, default=[0, 1128], 
    help="y range")
  args = parser.parse_args()

  main(args)
