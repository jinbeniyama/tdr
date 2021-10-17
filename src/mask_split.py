#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Split 3-d fits to 2-d fits while cutting and masking edge region 
(i.e. non-sensitive pixels and edge region).

1. Cut non-sensitive resions
from 
Toral pixels     2220 x 1360
to 
Sensitive pixels 2160 x 1280

2. Mask edge region (optional)
Mask not well corrected pixels 
(selected from 2021.Feb and Mar engineering data by J.BENIYAMA)
x 0:10
y 0:10
x 0:150 y 0:150
x 0:150 y 1130:1280
x 2010:2160 y 0:150
x 2010:2160 y 1130:1280
"""
from argparse import ArgumentParser as ap
import astropy.io.fits as fits
import sys
import numpy as np
import os
import datetime 


__naxis3_keywords = (
  'NAXIS3', 'CTYPE3', 'CRPIX3', 'CRVAL3', 'CUNIT3',
  'CD1_3', 'CD2_3', 'CD3_3', 'CD3_2', 'CD3_1',
)


def TriCCSmask():
  """Return boolian mask to be masked. This is easy way to paste correct wcs 
  while suppessing not well corrected pixels. (in J.BENIYAMA's knowledge)
  """
  # Total sensitive pixels
  nx, ny = 2160, 1280
  arr_temp = np.zeros(shape=(ny, nx))
  arr_temp[0:10,:] = 1
  arr_temp[ny-10:ny,:] = 1
  arr_temp[:,0:10] = 1
  arr_temp[:,nx-10:nx] = 1
  arr_temp[0:150, 0:150] = 1
  arr_temp[0:150, nx-150:nx] = 1
  arr_temp[ny-150:ny, 0:150] = 1
  arr_temp[ny-150:ny, nx-150:nx] = 1
  return arr_temp


def main(args):
  """This is the main function called by the `mask_split` script.

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
  basename = fitsname.split("fits")[0]

  # Temporally determine a used filter from file name.
  # Of course, header information is useful as well.
  band = basename.split("TRCS")[1][7]

  # Select sensitive pixels
  if band=="2":
    xmin, xmax = 1, 2160
    ymin, ymax = 1, 1280

  elif (band=="0") or (band=="1"):
    xmin, xmax = 61, 2220
    ymin, ymax = 1, 1280

  # Open a 3-d fits
  # Header keywords are optimized for TriCCS
  hdu = fits.open(args.fits)
  hdr = hdu[0].header
  tframe = hdr["TFRAME"]
  t_exp = hdr["EXPTIME"]
  t_exp = -t_exp

  # Obtain an exposure ending/starting time 
  try:
    ## For old one (ending)
    ### UTC     = '2021-03-08T16:11:32.555790' / exposure ending date and time
    t0 = hdr["UTC"] 
  except:
    ## For new one (starting)
    ### DATE-OBS= '2021-10-15'         / Observation start date (yyyy-mm-dd)
    ### UT-STR  = '14:40:07.40'        / UT at exposure start (HH:MM:SS.SS)
    t0 = f"{hdr['DATE-OBS']}T{hdr['UT-STR']}"

  t0_dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%f")
  t0_dt = t0_dt + datetime.timedelta(seconds=t_exp)
  data_temp = hdu[0].data
  mask = TriCCSmask()

  # Mask edge and split for 3-d fits
  if len(data_temp.shape)==3:
    nz, ny, nx = data_temp.shape
    # Remove useless header keywords
    hdr.set("NAXIS", 2)
    for key in __naxis3_keywords: hdr.remove(key, ignore_missing=True)
    hdu[0].header.add_history(
      f"[mask_split] original fits : {fitsname}")
    print(f"  Split to {nz} fits")
    for i in range(nz):
      t_dt_temp = t0_dt + datetime.timedelta(seconds=tframe*i)
      t_temp = datetime.datetime.strftime(t_dt_temp, "%Y-%m-%dT%H:%M:%S.%f")
      hdr["UTC"] = (t_temp, "exposure starting date and time")
      # Use sensitive pixels
      temp = data_temp[i, (ymin-1):ymax, (xmin-1):xmax]
      if args.mask:
        # Mask not well corrected pixels
        temp = np.where(mask==1, 1.0, temp)
      hdu[0].data = temp
      out = f"{basename}ms{i+1:04d}.fits"
      hdu.writeto(os.path.join(outdir, out), overwrite=True)

  # Mask edge only for 2-d fits
  if len(data_temp.shape)==2:
    # Use sensitive pixels
    temp = data_temp[(ymin-1):ymax, (xmin-1):xmax]
    if args.mask:
      # Mask not well corrected pixels
      temp = np.where(mask==1, 1.0, temp)
    hdu[0].data = temp
    hdu[0].header.add_history(
      f"[mask_split] original fits : {fitsname}")
    out = f"{filename}_c.fits"
    hdu.writeto(os.path.join(out), overwrite=True)


if __name__ == "__main__":
  parser = ap(
    description="Cut and mask not pixels and split to 2-d fits")
  parser.add_argument(
    "fits", type=str, 
    help="a reduced 3-d fits")
  parser.add_argument(
    "--mask", type=bool, default=None, 
    help="mask not well corrected pixes")
  args = parser.parse_args()

  main(args)
