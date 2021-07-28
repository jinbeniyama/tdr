#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Cut edge and Mask edge of the fits taken with the TriCCS.

mask region
for g,r band
x 0:10
y 0:10
x 0:150 y 0:150
x 0:150 y 1130:1280
x 2010:2160 y 0:150
x 2010:2160 y 1130:1280

for i/z band
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

from myio import get_filename


__naxis3_keywords = (
  'NAXIS3', 'CTYPE3', 'CRPIX3', 'CRVAL3', 'CUNIT3',
  'CD1_3', 'CD2_3', 'CD3_3', 'CD3_2', 'CD3_1',
)


def TriCCSmask():
  """Return boolian mask to be masked.
  This is necessary to paste correct wcs.
  """
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


if __name__ == '__main__':
  parser = ap(description='Cut edge of an image.')
  parser.add_argument(
    'fits', type=str, help='a raw fits image')
  args = parser.parse_args()
 
  filename = get_filename(args.fits)
  band = filename.split("TRCS")[1][7]
  if band=="2":
    xmin, xmax = 1, 2160
    ymin, ymax = 1, 1280

  elif (band=="0") or (band=="1"):
    xmin, xmax = 61, 2220
    ymin, ymax = 1, 1280
  outdir = "cutedge"
  if os.path.isdir(outdir):
    print(f"Already exists {outdir}")
  else:
    os.makedirs(outdir)

  filename = get_filename(args.fits)
  hdu = fits.open(args.fits)
  hdr = hdu[0].header
  tframe = hdr["TFRAME"]
  t_exp = hdr["EXPTIME"]
  t_exp = -t_exp
  t0 = hdr["UTC"]
  t0_dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%f")
  t0_dt = t0_dt + datetime.timedelta(seconds=t_exp)
  data_temp = hdu[0].data

  # cut and separate for 3-d fits
  if len(hdu[0].data.shape)==3:
    mask = TriCCSmask()
    nz, ny, nx = hdu[0].data.shape
    hdr.set('NAXIS', 2)
    for key in __naxis3_keywords: hdr.remove(key, ignore_missing=True)
    hdu[0].header.add_history(
      f'[cutandseparate] created from {filename}.')
    for i in range(nz):
      print(i)
      t_dt_temp = t0_dt + datetime.timedelta(seconds=tframe*i)
      t_temp = datetime.datetime.strftime(t_dt_temp, "%Y-%m-%dT%H:%M:%S.%f")
      hdr['UTC'] = (t_temp, 'exposure starting date and time')
      temp = data_temp[i, (ymin-1):ymax, (xmin-1):xmax]
      temp = np.where(mask==1, 1.0, temp)
      hdu[0].data = temp
      out = f"{filename}ms{i+1:04d}.fits"
      hdu.writeto(os.path.join(outdir, out), overwrite=True)

  # cut for 2-d fits
  if len(hdu[0].data.shape)==2:
    sys.exit()
    hdu[0].data = data_temp[(ymin-1):ymax, (xmin-1):xmax]
    hdu[0].header.add_history(
      f'[cut] created from {filename}.')
    out = f"{filename}_c.fits"
    hdu.writeto(os.path.join(out), overwrite=True)
