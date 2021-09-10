#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Extract 16 (4 times 4) patterns of GPS time from fits taken with TriCCS.
This script is optimized for data after ~2020-02-22T16:00(?).
Check and flip arrays for data before ~2020-02-22T16:00(?).
"""
import os
import numpy as np
import pandas as pd
from argparse import ArgumentParser as ap
import time
import astropy.io.fits as fits


def get_filename(file_path):
  """Get 'hoge' from a path 'huga/hoge.fits'.

  Parameter
  ---------
  file_path : str
    path of the file

  Return
  ------
  filename : str
    basename before extension of the path
  """
  filename = file_path.split("/")[-1].split(".")[0]
  return filename


def extract_GPStime(data):
  """Extract GPS time from 4 times 4 arrays below.
  All data in array are expressed by 16-bit integers.
  x-->
  |  1 |  2 |  3 |  4 | ^
  |  5 |  6 |  7 |  8 | |
  |  9 | 10 | 11 | 12 | |
  | 13 | 14 | 15 | 16 | y
    When X is an index of array(pixel), 
    Y and Z are indexes of a digit in X[Y:Z],
    'second' is expressed as
        2[19:16] + 4[15:0],
    'micro second' is expressed as 
        10[19:16] + 12[15:0].
  GPS time could be calculated by 'second' + 'micro second'.

  Parameter
  ---------
  data : arary-like
    4 times 4 arrays

  Returns
  -------
  sec : int
    second part of GPS time
  musec : int
    micro second part of GPS time
  """

  # Reproduce original value
  data = data*4.
  # Convert to 1-d arrays
  data_1d = data.flatten()
  
  # Second part ===============================================================
  sec1 = data_1d[11]
  sec2 = "{:016b}".format(int(data_1d[9]))
  # Extract last 4-bit values ([12:16]) related to GPS time.
  sec2 = sec2[12:16]
  sec2 = int(sec2, 2)
  # Multiple 2e16 for 16-bit shift
  sec = sec1 + sec2*2**16
  assert sec.is_integer(), f"Second is not an integer. {sec}"
  # Second part finish ========================================================


  # Micro second part =========================================================
  musec1 = data_1d[3]
  musec2 = "{:016b}".format(int(data_1d[1]))
  # Extract last 4-bit values ([12:16]) related to GPS time.
  musec2 = musec2[12:16]
  musec2 = int(musec2, 2)
  # Multiple 2e16 for 16-bit shift
  musec = musec1 + musec2*2**16
  assert musec.is_integer(), f"Micro second is not an integer. {musec}"
  # Micro second part finish ==================================================

  return int(sec), int(musec)
 

def main(args):
  """A main function.
  """
  assert args.fits[0:4]=="TRCS", "Input fits should be obtained with the TriCCS"

  filename = get_filename(args.fits)
  hdu = fits.open(args.fits)
  cube = hdu[0].data
  hdr = hdu[0].header
  nz, ny, nx = cube.shape
  print(f"  Cube shape : (nx,ny,nz)=({nx},{ny},{nz})")

  n_row = ny/4.0
  assert n_row.is_integer(), f"N_row is not multiple of 4."
  n_row = int(n_row)
  
  band = args.fits[11:12]
  print(f"  fits: {args.fits}")
  print(f"  band: {band} (0:g, r:1, i:2)")
  print(f"  N_row: {n_row}")
  print("==========================================")
  t_start = time.time()
  print(f"  Start GPS time extraction.")

  t0_list, t1_list = [], []
  musec0_list, musec1_list = [], []
  for z in range(nz):
    sec_pre, musec_pre, sec_total_pre = 0, 0, 0
    xmin, xmax = 0, 4
    ymin, ymax = (ny-4), ny
    data = cube[z]
    for i in range(n_row):
      data_temp = data[ymin:ymax, xmin:xmax]
      if band==2:
        # flip ud and lr for i/z band data
        data_temp = np.flipud(data_temp)
        data_temp = np.fliplr(data_temp)
      else:
        # flip ud for g,r band data
        data_temp = np.flipud(data_temp)
      sec, musec = extract_GPStime(data_temp)
      sec_total = sec + musec*1e-6
      if (i!=0) and (sec!=sec_pre):
        # Extract only moving up time.
        t0_list.append(sec_total_pre)
        t1_list.append(sec_total)
        musec0_list.append(musec_pre)
        musec1_list.append(musec)
      # Update values for next iteration.
      ymin -= 4
      ymax -= 4
      sec_pre = sec
      musec_pre = musec
      sec_total_pre = sec_total
  t_end = time.time()
  t_elapsed = t_end - t_start
  print(f"  Finish extraction. elapsed time : {t_elapsed}s")

  # Save as DataFrame.
  d = dict(t0=t0_list, t1=t1_list, musec0=musec0_list, musec1=musec1_list)
  df = pd.DataFrame(d.values(), index=d.keys()).T
  out = filename + "_GPS.csv"
  df.to_csv(out, sep=" ", index=False)


if __name__ == "__main__":
  parser = ap(description="Extract GPS time information")
  parser.add_argument(
    "fits", help='fits obtained with the TriCCS')
  args = parser.parse_args()
  
  main(args)
