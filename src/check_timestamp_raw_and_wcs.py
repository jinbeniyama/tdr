#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Check timestamp of TRCS00088900.fits.gz and rTRCS00088900ms0001w.fits
keyword is UTC0 in cutted fits
ex) UTC0    = '2021-11-26T09:54:32.400000' / exposure starting date and time

keyword is DATE-OBS and UTC0 in raw fits
ex) DATE-OBS= '2021-11-26'         / Observation start date (yyyy-mm-dd)
ex) UT-STR  = '09:54:32.40'        / UT at exposure start (HH:MM:SS.SS)
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import pandas as pd
import astropy.io.fits as fits
import datetime
from astropy.time import Time
import matplotlib.pyplot as plt



if __name__ == "__main__":
  parser = ap(description="Make a mater flat frame from flat frame(s)")
  parser.add_argument(
    "rlist", type=str, 
    help="fitslist before wcs pasting")
  parser.add_argument(
    "wlist", type=str, 
    help="fitslist after wcs pasting")
  parser.add_argument(
    "--rawdir", type=str, default="raw",
    help="directory for not wcs painted")
  parser.add_argument(
    "--wcsdir", type=str, default="wcs",
    help="directory for wcs painted")
  parser.add_argument(
    "--key_wcs", type=str, default="UTC0",
    help="time keyword for wcs pasted fits")
  parser.add_argument(
    "--out", type=str, default=None,
    help="output png name")
  args = parser.parse_args()
    

  t_raw = []
  rawdir = args.rawdir
  with open(args.rlist, "r") as f:
      lines = f.readlines()
      for line in lines:
          infits = line.split("\n")[0]
          infits = os.path.join(rawdir, infits)
          hdu = fits.open(infits)
          hdr = hdu[0].header
          # DATE-OBS= '2021-11-26'         / Observation start date (yyyy-mm-dd)
          # UT-STR  = '09:54:32.40'        / UT at exposure start (HH:MM:SS.SS)
          t0 = f"{hdr['DATE-OBS']}T{hdr['UT-STR']}"
          t_raw.append(t0)
  t_raw = Time(t_raw, format='isot', scale='utc')
  t_raw_jd = t_raw.jd

  t_wcs = []
  args = parser.parse_args()
    

  t_raw = []
  rawdir = args.rawdir
  with open(args.rlist, "r") as f:
      lines = f.readlines()
      for line in lines:
          infits = line.split("\n")[0]
          infits = os.path.join(rawdir, infits)
          hdu = fits.open(infits)
          hdr = hdu[0].header
          # DATE-OBS= '2021-11-26'         / Observation start date (yyyy-mm-dd)
          # UT-STR  = '09:54:32.40'        / UT at exposure start (HH:MM:SS.SS)
          t0 = f"{hdr['DATE-OBS']}T{hdr['UT-STR']}"
          t_raw.append(t0)
  t_raw = Time(t_raw, format='isot', scale='utc')
  t_raw_jd = t_raw.jd

  t_wcs = []
  wcsdir = args.wcsdir
  with open(args.wlist, "r") as f:
      lines = f.readlines()
      for line in lines:
          infits = line.split("\n")[0]
          infits = os.path.join(wcsdir, infits)
          hdu = fits.open(infits)
          hdr = hdu[0].header
          # ex) UTC0    = '2021-11-26T09:54:32.400000' / exposure starting date and time
          t0 = hdr[args.key_wcs]
          t_wcs.append(t0)
  t_wcs = Time(t_wcs, format='isot', scale='utc')
  t_wcs_jd = t_wcs.jd
   
  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot()
  ax.axes.yaxis.set_visible(False)
  ax.set_xlabel("JD")
  y0 = 0
  y1 = 0.5
  y2 = 1
  # Raw 
  for t in t_raw_jd:
      ax.vlines(t, y0, y1, ls="solid", label="Raw")
  for t in t_wcs_jd:
      ax.vlines(t, y1, y2, ls="dashed", label="Cutted")
  ax.legend()
  if args.out:
      out = args.out
  else:
      out = "checktime.png" 
  fig.savefig(out)
   
