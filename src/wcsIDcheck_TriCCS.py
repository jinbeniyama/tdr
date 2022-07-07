#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
wcs check for TriCCS.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import pandas as pd
import astropy.io.fits as fits
import subprocess

def main(args):
  """This is the main function called by the `makeflat` script.

  Parameters
  ----------
  args : argparse.Namespace
    Arguments passed from the command-line as defined below.
  """

  rid0 = args.rid0
  rid1 = args.rid1
  wid0 = args.wid0
  wid1 = args.wid1
  rdir = args.rdir
  wdir = args.wdir
  rlist = args.rlist
  wlist = args.wlist


  # Obtain fits ID list before wcs pasting
  # Assume handling cutted fits
  idlist_r = []
  for x in rlist:
      id_r = x[rid0:rid1+1]
      print(f"File name    : {x}")
      print(f"extracted id : {id_r}")
      idlist_r.append(id_r)

  idlist_w = []
  for x in wlist:
      id_w = x[wid0:wid1+1]
      print(f"File name    : {x}")
      print(f"extracted id : {id_w}")
      idlist_w.append(id_w)

  #print(idlist_r)
  #print(idlist_w)

  exist_or_not = [1 if x in idlist_w else 0 for x in idlist_r]
  df = pd.DataFrame({"ID":idlist_r, "flag":exist_or_not})
  # ID, 0/1 DataFrame
  #df.to_csv(f"idlist_{args.band}_N{len(df)}.csv", sep=" ", index=False)
  df.to_csv(f"idlist_{args.band}.csv", sep=" ", index=False)


if __name__ == "__main__":
  parser = ap(description="Make a mater flat frame from flat frame(s)")
  parser.add_argument(
    "--rlist", type=str, nargs="*",
    help="fitslist before wcs pasting")
  parser.add_argument(
    "--wlist", type=str, nargs="*",
    help="fitslist before wcs pasting")
  parser.add_argument(
    "--rid0", type=int,
    help="id start for fits before wcs pasting")
  parser.add_argument(
    "--rid1", type=int,
    help="id end for fits before wcs pasting")
  parser.add_argument(
    "--wid0", type=int,
    help="id start for fits after wcs pasting")
  parser.add_argument(
    "--wid1", type=int,
    help="id end for fits after wcs pasting")
  parser.add_argument(
    "--rdir", type=str, default="reduce",
    help="reduced fits directory")
  parser.add_argument(
    "--wdir", type=str, default="wcs",
    help="wcs pasted fits directory")
  parser.add_argument(
    "--band", type=str, default="0",
    help="band")
  args = parser.parse_args()
  
  main(args)
