#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
wcs check for TriCCS.
after wcsIDcheck_TriCCS.py

"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import pandas as pd
import astropy.io.fits as fits
import subprocess

def main(args):
  """This is the main function called by the `wcsIDshow_TriCCS` script.

  Parameters
  ----------
  args : argparse.Namespace
    Arguments passed from the command-line as defined below.
  """

  flist = args.flist
  N_band = len(flist)
  print(f"  N_band = {N_band}")

  flag_total = np.zeros(N_band)
  for idx,f in enumerate(flist):
      df = pd.read_csv(f, sep=" ")
      if idx==0:
          N_fits = len(df)
          print(f"  N_fits = {N_fits}")
          flag_total = np.zeros(N_fits)
      flag_total = [x + y for x, y in zip(flag_total, df["flag"])]
          
  s_temp = pd.Series(flag_total, name="flag_total")
  df_temp = pd.DataFrame(s_temp)
  # Use df 2
  df = pd.concat([df, df_temp], axis=1)
  N_3 = len(df[df["flag_total"]==3])
  N_2 = len(df[df["flag_total"]==2])
  N_1 = len(df[df["flag_total"]==1])
  N_0 = len(df[df["flag_total"]==0])

  N_all = N_3 + N_2 + N_1 + N_0
  print(f"N_3 : {N_3}")
  print(f"N_2 : {N_2}")
  print(f"N_1 : {N_1}")
  print(f"N_0 : {N_0}")
  print()
  print(f"N_all : {N_all}")

  # Output common fits list
  for idx,f in enumerate(flist):
      df = pd.read_csv(f, sep=" ")
      df = pd.concat([df, df_temp], axis=1)
      df_common = df[df["flag_total"]==3]
      ID_common = df_common["ID"].values.tolist()
      with open(f"idlist_{idx}_common.txt", "w") as f:
          for x in ID_common:
              f.write(f"{x}\n")

if __name__ == "__main__":
  parser = ap(description="Show wcs pasting status.")
  parser.add_argument(
    "flist", type=str, nargs="*",
    help="fitslist with ID and wcs pasting info")
  args = parser.parse_args()
  
  main(args)
