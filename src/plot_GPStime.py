#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot extracted GPS time information.
"""
import numpy as np
import pandas as pd
from argparse import ArgumentParser as ap
import matplotlib.pyplot as plt


def main(args):
  """Main function.
  """
  df_g = pd.read_csv(args.csv, sep=" ")
  
  fig = plt.figure(figsize=(10, 5))
  ax_g = fig.add_axes([0.2, 0.2, 0.5, 0.7])
  ax_g.set_xlabel("move up event")
  ax_g.set_ylabel("move up time [micro second]")
  title = f"{args.title} N={len(df_g)}"
  ax_g.set_title(title)

  # Use g-band
  frame = [i+1 for i in range(len(df_g))]
  # before moving up
  ax_g.plot(
    frame, df_g["musec0"], c="blue")
  ax_g.scatter(
    frame, df_g["musec0"], c="blue", label="before moving up")
  # After moving up
  ax_g.plot(
    frame, df_g["musec1"], c="red")
  ax_g.scatter(
    frame, df_g["musec1"], c="red", label="after moving up")
  # Mean 
  df_g["musec_mean"] = (df_g["musec0"] + df_g["musec1"])/2.0
  ax_g.plot(
    frame, df_g["musec_mean"], c="green")
  ax_g.scatter(
    frame, df_g["musec_mean"], c="green", label="(before+after)/2.")
  # Average and standard deviation
  average_g = np.mean(df_g["musec_mean"])
  average_g = np.round(average_g, 2)
  std_g = np.std(df_g["musec_mean"])
  std_g = np.round(std_g, 2)
  xmin, xmax = ax_g.get_xlim()
  ax_g.hlines(
    average_g, xmin, xmax, lw=1, color="green", linestyle="dashed", 
    label=f"Average: {average_g} us\nStandard Deviation: {std_g} us")
  plt.fill_between([xmin, xmax], average_g+std_g, average_g-std_g, 
    facecolor='green', alpha=0.3)
  ax_g.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0,
            markerscale=2)
  plt.savefig(args.out)


if __name__ == "__main__":
  parser = ap(description="Plot extracted GPS time information")
  parser.add_argument(
    "csv", type=str, 
    help="a csv of a GPS time")
  parser.add_argument(
    "-out", type=str, default="GPStime.png", 
    help="a name of an output image")
  parser.add_argument(
    "-title", type=str, default=None, 
    help="figure title")
  args = parser.parse_args()
   
  main(args)
