#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Check frame for flat fielding.
Useless?
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits


def main(args):
    """
    This is the main function.
  
    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """
    assert (args.width%2)==0, "Width should be even number"
    width = args.width
    rad = width/2.
    # Searched range x-49 to x+50 in total 100 pix
    x0 = int(args.x - rad)
    x1 = int(args.x + rad)
    y0 = int(args.y - rad)
    y1 = int(args.y + rad)
    print(f"x0, x1 = {x0}, {x1}")
    print(f"y0, y1 = {y0}, {y1}")
    
    # Read flat frame
    fi = os.path.basename(args.flat).split(".")[0]
    hdu_flat = fits.open(args.flat)
    src_flat = hdu_flat[0]
    flat = src_flat.data
    assert len(flat.shape)==2, "Flat should be 2-d fits."
  
  
    # Cut and calculate std after flat
    data_flat = hdu_flat[0].data[y0:y1, x0:x1]
    data_flat_list = data_flat.flatten()
    flat_mean, flat_std = np.mean(data_flat), np.std(data_flat)
    per_flat = flat_std/flat_mean*100.
    
    
    # Plot histograms
    import matplotlib.pyplot as plt  
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    ax1.set_title(fi)
    ax1.set_xlabel("Count [ADU]")
    ax1.set_ylabel("N")
    
    label_flat = (
        f"(Mean, Std) = ({flat_mean:.4f}, {flat_std:.4f})\n"
        f"Std/Mean = {per_flat:.3f} %\n"
        f"(x_c, y_c) = ({x0}-{x1}, {y0}-{y1})"
        )
    print(label_flat)
    bins = 100
    ax1.hist(data_flat_list, bins=bins, label=label_flat)
     
    ax1.legend()
    out = f"flat_hist_{fi}.png"
    plt.savefig(out)
 


if __name__ == "__main__":
    parser = ap(description="Check flat-fielding")
    parser.add_argument(
      "flat", type=str, 
      help="a raw flat frame")
    parser.add_argument(
      "-x", type=float, default=1000,
      help="central x")
    parser.add_argument(
      "-y", type=float, default=564,
      help="central y")
    parser.add_argument(
      "--width", type=float, default=100,
      help="width of pixels")
    args = parser.parse_args()

    main(args)
