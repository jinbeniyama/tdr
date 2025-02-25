#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Check two frames such as darks and flats.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt  


def main(args):
    """
    This is the main function.
  
    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """
    # Read flat frame 1
    hdu_flat1 = fits.open(args.flat1)
    src_flat1 = hdu_flat1[0]
    flat1 = src_flat1.data
    assert len(flat1.shape)==2, "Flat should be 2-d fits."

    # Read flat frame 2
    hdu_flat2 = fits.open(args.flat2)
    src_flat2 = hdu_flat2[0]
    flat2 = src_flat2.data
    assert len(flat2.shape)==2, "Flat should be 2-d fits."
  
    # Cut and calculate std after flat
    div = flat1/flat2
    # Sub
    sub = flat1 - flat2
    sub_min, sub_max, sub_std = np.min(sub), np.max(sub), np.std(sub)
    print(sub_min, sub_max, sub_std)
    sub_vmin, sub_vmax = -5*sub_std, 5*sub_std

    # Plot histograms
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_axes([0.15, 0.55, 0.75, 0.40])
    ax2 = fig.add_axes([0.15, 0.05, 0.75, 0.40])

    label = "frame1/frame2"
    im1 = ax1.imshow(div, vmin=args.vmin, vmax=args.vmax, cmap=args.cmap, label=label)
    fig.colorbar(im1, label="Normalized count")

    label = "frame1 - frame2"
    im2 = ax2.imshow(sub, cmap=args.cmap, label=label, vmin=sub_vmin, vmax=sub_vmax)
    fig.colorbar(im2, label="Subtracted count")
     

    ax1.legend()
    ax2.legend()

    plt.savefig(args.out)
 


if __name__ == "__main__":
    parser = ap(description="Compare frames.")
    parser.add_argument(
        "flat1", type=str, 
        help="a raw flat frame")
    parser.add_argument(
        "flat2", type=str, 
        help="a raw flat frame")
    parser.add_argument(
        "--vmin", type=float, default=0.99,
        help="Minimum value")
    parser.add_argument(
        "--vmax", type=float, default=1.01,
        help="Maxmimum value")
    parser.add_argument(
        "--cmap", type=str, default="coolwarm",
        help="Color map")
    parser.add_argument(
        "--out", type=str, default="frame_comp.jpg",
        help="Output file")
    args = parser.parse_args()

    main(args)
