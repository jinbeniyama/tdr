#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Cut and rotate 2-d fits image.
"""
from argparse import ArgumentParser as ap
import astropy.io.fits as fits
import numpy as np
import os
import datetime 


def main(args):
    """This is the main function called by the `cut_rot.py` script.

    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """
    
    fitsname = os.path.basename(args.fits)
    basename = fitsname.split(".fits")[0]

    # Temporally determine a used filter from file name.
    # Of course, header information is useful as well.
    band = basename.split("TRCS")[1][7]

    xmin, xmax = args.x0, args.x1
    ymin, ymax = args.y0, args.y1

    # Open a 3-d fits
    # Header keywords are optimized for TriCCS
    hdu = fits.open(args.fits)
    hdr = hdu[0].header
    tframe = hdr["TFRAME"]

    data = hdu[0].data
    data = data[ymin:ymax, xmin:xmax]
    # Rotate
    if args.rot:
        data = np.rot90(data)
        rot = " (w/rot90)"
    else:
        rot = " (wo/rot90)"

    hdu[0].data = data
    hdu[0].header.add_history(
        f"[cut_rot] original fits : {fitsname}{rot}")
    out = f"{basename}_cr.fits"
    hdu.writeto(os.path.join(out), overwrite=True)


if __name__ == "__main__":
    parser = ap(
        description="Convert video (3-d fits) to image (2-d fits).")
    parser.add_argument(
        "fits", type=str, 
        help="a reduced 2-d fits")
    parser.add_argument(
        "x0", type=int, help="x0")
    parser.add_argument(
        "x1", type=int, help="x1")
    parser.add_argument(
        "y0", type=int, help="y0")
    parser.add_argument(
        "y1", type=int, help="y1")
    parser.add_argument(
        "--rot", action="store_true", default=False, 
        help="rotate fits")
    args = parser.parse_args()

    main(args)
