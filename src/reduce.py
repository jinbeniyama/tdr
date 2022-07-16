#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Do standard reduction for an object frame.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits


def main(args):
    """This is the main function called by the `reduce` script.

    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """

    # Create a directory to save output fits
    outdir = "reduce"
    if os.path.isdir(outdir):
        pass
    else:
        os.makedirs(outdir)

    # Read object frame
    hdu_obj = fits.open(args.obj)
    hdr = hdu_obj[0].header

    # Read dark frame
    hdu_dark = fits.open(args.dark)
    src_dark = hdu_dark[0]
    dark = src_dark.data
    assert len(dark.shape)==2, "Dark should be 2-d fits."

    # Read flat frame
    hdu_flat = fits.open(args.flat)
    src_flat = hdu_flat[0]
    flat = src_flat.data
    assert len(flat.shape)==2, "Flat should be 2-d fits."

    # Do dark subtraction
    print(f"median count (before dark subtraction) {np.median(hdu_obj[0].data):.1f}")
    hdu_obj[0].data = hdu_obj[0].data - dark
    print(f"median count (after dark subtraction) {np.median(hdu_obj[0].data):.1f}")

    # Do flat-field correction
    hdu_obj[0].data = hdu_obj[0].data/flat

    # Add history
    hdr.add_history(
      f"[reduce] created from (dark) {os.path.basename(args.dark)}")
    hdr.add_history(
      f"[reduce] created from (flat) {os.path.basename(args.flat)}")

    if args.out is None:
      out = f"r{os.path.basename(args.obj)}"
    else:
      out = args.out

    # Write new fits 
    hdu_obj.writeto(os.path.join(outdir, out), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = ap(description="Do standard reduction for an object frame")
    parser.add_argument(
      "--obj", type=str, required=True,
      help="a raw fits image")
    parser.add_argument(
      "--dark", type=str, required=True,
      help="a raw dark frame")
    parser.add_argument(
      "--flat", type=str, required=True,
      help="a raw flat frame")
    parser.add_argument(
      "--out", type=str, default=None,
      help="output fits file")
    parser.add_argument(
      "-f", dest="overwrite", action="store_true",
      help="overwrite a fits image")
    args = parser.parse_args()

    main(args)
