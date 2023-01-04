#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Copy time for fits with broken time.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits


def main(args):
    """This is the main function called by the `copytime.py` script.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passed from the command-line as defined below.
    """

    # Read fits with broken time
    hdu1 = fits.open(args.fits_broken)
    hdr1 = hdu1[0].header

    # Read fits with correct time
    hdu2 = fits.open(args.fits_correct)
    hdr2 = hdu2[0].header

    # Check whether the time of fits_broken is broken
    t1 = hdr1["GEXP-STR"]
    t2 = hdr2["GEXP-STR"]
    assert t1 == "1970-01-01T00:00:00", "Possibly not broken?"

    # Print time info.
    print(f"Overwrite\n{t1}\nwith\n{t2}")

    # Overwrite GEXP-STR
    hdr1["GEXP-STR"] = hdr2["GEXP-STR"]
    
    # Add history
    hdr1.add_history(
      f"[copytime] original fits is {os.path.basename(args.fits_broken)}")
    hdr1.add_history(
      f"[copytime] time is copied with {os.path.basename(args.fits_correct)}")

    if args.out is None:
      out = f"t{os.path.basename(args.fits_broken)}"
    else:
      out = args.out

    # Write new fits in current directory
    hdu1.writeto(out, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = ap(description="Copy time info.")
    parser.add_argument(
      "fits_broken", type=str,
      help="a raw or reduced fits with broken time")
    parser.add_argument(
      "fits_correct", type=str,
      help="a raw or reduced fits with correct time")
    parser.add_argument(
      "--out", type=str, default=None,
      help="output fits file")
    parser.add_argument(
      "-f", dest="overwrite", action="store_true",
      help="overwrite a fits image")
    args = parser.parse_args()

    main(args)
