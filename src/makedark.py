#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Make a master dark frame from a 3-d fits
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits


def main(args):
    """
    This is the main function called by the `makedark` script.
    
    Note
    ----
    Fits header of first dark frame(s) is used in master dark frame.
  
    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """
    dark_list = []
    for idx,dk in enumerate(args.dark):
        print(f"  Read a dark frame {dk}")
        hdu = fits.open(dk)
        src = hdu[0]
        dark = src.data
  
        # Use header of first dark frame(s)
        if idx == 0:
            hdr_dark = src.header
            hdr_fits = os.path.basename(dk)
  
        # 3d dark frames
        if len(dark.shape)==3:
            nz, _, _ = dark.shape
            if nz == 1:
                print("  Use single 3-d dark fits directory")
                img = dark[0]
                dark_list.append(img)
            else: 
                for z in range(nz):
                    dark_list.append(dark[z])
        # 2d dark frames
        elif len(dark.shape)==2:
            print("  Use 2-d dark fits directly")
            img = dark
            dark_list.append(img)

    # Use median of all 2-d frames as master dark
    print(f"  Dimention of dark frames :{len(dark_list)}")
    masterdark = np.median(dark_list, axis=0)
    print(f"  Shape of a master dark :{masterdark.shape}")
  
    # Create 2-d master dark (Use a header of the first fits!)
    dark = fits.PrimaryHDU(data=masterdark, header=hdr_dark)
    # Add history
    hdr = dark.header
    hdr.add_history(
        f"[makedark] header info. is inherited from {hdr_fits}")
    for idx,dk in enumerate(args.dark):
        hdr.add_history(
        f"[makedark] original dark frames : {idx} {dk}")
  
    if args.out is None:
        out = f"d{os.path.basename(args.dark[0])}"
    else:
        out = args.out
  
    # Write new fits 
    dark.writeto(out, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = ap(description="Make a mater dark frame from a 3-d fits")
    parser.add_argument(
        "dark", type=str, nargs="*",
        help="fits files to compile a dark-frame")
    parser.add_argument(
        "--out", type=str, default=None,
        help="output fits file")
    parser.add_argument(
        "-f", dest="overwrite", action="store_true",
        help="overwrite a fits image")
    args = parser.parse_args()
  
    main(args)
