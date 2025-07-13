#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Make a master flat frame from a 3-d fits
2 types of normalization, median and mean, are implemented.

Various flat frames (i.e. field of view, the number of frames) are acceptible.
But all flat frames should have the same exposure time
since unique dark frame is permitted as an input argument.

The header info. of the 1st flat frame(s) is inherited.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits


def norm2dflat(flat, dark, norm):
    """Return a normalized flat frame from a flat frame.
  
    Parameters
    ----------
    flat : np.ndarray
      a 2-d flat frame
    dark : np.ndarray
      a 2-d dark frame
    norm : str
      normalization method (mean or median)
  
    Return
    ------
    img : np.ndarray
      a normalized 2-d flat frame
    """
    img = flat - dark
    if norm=="mean":
      img = img/np.mean(img)
    elif norm=="median":
      img = img/np.median(img)
    return img


def norm3dflat(flat, dark, norm):
    """Return a normalized flat frame from flat frames.

    Parameters
    ----------
    flat : np.ndarray
      a 3-d flat frames
    dark : np.ndarray
      a 2-d dark frame
    norm : str
      normalization method (mean or median)

    Return
    ------
    img : np.ndarray
      a normalized 2-d flat frame
    """
    nz, ny, nx = flat.shape
    cube = []
    for z in range(nz):
      img = flat[z] - dark
      if norm=="mean":
        img = img/np.mean(img)
      elif norm=="median":
        img = img/np.median(img)
      cube.append(img)
    return np.median(cube, axis=0)


def main(args):
    """This is the main function called by the `makeflat` script.

    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """
    # Exposure time of a single frame
    key_exp1 = "EXPTIME1"

    # Read master dark frame(s) with exposure tiems
    d_dark = dict()
    for d in args.dark:
        hdu_dark = fits.open(d)
        src_dark = hdu_dark[0]
        dark = src_dark.data
        assert len(dark.shape)==2, "A master dark frame should be 2-d fits."
        hdr_dark = src_dark.header
        texp = hdr_dark[key_exp1]
        print(f"  Read a dark frame {d} with exposure time of {texp} s")
        d_dark[texp] = dark
    
    # Read flat frames
    flat_list = []
    for idx,fl in enumerate(args.flat):
        print(f"  Read a flat frame {fl}")
        hdu_flat = fits.open(fl)
        src_flat = hdu_flat[0]
        # Exposure time
        hdr = src_flat.header
        texp = hdr[key_exp1]
        flat = src_flat.data
        # Save header and original fits file 
        if idx==0:
            hdr_flat = src_flat.header
            fits0 = os.path.basename(fl)
         
        # Set dark 
        dark_temp = d_dark[texp]

        # 2d flat frame
        if len(flat.shape)==2:
            flat_temp = flat - dark_temp
            flat_temp = flat_temp/np.median(flat_temp)
            flat_list.append(flat_temp)
        # 3d flat frames
        elif len(flat.shape)==3:
            nz = flat.shape[0]
            for z in range(nz):
                flat_temp = flat[z] - dark_temp
                flat_temp = flat_temp/np.median(flat_temp)
                flat_list.append(flat_temp)

    # Create a master flat from a flat cube
    print(f"  Dimention of flat list :{len(flat_list)}")
    masterflat = np.median(flat_list, axis=0)
    print(f"  Shape of a master flat :{masterflat.shape}")
    # Create 2-d master flat
    flat = fits.PrimaryHDU(data=masterflat, header=hdr_flat)
    # Add history
    hdr = flat.header
    hdr.add_history(
      f"[makeflat] header info. is inherited from {fits0}")
    for idx,d in enumerate(args.dark):
        hdr.add_history(
            f"[makeflat] dark frame {idx+1} : {os.path.basename(d)}")
        hdr.add_history(
            f"[makeflat] flat frame {idx+1} : {os.path.basename(fl)}")
    hdr.add_history(
      f"[makeflat] normalization : median")

    if args.out is None:
      out = f"f{os.path.basename(args.flat[0])}"
    else:
      out = args.out

    # Write new fits 
    flat.writeto(out, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = ap(description="Make a mater flat frame from flat frame(s)")
    parser.add_argument(
      "--flat", type=str, nargs="*", required=True,
      help="raw flat frames")
    parser.add_argument(
      "--dark", type=str, nargs="*", required=True,
      help="master dark frames")
    parser.add_argument(
      "--out", type=str, default=None,
      help="output fits file")
    parser.add_argument(
      "-f", dest="overwrite", action="store_true",
      help="overwrite a fits image")
    args = parser.parse_args()
    
    main(args)
