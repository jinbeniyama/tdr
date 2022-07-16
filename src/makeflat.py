#!/usr/bin/env python
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
    # Read a master dark frame
    hdu_dark = fits.open(args.dark)
    src_dark = hdu_dark[0]
    dark = src_dark.data
    assert len(dark.shape)==2, "A master dark frame should be 2-d fits."
    print(f"  Read a dark frame {args.dark}")
    
    # Read flat frames
    N_flat = len(args.flat)
    flatcube = []
    for idx,fl in enumerate(args.flat):
      hdu_flat = fits.open(fl)
      src_flat = hdu_flat[0]
      flat = src_flat.data
      if idx==0:
        hdr_flat = src_flat.header
        hdr_fits = os.path.basename(fl)
      
      # 2d flat frame
      if len(flat.shape)==2:
        img = norm2dflat(flat, dark, args.norm)
        # Add temporally cube 
        flatcube.append(img)
        
      # 3d flat frames
      elif len(flat.shape)==3:
        img = norm3dflat(flat, dark, args.norm)
        # Add temporally cube 
        flatcube.append(img)
      print(f"  Read a flat frame {fl}")

    # Create a master flat from a flat cube
    print(f"  Dimention of flat cube :{len(flatcube)}")
    masterflat = np.median(flatcube, axis=0)
    print(f"  Shape of a master flat :{masterflat.shape}")
    # Create 2-d master flat
    flat = fits.PrimaryHDU(data=masterflat, header=hdr_flat)
    # Add history
    hdr = flat.header
    hdr.add_history(
      f"[makeflat] header info. is inherited from {hdr_fits}")
    hdr.add_history(
      f"[makeflat] dark frame : {os.path.basename(args.dark)}")
    for idx,fl in enumerate(args.flat):
      hdr.add_history(
        f"[makeflat] flat frame {idx+1} : {os.path.basename(fl)}")
    hdr.add_history(
      f"[makeflat] normalization : {args.norm}")

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
      help="a raw flat frame")
    parser.add_argument(
      "--dark", type=str, required=True,
      help="a master dark frame")
    parser.add_argument(
      "--norm", type=str, default="mean",
      help="type of normalization (mean or median)")
    parser.add_argument(
      "--out", type=str, default=None,
      help="output fits file")
    parser.add_argument(
      "-f", dest="overwrite", action="store_true",
      help="overwrite a fits image")
    args = parser.parse_args()
    
    main(args)
