#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Stack 3-d fits by median cube and create a 2-d fits.
This script should be used after reduction (reduce.py).
Header keywords are optimized for Seimei/TriCCS.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits


hdr_kwd = {
  "tframe": "TFRAME", 
  "exp1": "EXPTIME1",
  "fps": "DATA-FPS"
  }


__naxis3_keywords = (
  'NAXIS3', 'CTYPE3', 'CRPIX3', 'CRVAL3', 'CUNIT3',
  'CD1_3', 'CD2_3', 'CD3_3', 'CD3_2', 'CD3_1',
)


def stack3dfits(cube, type):
    """Stack 3-di fits file.
    
    Parameters
    ----------
    cube : 3-d array-like
      3-dimentional array to be stacked
    type : str
      stacking type ('max', 'min', 'mean' or 'median')
  
    Return
    ------
    img : 2-d array-like
      compressed 2-dimentional array
    """
  
    if args.type=="max":
      img = np.max(cube, axis=0)
    elif args.type=="min":
      img = np.min(cube, axis=0)
    elif args.type=="mean":
      img = np.mean(cube, axis=0)
    elif args.type=="median":
      img = np.median(cube, axis=0)
  
    return img


def main(args=None):
    """This is the main function called by the `stackfits` script.
  
    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """
    parser = ap(description="Stack 3-d fits")
    parser.add_argument(
      "fits", type=str, 
      help="a raw 3-d fits")
    parser.add_argument(
      "type", type=str, choices=["max", "min", "mean", "median"], 
      help="stacking type")
    parser.add_argument(
      "-f", dest="overwrite", action="store_true",
      help="overwrite a fits image")
    args = parser.parse_args()
  
    # Extract filename
    filename = os.path.basename(args.fits)
 
    hdu = fits.open(args.fits)
    hdr = hdu[0].header
    cube = hdu[0].data
    assert len(cube.shape)==3, "Input should be 3-d fits."
  
    nz, ny, nx = cube.shape
  
    # Extract and update header keywords
    tframe = hdr[hdr_kwd["tframe"]]*nz
    exp1 = hdr[hdr_kwd["exp1"]]*nz
    fps = hdr[hdr_kwd["fps"]]/nz
     
    # Create 2-d image
    img = stack3dfits(cube, args.type)
  
    # Update header 
    # Remove useless keywords
    for key in __naxis3_keywords: hdr.remove(key, ignore_missing=True)
    hdr.set("NAXIS", 2)
    hdr.add_history(
      f"[stack] created from {filename}")
    hdr.add_history(
      f"[stack] type {args.type}")

    # Update header keywords
    key_tframe0 = "TFRAME"
    key_texp0   = "EXPTIME1"
    # Previously this was DATE-FPS (<< 2024-10)
    key_fps0    = "DAT-FPS"
    hdr[key_tframe0] = tframe
    hdr[key_texp0] = exp1
    hdr[key_fps0] = fps
  
    # Write new fits 
    hdu[0].data = img
    out = f"{args.type}{filename}"
    hdu.writeto(out, overwrite=args.overwrite)


if __name__ == '__main__':
    parser = ap(description="Stack 3-d fits")
    parser.add_argument(
      "fits", type=str, 
      help="a raw 3-d fits")
    parser.add_argument(
      "type", type=str, choices=["max", "min", "mean", "median"], 
      help="stacking type")
    parser.add_argument(
      "-f", dest="overwrite", action="store_true",
      help="overwrite a fits image")
    args = parser.parse_args()
   
  
    main(args)
