#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Shift and stack multiple 2d fits files.
Three modes are available:
1. velocity (stack simply)
2. sidereal (raise SNR of stars)
3. nonsidereal (raise SNR of a target for non-sidereal tracking data)

OBS-TIME is necessary for 3. nonsidereal stacking.
Please check OBS-TIME keyword.
Params: location, velocity, time
Return: stacked fits, 3d fits 
"""
from astropy.io import fits as fits
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd

from tdr import (
    instinfo, obtain_fitstime, shift_w_velocity, 
    shift_sidereal, shift_nonsidereal)


if __name__ == "__main__":
    parser = ap(description="Shift and stack 2d fits")
    parser.add_argument(
        "mode", choices=["velocity", "sidereal", "nonsidereal"],
        help="mode of shift")
    parser.add_argument(
        "inst", choices=["ishigaki", "okayama", "akeno", "saitama", "triccs"],
        help="instrument")
    parser.add_argument(
        "--fits", type=str, nargs="*", 
        help="fits to be shifted")
    parser.add_argument(
        "--flist", type=str, 
        help="fits list text file")
    # For velocity
    parser.add_argument(
        "--v_ra", type=float, default=0,
        help="RA velocity in arcsec/ho (sky motion)")
    parser.add_argument(
        "--v_dec", type=float, default=0,
        help="DEC velocity in arcsec/ho (sky motion)")
    parser.add_argument(
        "--second", action="store_true", default=None,
        help="velocity in arcsec/s")
    # For nonsidereal
    parser.add_argument(
        "--target", type=str, default=None, 
        help="Target name like '1', 'Ceres' or '2015 RN35'")
    parser.add_argument(
        "--out", type=str, default=None,
        help="output fits name")
    parser.add_argument(
        "--outdir", type=str, default=".",
        help="output directory")
    args = parser.parse_args()

    # Create output directory
    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
    

    # Create a list as flist
    if args.fits:
        flist = args.fits
    elif args.flist:
        with open (args.flist, "r") as f:
            f = f.readlines()
            flist = []
            for line in f:
                fi = line.strip("\n")
                flist.append(fi)
    print(f"Number of fits = {len(flist)}")

    # Open the first fits
    fi0 = flist[0]
    hdu = fits.open(fi0)
    # Update header
    hdu[0].header.add_history("[shift2d] mode {args.mode}")
    for fi in flist:
        hdu[0].header.add_history("[shift2d] created from {fi}")

    # Shift with velocity
    if args.mode == "velocity":
        # Extract header keyword and pixel scale
        hdr_kwd, p_scale = instinfo(args.inst)
        print(f"Instrument  : {args.inst}")
        print(f"Header info : {hdr_kwd}")
        print(f"Pixel scale : {p_scale} arcsec/pixel")

        # Velocity in arcsec/hour (default format of JPL/HORIZONS)
        v_ra, v_dec = args.v_ra, args.v_dec

        # Velocity in arcsec/s 
        if args.second:
            v_ra, v_dec = v_ra/3600., v_dec/3600.
        else:
            pass

        # Velocity in pixel/s 
        v_ra, v_dec = v_ra/p_scale, v_dec/p_scale
        
        # Shift
        img_shift = shift_w_velocity(flist, hdr_kwd, v_ra, v_dec)
        # Update header
        hdu[0].header.add_history(
            "[shift2d] (v_ra, v_dec) = ({v_ra}, {v_dec}) pixel/s")

    # Shift with wcs
    elif args.mode == "sidereal":
        img_shift = shift_sidereal(flist)

    # Shift for a certain target
    elif args.mode == "nonsidereal":
        # Extract header keyword and pixel scale
        hdr_kwd, p_scale, loc = instinfo(args.inst)
        print(f"Instrument  : {args.inst}")
        print(f"Header info : {hdr_kwd}")
        print(f"Pixel scale : {p_scale} arcsec/pixel")
        img_shift = shift_nonsidereal(flist, hdr_kwd, args.target, loc)
        # Update header
        hdu[0].header.add_history(
            "[shift2d] target {args.target}")

    # Update data
    hdu[0].data = img_shift
    if args.out:
        out = args.out
    else:
        filename = os.path.basename(fi0).split(".")[0]
        out = f"sa{filename}.fits"
    out = os.path.join(outdir, out)
    hdu.writeto(out, overwrite=True)
