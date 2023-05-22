#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Shift and stack multiple 2d fits files.
Three modes are available:
1. velocity (stack simply with pixel coordinates, useful only for quick look maybe)
    wcs of output fits is the same with the first fits
2. sidereal (stack with wcs coordinates to raise SNR of 'stars')
    wcs of output fits is the same with the first fits
3. nonsidereal (stack with certain velocity to raise SNR of a moving object)
    wcs of output fits is the same with the first fits
    OBS-TIME is necessary for 3.
    Please check OBS-TIME keyword.

Params: location, velocity, time
Return: stacked fits, 3d fits 

Note:
    Exposure times of all input fits should be the same to calculate
    UTC0 and UTC correctly.
"""
from astropy.io import fits as fits
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
from astropy.time import Time
import datetime

from tdr import (
    instinfo, obtain_fitstime, shift_w_velocity, 
    shift_sidereal, shift_nonsidereal)


if __name__ == "__main__":
    parser = ap(description="Shift and stack 2d fits")
    parser.add_argument(
        "mode", choices=["velocity", "sidereal", "nonsidereal"],
        help="mode of shift")
    parser.add_argument(
        "--inst", 
        choices=["ishigaki", "okayama", "akeno", "saitama", "triccs"], 
        default="triccs",
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
        help="RA velocity in arcsec/hr (sky motion)")
    parser.add_argument(
        "--v_dec", type=float, default=0,
        help="DEC velocity in arcsec/hr (sky motion)")
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
    print(f"  Number of fits = {len(flist)}")


    # Extract time to update it
    key_utc = "UTC"
    jd_list = []
    for fi in flist:
        hdu = fits.open(fi)
        hdr = hdu[0].header
        utc = hdr[key_utc]
        t = Time(str(utc), format='isot', scale='utc')
        # Convert to JD
        jd = t.jd
        jd_list.append(jd)
    # Calculate median jd
    jd_median = np.median(jd_list)
    jd_mean = np.mean(jd_list)
    assert jd_median == jd_mean, "Update the code."
    


    # Open the first fits
    fi0 = flist[0]
    hdu = fits.open(fi0)
    # Update header
    hdu[0].header.add_history(f"[shift2d] Mode {args.mode}")
    for fi in flist:
        hdu[0].header.add_history(f"[shift2d] Created from {fi}")

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
        str_mode = "vel"

    # Shift with wcs
    elif args.mode == "sidereal":
        img_shift = shift_sidereal(flist)
        str_mode = "sid"

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
            f"[shift2d] Target {args.target}")
        str_mode = "nonsid"

    # Update data
    hdu[0].data = img_shift
    Nfits = len(flist)
    hdu[0].header.add_history(f"[shift2d] Number of stacked images {Nfits}")

    # Update Time
    hdu[0].header["JD_MED"] = (jd_median, "Median JD")
    hdu[0].header["JD_MEAN"] = (jd_mean, "Mean JD")
    # Convert to UTC
    t_median = Time(str(jd_median), format='jd', scale='utc')
    utc_median = t_median.datetime
    utc_median = datetime.datetime.strftime(utc_median, "%Y-%m-%dT%H:%M:%S.%f")
    t_mean = Time(str(jd_mean), format='jd', scale='utc')
    utc_mean = t_mean.datetime
    utc_mean = datetime.datetime.strftime(utc_mean, "%Y-%m-%dT%H:%M:%S.%f")
    hdu[0].header["UTC_MED"] = (utc_median, "Median UTC")
    hdu[0].header["UTC_MEAN"] = (utc_mean, "Mean UTC")
    # Use median
    # Note: The time between frames should be the same.
    #       Otherwise care should be taken when using the inherited time here.
    hdu[0].header[key_utc] = utc_median

    # Update single exposure time (EXPTIME1) and time between frames (TFRAME)
    key_exp1, key_tframe = "EXPTIME1", "TFRAME"
    t_exp = hdu[0].header[key_exp1]
    t_exp_total = Nfits * t_exp
    hdu[0].header[key_exp1] = t_exp_total
    hdu[0].header[key_tframe] = t_exp_total


    if args.out:
        out = args.out
    else:
        filename = os.path.basename(fi0).split(".")[0]
        out = f"{str_mode}_{filename}.fits"
    out = os.path.join(outdir, out)
    hdu.writeto(out, overwrite=True)
