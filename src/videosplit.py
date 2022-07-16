#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Split 3-d fits to n 3-d fits.
This script is useful, for example, to avoid memory error.
(We recommend to check the memory before observations and set appropreate 
NAXIS3.)
"""
from argparse import ArgumentParser as ap
import astropy.io.fits as fits
import numpy as np
import os
import datetime 


def main(args):
    """This is the main function called by the `videosplit.py` script.

    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """
    
    # Create a directory to save output fits
    fitsname = os.path.basename(args.fits)
    basename = fitsname.split(".fits")[0]


    # Open a 3-d fits
    # Header keywords are optimized for TriCCS
    hdu = fits.open(args.fits)
    hdr = hdu[0].header
    tframe = hdr["TFRAME"]

    # Check this part carefully !!! ===========================================
    # Obtain an exposure starting time 
    try:
        # For old one (UTC is exposure 'starting' time, not 'ending') -2021.08
        # (K. Matsubayashi, private communications in 2022.07)
        # header log (ending) is wrong !!!
        # UTC     = '2021-03-08T16:11:32.555790' / exposure ending date and time
        t0 = hdr["UTC"] 
        t0_dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%f")
        kwd_time = "UTC"
    except:
        # For new one 2021.08-
        # GEXP-STR= '2021-10-27T20:10:50.371615' / GPS time at exposure start
        t0 = hdr["GEXP-STR"] 
        t0_dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%f")
        kwd_time = "GEXP-STR"

        # Some headers below are not recomended !
        # DATE-OBS= '2021-10-15'         / Observation start date (yyyy-mm-dd)
        # UT-STR  = '14:40:07.40'        / UT at exposure start (HH:MM:SS.SS)
        # UT_STR is not GPS time, but telescope time. 
        # (~ 1 sec differs from GPS time)
    # =========================================================================

    data_temp = hdu[0].data


    # Mask edge and split for 3-d fits
    if len(data_temp.shape)==3:
        nz, ny, nx = data_temp.shape
        N_fits = args.n
        fits_per_video = int(nz/args.n)
        print(f"  Split to {N_fits} fits")
        print(f"  nz : {nz} to {fits_per_video}")
        hdu[0].header.add_history(
            f"[videosplit] original fits : {fitsname}")
        for n in range(N_fits):
            # epalsed time
            t_elapse = tframe*n*fits_per_video
            # Starting time of exposure
            t0_dt_temp = t0_dt + datetime.timedelta(seconds=t_elapse)
            # Central time of exposure
            t0_temp = datetime.datetime.strftime(t0_dt_temp, "%Y-%m-%dT%H:%M:%S.%f")

            # Save as kwd_time     
            # UTC and UTC0 is calculated by video2image.py
            hdr[kwd_time] = (t0_temp, "exposure starting date and time")

            # Use sensitive pixels
            temp = data_temp[n*fits_per_video:(n+1)*fits_per_video, :, :]
            hdu[0].data = temp
            # s: split
            out = f"{basename}_s{n+1:03d}.fits"
            hdu.writeto(os.path.join(out), overwrite=True)


if __name__ == "__main__":
    parser = ap(
        description="Convert video (3-d fits) to n 3-d fits.")
    parser.add_argument(
        "fits", type=str, 
        help="a raw 3-d fits")
    parser.add_argument(
        "n", type=int, 
        help="split number")
    args = parser.parse_args()

    main(args)
