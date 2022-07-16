#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Split 3-d fits to 2-d fits while cutting non-sensitive pixels 
(and masking edge region).

1. Cut non-sensitive resions
from 
Toral pixels     2220 x 1360
to 
Sensitive pixels 2160 x 1280

2. Mask edge region (optional)
Mask not well corrected pixels 
(selected from 2021.Feb and Mar engineering data by J.BENIYAMA)
x 0:10
y 0:10
x 0:150 y 0:150
x 0:150 y 1130:1280
x 2010:2160 y 0:150
x 2010:2160 y 1130:1280
"""
from argparse import ArgumentParser as ap
import astropy.io.fits as fits
import numpy as np
import os
import datetime 

__naxis3_keywords = (
    'NAXIS3', 'CTYPE3', 'CRPIX3', 'CRVAL3', 'CUNIT3',
    'CD1_3', 'CD2_3', 'CD3_3', 'CD3_2', 'CD3_1',
)


def TriCCSmask():
    """
    Return boolian mask to be masked. 
    This is easy way to paste correct wcs while suppessing 
    false detections from not well corrected pixels. (in J.BENIYAMA's knowledge)
    """
    # Total sensitive pixels
    nx, ny = 2160, 1280
    arr_temp = np.zeros(shape=(ny, nx))
    arr_temp[0:10,:] = 1
    arr_temp[ny-10:ny,:] = 1
    arr_temp[:,0:10] = 1
    arr_temp[:,nx-10:nx] = 1
    arr_temp[0:150, 0:150] = 1
    arr_temp[0:150, nx-150:nx] = 1
    arr_temp[ny-150:ny, 0:150] = 1
    arr_temp[ny-150:ny, nx-150:nx] = 1
    return arr_temp


def stack3dfits(cube, stype):
    """Stack 3-di fits file.
    
    Parameters
    ----------
    cube : 3-d array-like
      3-dimentional array to be stacked
    stype : str
      stacking type ('max', 'min', 'mean' or 'median')
  
    Return
    ------
    img : 2-d array-like
      compressed 2-dimentional array
    """
  
    if args.stype=="max":
      img = np.max(cube, axis=0)
    elif args.stype=="min":
      img = np.min(cube, axis=0)
    elif args.stype=="mean":
      img = np.mean(cube, axis=0)
    elif args.stype=="median":
      img = np.median(cube, axis=0)
  
    return img


def main(args):
    """This is the main function called by the `video2image` script.

    Parameters
    ----------
    args : argparse.Namespace
      Arguments passed from the command-line as defined below.
    """
    
    # Create a directory to save output fits
    outdir = "cutedge"
    if os.path.isdir(outdir):
        pass
    else:
        os.makedirs(outdir)
    
    fitsname = os.path.basename(args.fits)
    basename = fitsname.split(".fits")[0]

    # Temporally determine a used filter from file name.
    # Of course, header information is useful as well.
    band = basename.split("TRCS")[1][7]

    # Select sensitive pixels
    if band=="2":
        xmin, xmax = 1, 2160
        ymin, ymax = 1, 1280

    elif (band=="0") or (band=="1"):
        xmin, xmax = 61, 2220
        ymin, ymax = 1, 1280

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
    except:
        # For new one 2021.08-
        # GEXP-STR= '2021-10-27T20:10:50.371615' / GPS time at exposure start
        t0 = hdr["GEXP-STR"] 
        t0_dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%f")

        # Some headers below are not recomended !
        # DATE-OBS= '2021-10-15'         / Observation start date (yyyy-mm-dd)
        # UT-STR  = '14:40:07.40'        / UT at exposure start (HH:MM:SS.SS)
        # UT_STR is not GPS time, but telescope time. 
        # (~ 1 sec differs from GPS time)
    # =========================================================================

    data_temp = hdu[0].data
    mask = TriCCSmask()

    # Mask edge and split for 3-d fits
    if len(data_temp.shape)==3:
        nz, ny, nx = data_temp.shape
        # Remove useless header keywords
        hdr.set("NAXIS", 2)
        for key in __naxis3_keywords: hdr.remove(key, ignore_missing=True)
        hdu[0].header.add_history(
            f"[video2image] original fits : {fitsname}")

        
        # With stacking
        if args.stack:
            N_stack = int(args.stack)
            N_fits = int(nz/N_stack)
            print(f"  Split to {N_fits} fits")
            hdu[0].header.add_history(
                f"[video2image] stacking number : {args.stack}")
            hdu[0].header.add_history(
                f"[video2image] stacking type : {args.stype}")
            for i in range(N_fits):
                # Starting time of exposure
                t0_dt_temp = t0_dt + datetime.timedelta(seconds=tframe*i*N_stack)
                # Central time of exposure
                t_mid_dt_temp = t0_dt_temp + datetime.timedelta(seconds=tframe*N_stack*0.5)
                t0_temp = datetime.datetime.strftime(t0_dt_temp, "%Y-%m-%dT%H:%M:%S.%f")
                t_mid_temp = datetime.datetime.strftime(t_mid_dt_temp, "%Y-%m-%dT%H:%M:%S.%f")
                hdr["UTC0"] = (t0_temp, "exposure starting date and time")
                hdr["UTC"] = (t_mid_temp, "central exposure date and time")
                # Update TFRAME
                hdr["TFRAME"] = (tframe*N_stack, "frame interval in seconds")
                hdr["TFRAME0"] = (tframe, "original frame interval in seconds")

                # Extract (N_stack) 3-d fits removing unsensitive pixels
                temp = data_temp[i*N_stack:(i+1)*N_stack, (ymin-1):ymax, (xmin-1):xmax]
                # Stacking 
                temp = stack3dfits(temp, stype=args.stype)
                if args.mask:
                    # Mask not well corrected pixels
                    temp = np.where(mask==1, 1.0, temp)
                hdu[0].data = temp
                out = f"{basename}s{N_stack}ms{i+1:04d}.fits"
                hdu.writeto(os.path.join(outdir, out), overwrite=True)


        else:
            # Without stacking
            print(f"  Split to {nz} fits")
            for i in range(nz):
                # Starting time of exposure
                t0_dt_temp = t0_dt + datetime.timedelta(seconds=tframe*i)
                # Central time of exposure
                t_mid_dt_temp = t0_dt_temp + datetime.timedelta(seconds=tframe*0.5)
                t0_temp = datetime.datetime.strftime(t0_dt_temp, "%Y-%m-%dT%H:%M:%S.%f")
                t_mid_temp = datetime.datetime.strftime(t_mid_dt_temp, "%Y-%m-%dT%H:%M:%S.%f")
                hdr["UTC0"] = (t0_temp, "exposure starting date and time")
                hdr["UTC"] = (t_mid_temp, "central exposure date and time")

                # Remove unsensitive pixels
                temp = data_temp[i, (ymin-1):ymax, (xmin-1):xmax]
                if args.mask:
                    # Mask not well corrected pixels
                    temp = np.where(mask==1, 1.0, temp)
                hdu[0].data = temp
                out = f"{basename}ms{i+1:04d}.fits"
                hdu.writeto(os.path.join(outdir, out), overwrite=True)


    # Only masking edge for 2-d fits
    if len(data_temp.shape)==2:
        # Use sensitive pixels
        temp = data_temp[(ymin-1):ymax, (xmin-1):xmax]
        if args.mask:
            # Mask not well corrected pixels
            temp = np.where(mask==1, 1.0, temp)
        hdu[0].data = temp
        hdu[0].header.add_history(
            f"[video2image] original fits : {fitsname}")
        out = f"{filename}_c.fits"
        hdu.writeto(os.path.join(out), overwrite=True)


if __name__ == "__main__":
    parser = ap(
        description="Convert video (3-d fits) to image (2-d fits).")
    parser.add_argument(
        "fits", type=str, 
        help="a reduced 3-d fits")
    parser.add_argument(
        "--mask", action="store_true", default=False, 
        help="mask not well corrected pixes")
    parser.add_argument(
        "--stack", type=int, default=None, 
        help="stacking number")
    parser.add_argument(
        "--stype", type=str, default="median", 
        help="stacking type")
    args = parser.parse_args()

    main(args)
