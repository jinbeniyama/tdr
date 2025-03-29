#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paste WCS (World Coordinate System) with `astrometry.net`.

Example
-------
# Paste wcs to all fits listed in fitslist.txt at once with
# p_scale 0.34 (arcsec/pixel) and binning n=3.
# Pasted fits files are save in `wcs` directory.
>>> paste_wcs 0.34 --flist fitslist.txt -n 3 --outdir wcs

"""
import os
import subprocess
import numpy as np
from argparse import ArgumentParser as ap
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy import units as u

# Change for your emvironment
__solve_field = "/opt/homebrew/bin/solve-field"


def get_astrometry_scale(header, p_scale):
    """ 
    Calculate index scales for `astrometry.net`.

    Parameter
    ---------
    header : Header 
      fits header object

    Returns
    -------
      l_scale : float 
        lower bound of image scale in arcmin
      u_scale : float
        upper bound of image scale in arcmin
    """

    nx1 = header.get("NAXIS1")
    nx2 = header.get("NAXIS2")
    s, d = np.min([nx1,nx2]), nx1+nx2
    # pixel scale in units of arcmin/pix!!
    p_scale = p_scale/60.0 
    l_scale, u_scale = s*p_scale, d*p_scale
    return l_scale, u_scale


def solve_field(
    infits, fitsdir, outdir, nsample, depth, radius, odds, sigma, p_scale):
    """ 
    Calibrate astrometry infomation of the 2-d fits.

    Parameters
    ----------
    infits : str
      Input fits file to be solved
    fitsdir : str
      fits directory of wcs to be pasted
    outdir : str
      output directory of wcs pasted fits
    nsample : int
      number of down sampling
    depth : int
      star brightness
    radius : float
      search radius
    oods : float
      need odds to be solved
    sigma : float
      detection threshold (not used)
    p_scale : float
      pixel scale in arcsec
    """

    # Open the fits
    infits = os.path.join(fitsdir, infits)
    src = fits.open(infits)[0]
    hdr = src.header
    ra  = hdr.get("RA")
    dec = hdr.get("DEC")
    radec = f"{ra} {dec}"
    c = SkyCoord(radec, unit=(u.hourangle, u.degree))
    # Convert to degree
    ra, dec = c.ra.degree, c.dec.degree
 
    # Directory for intermediate product
    out_med = "wcsout"

    # Obtain file identification part
    frame_id = infits.split("/")[-1].split(".")[0]

    # wcs pasted fits file
    outfits = f"{frame_id}w.fits"
    outfits = os.path.join(outdir, outfits)

    # Calculate scales with the pixel scale
    scale_low, scale_high = get_astrometry_scale(hdr, p_scale)

    cmd = [__solve_field, infits, "--no-plots", "--no-tweak", 
          "--depth", str(depth), "--ra", str(ra), "--dec", str(dec), 
          "--radius", str(radius), "--scale-low", str(scale_low), 
          "--scale-high", str(scale_high), "--scale-units", "arcminwidth", 
          "--sigma", str(sigma), "--new-fits", outfits, 
          "--downsample", str(nsample), "--overwrite", "--dir", str(out_med), 
          "--odds-to-solve", str(odds), "--no-verify",
          ]
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        print("CalledProcessError due to the lack of stars?")
    except ZeroDivisionError:
        print("ZeroDivisionError due to the lack of stars?")


if __name__ == "__main__":
    parser = ap(description="Paste wcs to fits obtained with Seimei/TriCCS")
    parser.add_argument(
        "--fits", default=None, 
        help="fits to be solved")
    parser.add_argument(
        "--flist", default=None,
        help="fitslist text to be solved")
    parser.add_argument(
        "--p_scale", type=float, default=0.35,
        help="pixel scale in arcsec/pix (for Seimei/TriCCS)")
    parser.add_argument(
        "--outdir", type=str, default=".",
        help="output directory")
    parser.add_argument(
        "--fitsdir", type=str, default=".",
        help="fits directory")
    parser.add_argument(
        "-n", "--nsample", type=int, default=3)
    parser.add_argument(
        "-d", "--depth", type=int, default=100,
        help="search depth")
    parser.add_argument(
        "-o", "--output", action="store", default=None)
    parser.add_argument(
        "-r", "--radius", type=float, default=0.5,
        help="search radius in degree")
    parser.add_argument(
        "--odds", type=float, default=1e11,
        help="need odds")
    parser.add_argument(
        "--sigma", type=float, default=3,
        help="ditection threshold")
    args = parser.parse_args()

    assert (args.fits or args.flist), "Invalid Input data"
    
    # Create the directory
    os.makedirs(args.outdir, exist_ok=True)

    # For 2-d fits
    # At once
    if args.flist:
        with open(args.flist, "r") as f:
            for line in f:
                tfits = line.strip("\n")
                solve_field(
                    tfits, args.fitsdir, args.outdir, args.nsample, args.depth, 
                    args.radius, args.odds, args.sigma, args.p_scale)
    # 1 fits file
    elif args.fits:
        solve_field(
            args.fits, args.fitsdir, args.outdir, args.nsample, args.depth, 
            args.radius, args.odds, args.sigma, args.p_scale)
