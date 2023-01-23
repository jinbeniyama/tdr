#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Common functions for triccs data reduction.
"""
from astropy.io import fits as fits
from astropy.wcs import WCS as wcs
from astropy.time import Time
import datetime
import numpy as np


def instinfo(inst):
    """Get header keywords and pixel scale.

    Parameter
    ---------
    inst : str
        instrument (e.g., murikabushi, saitama)

    Returns
    -------
    hdr_kwd : dictionary
        dictionary of header keywords
    p_scale : float
        pixel scale
    """

    if inst == "murikabushi":
        # Example
        # DATE-OBS= '2020-10-28'         /  [yyyy-mm-dd] Observation start date
        # UT      = '14:10:15.84'        / [HH:MM:SS.SS] Universal Time at start
        # EXPTIME =               30.000 / [sec] Exposure time
        # GAIN    =                 1.70 / [e-/ADU] CCD gain
        hdr_kwd = dict(
            datetime=None, date="DATE-OBS", time="UT",
            exp="EXPTIME", gain="GAIN")
        # 0.72 arcsec/pixel
        p_scale = 0.72

    # TODO: Update
    if inst == "saitama":
        # Example
        # DATE-OBS= '2021-02-11T11:17:41' / Date & Time of start of obs.in UTC
        # EXPTIME =               10.000 / Exposure Time [sec]
        hdr_kwd = dict(
            datetime="DATE-OBS", date=None, time=None,
            exp="EXPTIME", gain=None)
        # 0.73 arcsec/pixel
        p_scale = 0.73

    # TODO: Update
    if inst == "akeno":
        hdr_kwd = dict(
            datetime=None, date="DATE-OBS", time="UT-CEN",
            exp="EXPTIME", gain=None)

    if inst == "triccs":
        hdr_kwd = dict(
            datetime="UTC0", date="DATE-OBS", time="UT-CEN",
            exp="TFRAME", gain=None)
        # Note: Original keyword for datetime is 'GEXP-STR'.
        #       We use updated keywords 'UTC0' for (starting) datetime here.
        #       FYI, 'UTC' is central time.
        # 0.35 arcsec/pixel
        p_scale = 0.35

    return hdr_kwd, p_scale


def obtain_fitstime(flist, hdr_kwd):
    """
    Extract central time of fits.
    Predict orbit of moving object by fitting.
    The function is for 2-d fits and optimized for murikabushi and TriCCS.

    Parameters
    ----------
    hdr_kwd : dict
        header keyword
    flist : array-like
        fits file list

    Return
    ------
    t_list : array of float
        list of central times in unix time (second)
    """

    t_list = []
    for idx, fi in enumerate(flist):
        src = fits.open(fi)[0]
        hdr = src.header
          
        if hdr_kwd["datetime"]:
            exp_start = hdr[hdr_kwd["datetime"]]
        else:
            exp_start = f"{hdr[hdr_kwd['date']]}T{hdr[hdr_kwd['time']]}"

        exp_frame = hdr[hdr_kwd["exp"]]
        try:
            obs_start = datetime.datetime.strptime(
                exp_start, "%Y-%m-%dT%H:%M:%S.%f")
            obs_center = obs_start + datetime.timedelta(seconds=exp_frame/2.0)
            obs_center = datetime.datetime.strftime(
                obs_center, "%Y-%m-%dT%H:%M:%S.%f")
        except ValueError as e: 
            #print("  Value error happens.:")
            #print(e)
            obs_start = datetime.datetime.strptime(
                exp_start, "%Y-%m-%dT%H:%M:%S")
            obs_center = obs_start + datetime.timedelta(seconds=exp_frame/2.0)
            obs_center = datetime.datetime.strftime(
                obs_center, "%Y-%m-%dT%H:%M:%S")
        
        if idx==0:
            print(f"  Observation starting time: {obs_center}")
        t = Time(str(obs_center), format='isot', scale='utc')
        t_mjd = t.mjd
        # seconds since 1970.0 (UTC)
        t_s = t.unix
        t_list.append(t_s)
    return t_list


def shift_w_velocity(flist, hdr_kwd, v_ra, v_dec):
    """
    Shift 2d-fits files with velocity in pixel/s.

    Parameters
    ----------
    flist : array-like
        list of fits files
    hdr_kwd : dictionary
        header keywords
    v_ra, v_dec : float
        velocity in arvsec/pixel

    Return
    ------
    img_shift : array-like
        shifted 2-d image
    """
    
    # Obtain list of time in unix time
    t_unix = obtain_fitstime(flist, hdr_kwd)

    cube = []
    for idx, fi in enumerate(flist):
        hdu = fits.open(fi)
        hdr = hdu[0].header
        data = hdu[0].data

        # Calculate shift length in pixel using 
        # time and velocity
        t_elapse = t_unix[idx]
        ax1 = int(t_elapse*v_ra)
        ax2 = int(t_elapse*v_dec)
        # Shift
        print(f"  Shift x,y = {ax1} pix, {ax2} pix")
        tmp = np.roll(data, ax2, axis=0)
        tmp = np.roll(tmp, ax1, axis=1)
        # Save
        cube.append(tmp)

    img_shift = np.median(cube, axis=0)
    return img_shift


def shift_sidereal(flist):
    """
    Shift 2d-fits files with velocity in pixel/s.

    Parameters
    ----------
    flist : array-like
        list of fits files
    hdr_kwd : dictionary
        header keywords

    Return
    ------
    img_shift : array-like
        shifted 2-d image
    """

    # Obtain wcs info of the first fits
    fi0 = flist[0]
    hdu = fits.open(fi0)
    hdr = hdu[0].header
    ny, nx = hdu[0].data.shape
    w = wcs(header=hdr)
    # Central x and y of the first fits in pixel 
    x0, y0 = nx/2., ny/2.
    # Central ra and dec of the first fits
    cra0, cdec0  = w.all_pix2world(x0, y0, 0)

    cube = []
    for idx, fi in enumerate(flist):
        # Obtain pixel coordinates of the standard ra and dec (cra0 and cdec0)
        hdu = fits.open(fi)
        hdr = hdu[0].header
        w = wcs(header=hdr)
        x, y  = w.all_world2pix(cra0, cdec0, 0)
        # Shift length
        dx, dy = x0-x, y0-y
        dx, dy = int(dx), int(dy)
        print(
            f"  fits {idx+1:03d}: (cx, cy) = ({x:.1f}, {y:.1f}), "
            f"(dx, dy) = ({dx:02d}, {dy:02d}) pix")

        # Shift
        tmp = np.roll(hdu[0].data, dy, axis=0)
        tmp = np.roll(tmp, dx, axis=1)
        # Save
        cube.append(tmp)

    img_shift = np.median(cube, axis=0)
    return img_shift


def shift_nonsidereal(flist, hdr_kwd, target):
    """
    Shift 2d-fits files with velocity in pixel/s.

    Parameters
    ----------
    flist : array-like
        list of fits files
    hdr_kwd : dictionary
        header keywords
    target : str
        target name

    Return
    ------
    img_shift : array-like
        shifted 2-d image
    """
    return img_shift

