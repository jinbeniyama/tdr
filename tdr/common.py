#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Common functions for triccs data reduction.
"""
from astropy.io import fits as fits
from astropy.wcs import WCS as wcs
from astropy.time import Time
from scipy import interpolate
import datetime
import numpy as np
from astroquery.jplhorizons import Horizons


def utc2jd(utc):
  """
  utc : %Y-%m-%dT%H:%M:%S.%f
  """
  t = Time(str(utc), format='isot', scale='utc')
  return t.jd


def calc_JPLephem(
    asteroid, date0, date1, step, obscode=None, loc=None, air=False):
    """
    Calculate asteroid ephemeris.
  
    Parameters
    ----------
    asteroid : str
      asteroid name like "Ceres", "2019 FA" (should have space)
    date0 : str
      ephemeris start date like "2020-12-12"
    date1 : str
      ephemeris end date like "2020-12-12"
    step : str
      ephemeris step date like '1d' for 1-day, '30m' for 30-minutes
    obscode : str, optional
      IAU observation code. 
      371:Okayama Astronomical Observatory
      381:Kiso observatory
    loc : dict, optional
      dictionary like {lon=134.3356, lat=35.0253, elevation=449}
    air : bool, optional
      whether consider air pressure
  
    Return
    ------
    ephem : astropy.table.table.Table
      calculated ephemeris
    """
    if obscode:
        obj = Horizons(id=asteroid, location=obscode,
            epochs={'start':date0, 'stop':date1, 'step':step})
    elif loc:
        obj = Horizons(id=asteroid, location=loc,
            epochs={'start':date0, 'stop':date1, 'step':step})
    eph = obj.ephemerides(refraction=air)
    return eph


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
    loc : str
        MPC observatory code
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
        loc = None

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
        loc = None

    # TODO: Update
    if inst == "akeno":
        hdr_kwd = dict(
            datetime=None, date="DATE-OBS", time="UT-CEN",
            exp="EXPTIME", gain=None)
        loc = None

    if inst == "triccs":
        hdr_kwd = dict(
            datetime="UTC0", date="DATE-OBS", time="UT-CEN",
            exp="TFRAME", gain=None)
        # Note: Original keyword for datetime is 'GEXP-STR'.
        #       We use updated keywords 'UTC0' for (starting) datetime here.
        #       FYI, 'UTC' is central time.
        # 0.35 arcsec/pixel
        p_scale = 0.35
        loc = 371

    elif inst == "DECam":
        hdr_kwd = dict(
            datetime="DATE-OBS", date="DATE-OBS", time="DATE-OBS",
            exp="EXPTIME", gain=None)
        # 0.27 arcsec/pixel from fitsheader 
        p_scale = 0.27
        # Cerro Tololo-Dark Energy Camera (DECam) of The Dark Energy Survey, also see #807
        # from wikipedia
        loc = "W84"

    return hdr_kwd, p_scale, loc


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


def shift_sidereal(flist, stackmode="median"):
    """
    Shift 2d-fits files with velocity in pixel/s.

    Parameters
    ----------
    flist : array-like
        list of fits files
    stackmode : str, optional
        stacking mode (median or mean)

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
        img = hdu[0].data
        ny0, nx0 = img.shape
        if (ny0 == ny) & (nx0 == nx):
            pass
        else:
            # Skip unusual shape (for DECam)
            continue

        tmp = np.roll(img, dy, axis=0)
        print(tmp.shape)
        tmp = np.roll(tmp, dx, axis=1)
        # Save
        cube.append(tmp)

        if idx == 0:
            # Save shape since some DECam image has unusual shape......
            ny, nx = tmp.shape

    
    if stackmode == "median":
        img_shift = np.median(cube, axis=0)
    elif stackmode == "mean":
        img_shift = np.mean(cube, axis=0)
    else:
        assert False, "Invalid stackmode."
    return img_shift


def shift_nonsidereal(flist, hdr_kwd, target, loc, stackmode="median"):
    """
    Shift 2d-fits files a target name.
    Move target to a center.

    Parameters
    ----------
    flist : array-like
        list of fits files
    hdr_kwd : dictionary
        header keywords
    target : str
        target name
    loc : str
        Minor Planet Center observatory code
    stackmode : str, optional
        stacking mode (median or mean)

    Return
    ------
    img_shift : array-like
        shifted 2-d image
    """
    
    # Obtain the ephemeris of the target during the observations ==============
    hdu0 = fits.open(flist[0])
    hdu1 = fits.open(flist[-1])
    hdr0 = hdu0[0].header
    hdr1 = hdu1[0].header
    t_exp = hdr1[hdr_kwd["exp"]]
    # Starting time
    t0_str = hdr0[hdr_kwd["datetime"]]
    t1_str = hdr1[hdr_kwd["datetime"]]
    try:
        t1_str = datetime.datetime.strptime(t1_str, "%Y-%m-%dT%H:%M:%S.%f")
    except:
        t1_str = datetime.datetime.strptime(t1_str, "%Y-%m-%dT%H:%M:%S")
    # Ending time
    # +60 is necesary for ephemeris with 1min step to avoid 
    #   'ValueError: A value in x_new is above the interpolation range.'
    t1_end = t1_str + datetime.timedelta(minutes=(+30))
    t1_end = datetime.datetime.strftime(t1_end, "%Y-%m-%dT%H:%M:%S.%f")
    print(t1_str)
    print(t1_end)
    eph = calc_JPLephem(target, t0_str, t1_end, "1m", obscode=loc)
    jd_list, ra_list, dec_list = eph["datetime_jd"].tolist(), eph["RA"].tolist(), eph["DEC"].tolist()
    print(f"{np.min(jd_list)}")
    print(f"{np.max(jd_list)}")
    # To calculate ra and dec from jd
    f_ra = interpolate.interp1d(jd_list, ra_list, kind='linear')
    f_dec = interpolate.interp1d(jd_list, dec_list, kind='linear')
    # Obtain the ephemeris of the target during the observations ==============



    cube = []
    for idx, fi in enumerate(flist):
        # Obtain pixel coordinates of the standard ra and dec (cra0 and cdec0)
        hdu = fits.open(fi)
        hdr = hdu[0].header
        w = wcs(header=hdr)
        # Starting
        t_str = hdr[hdr_kwd["datetime"]]
        try:
            t_str = datetime.datetime.strptime(t_str, "%Y-%m-%dT%H:%M:%S.%f")
        except:
            t_str = datetime.datetime.strptime(t_str, "%Y-%m-%dT%H:%M:%S")
        t_exp = hdr[hdr_kwd["exp"]]
        # Central time
        t_cen = t_str + datetime.timedelta(seconds=t_exp/2.0)
        t_cen = datetime.datetime.strftime(t_cen, "%Y-%m-%dT%H:%M:%S.%f")
        print(t_cen)
        t_jd = utc2jd(t_cen)
        print(t_jd)
        
        # Obtain coordinates of a minor planet at the central time
        ra, dec = f_ra(t_jd), f_dec(t_jd)

        # Standard pixel coordinates
        if idx == 0:
            x0, y0 = w.all_world2pix(ra, dec, 0)

        # Convert equatorial coordinates to pixel coordinates
        x, y  = w.all_world2pix(ra, dec, 0)

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

    if stackmode == "median":
        img_shift = np.median(cube, axis=0)
    elif stackmode == "mean":
        img_shift = np.mean(cube, axis=0)
    else:
        assert False, "Invalid stackmode."
    return img_shift


