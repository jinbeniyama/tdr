#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create video from fits data.
(gif movie is heavy.)

ffmpeg should be installed by e.g.,
> brew install ffmpeg
"""
import os
import numpy as np
from argparse import ArgumentParser as ap
import astropy.io.fits as fits
from astropy.wcs import WCS as wcs
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and
        new axes must divide old ones.

    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray

if __name__ == "__main__":
    parser = ap(description="Create video from 2d fits files.")
    parser.add_argument("fits", nargs="*", help="fits")
    parser.add_argument(
        "--xr", nargs=2, type=int, default=None, help="plot x range in pixel")
    parser.add_argument(
        "--yr", nargs=2, type=int, default=None, help="plot y range in pixel")
    parser.add_argument(
        "--fps", type=float, default=30,
        help="frames per second")
    parser.add_argument(
        "--dpi", type=int, default=20,
        help="dots per inch")
    parser.add_argument(
        "--rot", type=int, default=0,
        help="Rotate image by hand in degree")
    parser.add_argument(
        "--lr", action="store_true", default=False,
        help="Flip image by hand")
    parser.add_argument(
        "--ud", action="store_true", default=False,
        help="Flip image by hand")
    parser.add_argument(
        "--cmap", type=str, default="Reds",
        help="Color scheme")
    parser.add_argument(
        "--out", default="video_from_fits.mp4", help="output file")
    args = parser.parse_args()

    if args.xr:
        xmin, xmax = args.xr
    else:
        xmin, xmax = 0, 2000

    if args.yr:
        ymin, ymax = args.yr
    else:
        ymin, ymax = 0, 2000

    im_list = []

    # Read data to fix the image size
    hdu = fits.open(args.fits[0])
    data = hdu[0].data
    img = data[ymin:ymax, xmin:xmax]
    ny, nx = img.shape

    # TODO: update
    if args.rot:
        rat = nx/ny
    else:
        rat = ny/nx
    wx, wy = 10, 10*rat

    fig = plt.figure(figsize=(wx, wy))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)



    # For multiple 2-d fits
    for idx, infits in enumerate(args.fits):
      print("Input fits")
      print(infits)
      hdu = fits.open(infits)
      hdr = hdu[0].header
      data = hdu[0].data
      img = data[ymin:ymax, xmin:xmax]
      ny, nx = img.shape
      w = wcs(header=hdr)
      
      # TODO Update
      if args.lr:
          img = np.fliplr(img)
      if args.ud:
          img = np.flipud(img)
      # Rotate image by hand 
      if args.rot:
          # Rotate clockwise in python image, not ds9 image
          img = np.rot90(img, 1)
      else:
          # Assuming the vertical direction of the image is NS.
          # These are valid for usual data,
          # whereas not valid for TriCCS data on 2022-12 (due to bugs?)

          # ra, dec of (0, 0)
          ra0, dec0  = w.all_pix2world(0, 0, 0)
          # ra, dec of (nx, ny)
          ra1, dec1  = w.all_pix2world(nx, ny, 0)
          print(f"    Coordinates of {infits}")
          print(f"      (ra0, dec0) = ({ra0:.2f}, {dec0:.2f})")
          print(f"      (ra1, dec1) = ({ra1:.2f}, {dec1:.2f})")
          # The coor. of figure created with matplotlib are as follows.
          # (ny,0)   .... (nx, ny)
          # ...          .......
          # (0, 0) .... (nx, 0)
          if dec0 < dec1:
              northistop = True
          else:
              northistop = False

          if ra0 > ra1:
              eastisleft = True
          else:
              eastisleft = False

          # always North is to the top, East is to the left
          if northistop and eastisleft:
              pass
          elif northistop and (not eastisleft):
              img = np.fliplr(img)
          elif (not northistop) and eastisleft:
              img = np.flipud(img)
          elif (not northistop) and (not eastisleft):
              img = np.flipud(img)
              img = np.fliplr(img)

      img_median = np.median(img[img!=0])
      sigma = 10
      vmin, vmax = -5*sigma, 5*sigma
      # Subtract background
      img -= img_median
      im = ax.imshow(img, vmin=vmin, vmax=vmax, cmap=args.cmap)
      im_list.append([im])
          

    ani = animation.ArtistAnimation(fig, im_list, interval=1)
    
    out = args.out
    outtype = out.split(".")[-1]
    if outtype == "gif":
        writer = "imagemagick"
    elif outtype == "mp4":
        writer = "ffmpeg"
    ani.save(out, writer=writer, fps=args.fps, dpi=args.dpi)

