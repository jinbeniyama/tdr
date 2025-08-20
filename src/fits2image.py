#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert a 2D FITS image to a static image (jpg/png).
"""

import numpy as np
from argparse import ArgumentParser as ap
import astropy.io.fits as fits
from astropy.wcs import WCS as wcs
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = ap(description="Convert 2D FITS to image file.")
    parser.add_argument("fits", help="input 2D fits file")
    parser.add_argument(
        "--xr", nargs=2, type=int, default=None, help="plot x range in pixel")
    parser.add_argument(
        "--yr", nargs=2, type=int, default=None, help="plot y range in pixel")
    parser.add_argument(
        "--rot", type=int, default=0, help="Rotate image by hand in degree")
    parser.add_argument(
        "--lr", action="store_true", default=False, help="Flip image left-right")
    parser.add_argument(
        "--ud", action="store_true", default=False, help="Flip image up-down")
    parser.add_argument(
        "--cmap", type=str, default="gray_r", help="Color scheme")
    parser.add_argument(
        "--dpi", type=int, default=150, help="Dots per inch")
    parser.add_argument(
        "--out", default="fits_image.jpg", help="output image file (jpg/png)")
    args = parser.parse_args()

    # Read FITS
    hdu = fits.open(args.fits)[0]
    hdr = hdu.header
    data = hdu.data

    # Cut out range
    if args.xr:
        xmin, xmax = args.xr
    else:
        xmin, xmax = 0, data.shape[1]
    if args.yr:
        ymin, ymax = args.yr
    else:
        ymin, ymax = 0, data.shape[0]

    img = data[ymin:ymax, xmin:xmax]

    # Flips and rotations
    if args.lr:
        img = np.fliplr(img)
    if args.ud:
        img = np.flipud(img)
    if args.rot:
        img = np.rot90(img, args.rot // 90)

    # Background subtraction (optional)
    img_median = np.median(img[img != 0])
    img = img - img_median

    # Scale (simple: ±5σ around 0)
    sigma = np.std(img)
    vmin, vmax = -5 * sigma, 5 * sigma

    # Plot
    fig, ax = plt.subplots()
    ax.axis("off")
    im = ax.imshow(img, vmin=vmin, vmax=vmax, cmap=args.cmap, origin="lower")

    # Save
    plt.savefig(args.out, dpi=args.dpi, bbox_inches="tight", pad_inches=0)
    print(f"Saved {args.out}")
