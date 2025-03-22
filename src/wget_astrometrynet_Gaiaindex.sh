#! /bin/bash
# Obtain index file for astrometry
# See http://astrometry.net/doc/readme.html#getting-index-files and http://data.astrometry.net/ for details.
# This index files work well for Seimei/TriCCS data (FoV ~ 12.6 arcmic x 7.5 arcmin)
# index-5205-*
for ((i=0; i<48; i++)); do
    I=$(printf %02i $i)
    wget https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5205-$I.fits
done
