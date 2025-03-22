#! /bin/bash
# Obtain index file for astrometry
# See http://astrometry.net/doc/readme.html#getting-index-files and http://data.astrometry.net/ for details.
# index5203-xx files work well for Seimei/TriCCS data (FoV ~ 12.6 arcmic x 7.5 arcmin)

# index-5205-*
#for ((i=0; i<48; i++)); do
#    I=$(printf %02i $i)
#    wget https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5205-$I.fits
#done

# index-5204-*
#for ((i=0; i<48; i++)); do
#    I=$(printf %02i $i)
#    wget https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5204-$I.fits
#done

# index-5203-*
for ((i=0; i<48; i++)); do
    I=$(printf %02i $i)
    wget https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5203-$I.fits
done
