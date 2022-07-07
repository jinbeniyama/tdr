#!/bin/bash

#wcscheck_TriCCS_auto.sh cutedge wcs 5 18 5 18
# 2022-04-11 
# RID0 = 5, IRD1 = 18 is good
# rTRCS00135172ms0293.fits
# WID0 = 5, WRD1 = 18 is good
# rTRCS00135172ms0293w.fits

HOME=`pwd`
RDIR=${1}
WDIR=${2}

RID0=${3}
RID1=${4}
WID0=${5}
WID1=${6}

echo "Parameters RAW DIRECTORY: ${RDIR}"
echo "           WCS DIRECTORY: ${WDIR}"
echo "           (RID0, RID1) = (${RID0}, ${RID1})"
echo "           (WID0, WID1) = (${WID0}, ${WID1})"

# Save flist before wcs pasting assuming g,r,i/z bands
cd ${HOME}/${RDIR}
for i in `seq 0 2`
do 
    RLIST[i]=`ls | grep ${i}ms`
done


# Save flist after wcs pasting assuming g,r,i/z bands
cd ${HOME}/${WDIR}
for i in `seq 0 2`
do 
    WLIST[i]=`ls | grep ${i}ms`
done


# Create diff list 
# Saved as idlist_0.csv etc.
cd ${HOME}
for i in `seq 0 2`
do
    wcsIDcheck_TriCCS.py --rlist ${RLIST[${i}]} --wlist ${WLIST[${i}]} --rid0 ${RID0} --rid1 ${RID1} --wid0 ${WID0} --wid1 ${WID1} --rdir ${RDIR} --wdir ${WDIR} --band $i
done

# Show diff
wcsIDshow_TriCCS.py idlist_0.csv idlist_1.csv idlist_2.csv

# move common fits to wcs_common
wcsIDmove_TriCCS.sh idlist_0_common.csv idlist_1_common.csv idlist_2_common.csv
