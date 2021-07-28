# TriCCS Data Reduction (TDR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](<mailto:beniyama@ioa.s.u-tokyo.ac.jp>)


## Overview
Data reduction of data taken with Seimei/TriCCS could be done in this repository.  
Though being optimized for Seimei/TriCCS, 
you can apply it for imaging data taken with other high-speed camera
(Tomo-e Gozen etc.).

### Procedures
1. Calibration (dark subtraction, flat-field correction)


2. Stack fits
mean, median etc.


3. Split fits 
(only for photometry using `Moving Object Photometry (movphot)`)

4. Common ID search

## Installing (in preparation)
```
pip install tdr
```


## Usage
Here g-band data taken with Seimei/TriCCS is condidered.
All bands data could be analyzed by the same way.

Fits data taken with TriCCS have format like `TRCS00005180.fits`.

First 4 characters `TRCS` means the instrument *TriCCS*,
next 7 characters are the exposure ID and 
the last 1 character is band identical number (`0` for g-band, `1` for r-band and `2` for i/z-band).

After each reduction stage, a prefix is added to the filename.
History can be checked in fitsheader as well.


### 1. Calibration
Here, consider the situation:
exposure time for object frame is 10 s,
exposure time for flat frame is 1 s 
and 1s, 10s dark frames are obtained (all in g-band).

The obtained fits are 

1. object `TRCS00000010.fits` (10 s)
2. dark for flat `TRCS00000020.fits` (1s)
3. dark for object `TRCS00000030.fits` (10s)
4. flat `TRCS00000040.fits` (1s)
.

The, dark subtraction and flat fielding are done as following.


First, create master dark frame.

The master dark has prefix `d` like `dTRCS00000020.fits`.

```
[usage]
# Create master dark
makedark (3-d dark)

[example]
# Create master dark for flat
makedark TRCS00000020.fits
# Create master dark for object
makedark TRCS00000030.fits
```

Next, create master normalized flat frame using master dark for flat frame.

The master flat has prefix `f` like `fTRCS00000040.fits`.

```
[usage]
# Create master flat
makeflat --flat (3-d flat) --dark (2-d master dark)

[example]
# Create master flat
makeflat --flat TRCS00000040.fits --dark dTRCS00000020.fits
```

Finally, reduce object frame using both master dark and flat frames. 

The reduced object frame has prefix `r` like `rTRCS00000010.fits`.

```
[usage]
# Do dark subtraction and flat-field correction
reduce --obj (3-d object) --dark (2-d master dark) --flat (2-d master flat)

[example]
# Do dark subtraction and flat-field correction
reduce --obj TRCS00000010.fits --dark dTRCS00000020.fits --flat fTRCS00000040.fits
```

### 2. Stacking
Output fits has format like `meanrTRCS00000010.fits` (mean) and
`medianrTRCS00000010.fits` (median).

```
[usage]
# Mean stacking
stackfits (3-d reduced fits) --mean
# Median stacking
stackfits (3-d reduced fits) --median

[example]
# Mean stacking
stackfits rTRCS00000010.fits  --mean
# Median stacking
stackfits rTRCS00000010.fits  --median
```


### 3. Splitting
If you are going to use `Moving Object Photometry (movphot)` for photometry,
3-d fits cube should be splitted to multiple 2-d fits.

```
[usage]
# Cut and create multiple 2-d fits
cut3dfits (3-d fits)

[example]
# Cut and create multiple 2-d fits
cut3dfits rTRCS00000010.fits
```
Output fits are as follows (when number of frame is 3).
```
rTRCS00000010c0001.fits
rTRCS00000010c0002.fits
rTRCS00000010c0003.fits
``` 

### 4. Common ID search
If wcs pasting failed for some fits,
it is necessary to extract common ID fits.

```
[usage]
# Extract common fits ID from g,r,i bands list (glist.txt, rlist.txt, ilist.txt)
cat glist.txt | awk '{print substr($0,17,3)}' > gID.txt
cat rlist.txt | awk '{print substr($0,17,3)}' > rID.txt
cat ilist.txt | awk '{print substr($0,17,3)}' > iID.txt
# Create common ID list 
commonIDsearch.py gID.txt rID.txt iID.txt --pre rTRCS00001260ce0 
--post w.fits > glist_common.txt
```

## Dependencies
This library is depending on `NumPy`.
Scripts are developed on `Python 3.7.10` and `NumPy 1.19.2`.


## LICENCE
This software is released under the MIT License.
