# TriCCS Data Reduction (TDR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](<mailto:beniyama@ioa.s.u-tokyo.ac.jp>)


## Overview
Data reduction of data taken with Seimei/TriCCS could be done in this repository.  
Though optimized for Seimei/TriCCS, 
you can apply it for imaging data taken with other high-speed cameras
(Tomo-e Gozen etc.).


Stackmode in shift2d, mean or median, may affect the photometric results. Please carefully check it.
J.B. used median mode for 2015 RN35 (Beniyama et al., 2023c) 
and mean mode for 2001 QJ142 (Beniyama et al. 2024).

### Procedures
1. Calibration (dark subtraction, flat-field correction)

2. Stacking fits by mean, median, etc.

3. Masking and splitting
(only for photometry using `Moving Object Photometry (movphot)` ([bitbucket](https://bitbucket.org/jin_beniyama/movphot/src/master/)) )

4. Common ID search

X. Spectroscopy with Pypelt


## Installing
```
# Install from PyPI
pip install tdr
```
or
```
# Install with `setup.py`
git clone https://jin_beniyama@bitbucket.org/jin_beniyama/triccs-data-reduction.git
python setup.py install
```


## Usage
Here g-band data taken with Seimei/TriCCS is considered.
All band data could be analyzed in the same way.

Fits data taken with TriCCS have a format like `TRCS00005180.fits`.

First 4 characters `TRCS` means the instrument *TriCCS*,
the next 7 characters are the exposure ID 

and the last 1 character is band identical number 
(`0` for g-band, `1` for r-band and `2` for i/z-band).

After each reduction stage, a prefix is added to the filename.
History can be checked in fits header as well.


### 1. Calibration
Here, consider the situation:
exposure time for an object frame is 10 s,
for a flat frame is 1 s,
and for dark frames are 1s and 10s (all in g-band) 
like below.

1. object `TRCS00000010.fits` (10 s)
2. dark for flat `TRCS00000020.fits` (1s)
3. dark for object `TRCS00000030.fits` (10s)
4. flat `TRCS00000040.fits` (1s, twilight field of view 1) 
and `TRCS00000050.fits` (1s, twilight field of view 2), 

Dark subtraction and flat-field correction are done as follows.


First, create master dark frame, which has prefix `d` like `dTRCS00000020.fits`.

The maximum count frame is not used for stacking,
which leads to avoiding cosmic rays or fast-moving object contamination.
```
[usage]
# Create master dark
makedark (3-d dark)

[example]
# Create master dark for flat
makedark.py TRCS00000020.fits
# Create master dark for object
makedark.py TRCS00000030.fits
```

Next, create a master normalized flat frame using master dark for flat frames
, which has prefix `f` like `fTRCS00000040.fits`.

```
[usage]
# Create master flat
makeflat.py --flat (3-d flat frames) --dark (2-d master dark)

[example]
# Create master flat
makeflat.py --flat TRCS00000040.fits TRCS00000050.fits --dark dTRCS00000020.fits
```

Reduce an object frame using both master dark and flat frames. 
The reduced object frame has prefix `r` like `rTRCS00000010.fits`.

```
[usage]
# Do dark subtraction and flat-field correction
reduce.py --obj (3-d object) --dark (2-d master dark) --flat (2-d master flat)

[example]
# Do dark subtraction and flat-field correction
reduce.py --obj TRCS00000010.fits --dark dTRCS00000020.fits --flat fTRCS00000040.fits
```

Finally, make 2-d fits from 3-d fits.
The resulting object frame is like `rTRCS00000010ms0005.fits`, in which 
`ms0005` means that this is the 5th frame if the original reduced fits.

```
[usage]
video2image.py (3-d reduced object)

[example]
video2image.py rTRCS00000010.fits
```


### 2. Stacking
#### 2.1. Pixed based stacking
Do stacking (i.e., make 2-d fits from 3-d fits) using pixel cooridnates.
Output fits have format like 
`maxrTRCS00000010.fits` (max),
`minrTRCS00000010.fits` (min),
`meanrTRCS00000010.fits` (mean) and
`medianrTRCS00000010.fits` (median).

```
[usage]
# Maximum stacking
stackfits.py (3-d reduced fits) max
# Minimum stacking
stackfits.py (3-d reduced fits) min
# Mean stacking
stackfits.py (3-d reduced fits) mean
# Median stacking
stackfits.py (3-d reduced fits) median

[example]
# Maximum stacking
stackfits.py rTRCS00000010.fits  max
# Minimum stacking
stackfits.py rTRCS00000010.fits  min
# Mean stacking
stackfits.py rTRCS00000010.fits  mean
# Median stacking
stackfits.py rTRCS00000010.fits  median
```

#### 2.2. WSC based stacking
Do stacking (i.e., make 2-d fits from 3-d fits) using World Coordinate System (WCS).
I am sorry that you have to prepare corrected fits by yourselves using such as astrometry.net.

```
[usage]
# Maximum stacking
stackfits.py (3-d reduced fits) max
# Minimum stacking
stackfits.py (3-d reduced fits) min

[example]
# Maximum stacking
stackfits.py rTRCS00000010.fits  max
# Minimum stacking
stackfits.py rTRCS00000010.fits  median
```


### 3. Masking and splitting
If you are going to use `Moving Object Photometry (movphot)`[(bitbucket)](https://bitbucket.org/jin_beniyama/movphot/src/master/) for photometry,
3-d fits cube should be split into multiple 2-d fits.
Cutting non-sensitive pixels and masking nor well corrected pixels are done
as well.

![pixel map](/TriCCS_pixel_map_20210817.jpg){width=30%}

```
[usage]
# Split fits into multiple 2-d fits while cutting non-sensitive pixels
video2image.py (reduced 3-d fits)
# + Masking not well corrected pixels (count is set to 1 for masked pixels)
video2image.py (reduced 3-d fits) --mask

[example]
# Split fits into multiple 2-d fits while cutting non-sensitive pixels
video2image.py rTRCS00000010.fits 
# + Masking not well corrected pixels (count is set to 1 for masked pixels)
video2image.py rTRCS00000010.fits  --mask
```
The masked and splitted frames hav suffix `ms` like `rTRCS00000010ms0001.fits`.
Output fits examples are as follows (when the number of frames is 3).

- `rTRCS00000010ms0001.fits`
- `rTRCS00000010ms0002.fits`
- `rTRCS00000010ms0003.fits`


### 4. Common ID search
If the wcs pasting failed for some fits,
it is necessary to extract common ID fits.


```
[example]
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
