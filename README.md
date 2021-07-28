# TriCCS Data Reduction (TDR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](<mailto:beniyama@ioa.s.u-tokyo.ac.jp>)


## Overview
Data reduction of data taken with Seimei/TriCCS could be done in this repository.
Though being optimized for Seimei/TriCCS data, 
you can apply it for imaging data taken with other high-speed camera
(Tomo-e Gozen etc.).

### Standard Procedure
1. dark subtraction

2. flat fielding

3. split fits to 2d (for photometry using `Moving Object Photometry (movphot)`)


## Usage
Here 3-bands observation by Seimei/TriCCS is assumed.
All bands data could be analyzed by the same way.


### 1. dark subtraction

```
[template]
# Create master dark

# Do dark subtraction

```

```
[example]
# Create master dark

# Do dark subtraction
```


### 2. flat fielding
```
[template]
# Create master flat

# Do flat fielding

```

```
[example]
# Create master flat

# Do flat fielding
```


### 3. split fits to 2d
If you are going to use `Moving Object Photometry (movphot)` for photometry,
3-d fits cube should be splitted to multiple 2-d fits.

```
[template]
# Create master flat

# Do flat fielding

```

```
[example]
# Create master flat

# Do flat fielding
```





## Dependencies
This library is depending on `NumPy`.
Scripts are developed on `Python 3.7.10` and `NumPy 1.19.2`.


## LICENCE
This software is released under the MIT License.
