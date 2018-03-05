# FFIToolsPub

This is the public version of FFITools. FFITools is a software tool package for reduce wide field images and transit search. 

## Dependence: 

external softwares:

* fitsh 0.9.2 

* vartools 1.33

python modules:

* astrometry.net

* astropy>1.3

* numpy==1.11.1 

* scipy==0.18.1

* multiprocessing==0.70a1 

* matplotlib==1.5.3

* astroquery==0.3.2

* astropy>1.3

## Usage: 

Display a fits file: 

```
patools-display tess2019129080826-4-4-0016-s_ffic.fits
```

Display a fits file and project a star list based on its x, y coordinates (tess2019129080826-4-4-0016-s_ffic.fistar is a file result from source extraction and has column 2/3 for x/y coordinates ): 

```
patools-display tess2019129080826-4-4-0016-s_ffic.fits -x 2,3 -i tess2019129080826-4-4-0016-s_ffic.fistar
```

The following steps need to be operated in sequences: 

source extraction (in etc directory): 

```
patools-sampler fistar -c ETE6-cam4-ccd4.cfg --debug --logfile -
```

astrometry 

```
patools-sampler astrometry -c ETE6-cam4-ccd4.cfg --debug --logfile - 
```

photometry

```
patools-sampler phot -c ETE6-cam4-ccd4.cfg --debug --logfile -
```

# Authors

Chelsea Huang

Andras Pal
...
