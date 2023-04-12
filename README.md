# KT19 Processing and Visualization Software 
This repository contains software (src/ and tools/) for processing KT19 data colled by the Modular Aerial Sensing System (MASS) along with processing documentation. Data is included to help verify software is working properly. 

## Authors 
* [Nicholas Statom](https://airsea.ucsd.edu/people/) <<nstatom@ucsd.edu>>
* [Luke Colosi](https://lcolosi.github.io/)<<lcolosi@ucsd.edu>>

## How to use this repository

 using the MatLab scripts from this repository and processed [data](insert collection url) published to the UCSD library digital collections. To do so, follow these steps:

1. Make a local copy of this repository by either cloning or downloading it.

2. Download the processed [data](insert collection url), untar the files, and move all directories to `data` in the project root. After doing so, your directory tree should look like this:

```
WaveSpectrum/
├── data
│   ├── DelMAR2020
│   ├── SMODE2021
│   └── BATHY
├── figs
├── src
└── tools
```

3. Make sure that you have downloaded the R2022a version of MATLAB. Other versions may work, but small difference in keyword argument and outputs of functions may arise.   

4. If you follow the steps above you should be able to reproduce all figures, by running `figXX.m` from the `src` directory without having to adjust any paths.

**Note on MatLab software version:** Code was developed using MatLab R2022a. Version released before and after this may run into some errors. 

**Acknowledgements** : Code for processing KT19 data was originally written by Nick Statom. LC made adaptations to this code.
