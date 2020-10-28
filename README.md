# Europe-wide phenology from Sentinel-2: Calibration with eddy covariance, PhenoCam and ground phenology observations

## Authors
Feng Tian, Zhanzhang Cai, Hongxiao Jin, Koen Hufkens, Helfried Scheifinger, Torbern Tagesson, Xiaoye Tong, Jonas Ardö, & Lars Eklundh

## Description
This repository contains all the scripts used for the calibration of Sentinel-2 vegetation phenology across Europe.

Sentinel-2 reflectance data were downloaded from Google Earth Engine: https://code.earthengine.google.com/c4d374f32991890c0b9d3395fe87429d. VIs are then calculated in R, followed by smoothing and phenometrics (SOS and EOS) extraction with the TIMESAT software. The final smoothed curves and phenometrics are stored as csv files. 

Contacts: Feng Tian (tian.feng[at]whu.edu.cn) and Lars Eklundh (lars.eklundh[at]nateko.lu.se)
