LIC visualization of polarization vectors on ipole images.

Copyright 2020 Monika Moscibrodzka 

ipole-lic  version 1.0  (released July, 2020)

The files in this directory are part of ipole. 

ipole is a program that solves the equations of polarized radiative transfer
in covariant form, using ray tracing. It produces images in Stokes I,Q,U,V
at the location of a "camera."  
ipole-LIC makes visualization of polarization vectors on ipole images.
 

You are morally obligated to cite the following papera in any
scientific literature that results from use of any part of ipole:

[1] Moscibrodzka, Monika, and Gammie, Charles F.,
    Monthly Notices of the Royal Astronomical Society, Volume 475, Issue 1, p.43-54
    
ipole-lic is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


How to use the code?
make
./lic ipole.dat
three new images in ppm forma will be created
ipole_fnu.ppm - image of Stokes I in linear (or other) scale
ipole_lfnu.ppm - image of Stokes I in logarithmic scale
ipole_lic.ppm - lic image
