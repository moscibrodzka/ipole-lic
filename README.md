Line Integral Convolution by Cabral 1993.

written by Monika Moscibrodzka

w/ some functions are borrowed from ipole (Moscibrodzka and Gammie 2018, Noble et et. 2007)


Instruction to make a LIC map:

1. convert fits file with image into ASCII file .dat (i,j,I,Q,U,V)

   python script using ehtim provided, here one can also blur the image

2. in decs.h

   set resolution of the image (res: 256x256 by default)
   
   set FOV (fov: 128 muas by default)
   
   choose color scale of the output ppm files
   
        #define RAINBOW 1
        #define AFMHOT 0
        #define BW 0


3. compile program typing: make

4. run: ./lic file_w_image.txt

5. program will produce 2 images in ppm format

   image_fnu.ppm - total intensity map
      
   image_lic.ppm - LIC map
