# read in fits file and converts it to txt file that can be read in in lic.c 
import ehtim as eh

# reading in the averaged images, these are already blurred so no extra blurring is required
im=eh.image.load_fits('M87_lo_3601_polarimetric_average_image_256.fits',polrep='stokes')

print('res:',im.xdim)
print('FOV=',im.fovx()/eh.RADPERUAS,'uas')

# dump txt
im.save_txt("image.txt")
