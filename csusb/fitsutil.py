from __future__ import division


import csusb
import os, pyfits
import numpy as np


def fitsread(imgname, header = False):
    """
    Read CSUSB telescope FITS image cube.
    
    Parameters
    ----------
    image : string
        FITS image name
        
    header : boolean
        Return FITS image header?
        
    Returns
    -------
    img_data : numpy array
        2D or 3D numpy array
    """
    try:
        if header:
            img_data, header = pyfits.getdata(imgname, ignore_missing_end = True, header = True)
            return img_data, header
        else:
            img_data = pyfits.getdata(imgname, ignore_missing_end = True)
            return img_data
    except IOError:
        print "FITSREAD: Unable to open FITS image %s" %imgname
    
    return
    
    
def fitswrite(img, imgname, **kwargs):
    """
    Write FITS image to disk.
    
    Parameters
    ----------
    img : numpy array
        2D or 3D numpy array
        
    imgname : string
        Name of the output FITS image
        
    Optional Keywords
    -----------------
    header : pyFITS header
        FITS header object
        
    Return
    ------
    None
    """
    try:
        if kwargs.has_key('header'):
            hdu = pyfits.PrimaryHDU(img, header = kwargs['header'])
        else:
            hdu = pyfits.PrimaryHDU(img)
        hdu.writeto(imgname)
    except IOError:
        print "FITSWRITE: Unable to write FITS image %s. Stopping." %imgname
    
    return
    
    
def fitshead(imgname):
    """
    Read CSUSB telescope FITS image header.
    
    Parameters
    ----------
    image : string
        FITS image name
        
    Returns
    -------
    img_header : python dictionary
        Dictionary of image header keywords
    """
    try:
        img_header = pyfits.getheader(imgname, ignore_missing_end = True)
        return img_header
    except IOError:
        print "FITSHEAD: Unable to open FITS image %s. Stopping." %imgname
    
    return
    
    
def fitscombine(inputimgs, combine = "average", nframes = 100, outfile = ""):
    """
    Combine FITS image frames of CSUSB telescope instrument 3D image cubes.
    
    Parameters
    ----------
    imgname : string
        FITS image cube name
        
    combine : string
        Combine type - average, sum or median
    
    nframes : int
        Number of frames to combine
        
    outfile : string
        Name of the output FITS image
        
    Returns
    -------
    None
    """
    # Read the image names to combine
    imgnames = []
    if inputimgs[0] == "@":
        with open(inputimgs[1:], "r") as fd:
            for line in fd.readlines():
                if len(line) > 1:
                    imgnames.append(line.replace("\n", ""))
    else:
        imgnames = inputimgs.split(",")
    
    # Determine image dimension 
    img_header = fitshead(imgnames[0])
    nx, ny = img_header["NAXIS1"], img_header["NAXIS2"]
    nframes = len(imgnames)
    
    # Generate a numpy 3d array
    img_cube = np.zeros([nframes,ny,nx], dtype = np.float32)
    for i in range(len(imgnames)):
        img_data = fitsread(imgnames[i], header = False)
        img_cube[i,:,:] = img_data        
    
    comb_img = csusb.imcombine(img_cube)            
        
    if outfile == "":
        fitswrite(comb_img, os.path.splitext(imgnames[0])[0] + "_comb.fits", header = img_header)
    else:
        fitswrite(comb_img, outfile, header = img_header)
        
    return
