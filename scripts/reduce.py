#!/usr/bin/env python

"""
    --------------------------------------------------------------------------
    Routine to reduce raw science images from CSUSB telescope.
    
    Usage: python reduce.py [options] image
        
                                
    Author:
        Navtej Saini

    Organization:
        Caltech, Pasadena, CA, USA

    Version:
        21 January 2016     0.1     Initial implementation 
    --------------------------------------------------------------------------        
"""

import os, sys
import numpy as np
from StringIO import StringIO
from optparse import OptionParser

import csusb


def process(sci_files, bias_files, dark_files, flat_files, threshold, output):
    """
    Entry point function to process science images.
    
    Parameters
    ----------
    sci_files : string
        Science image file names
        
    bias_files : string
        Either the master bias image or comma separated bias images
        
    dark_files : string
        Either master dark image or list of comma separated dark frames
        
    flat_files : string
        Either master flat image or list of comma separated flat frames
        
    threshold : float
        Threshold for normalized fat field (value between 0 and 1.0). 
        Default is 0.8.
        
    Returns
    -------
    None 
    """
    print "REDUCE: CSUSB Image Reduction Routine"
    
    threshold = float(threshold)
         
    # All science images
    if sci_files[0] == "@":
        infile = sci_files[1:]
        
        if not os.path.exists(infile):
            print "  Not able to locate file %s" %infile
        
        sci_img_fnames = []
        with open(infile, "r") as fd:
            for line in fd.readlines():
                if len(line) > 1:
                    sci_img_fnames.append(line.replace("\n", ""))
    else:
        sci_img_fnames = sci_files.split(",")

    # All bias images
    if bias_files[0] == "@":
        infile = bias_files[1:]
        
        if not os.path.exists(infile):
            print "  Not able to locate file %s" %infile
        
        bias_img_fnames = []
        with open(infile, "r") as fd:
            for line in fd.readlines():
                if len(line) > 1:
                    bias_img_fnames.append(line.replace("\n", ""))
    else:
        bias_img_fnames = bias_files.split(",")


    # All dark images
    if dark_files[0] == "@":
        infile = dark_files[1:]
        
        if not os.path.exists(infile):
            print "  Not able to locate file %s" %infile
        
        dark_img_fnames = []
        with open(infile, "r") as fd:
            for line in fd.readlines():
                if len(line) > 1:
                    dark_img_fnames.append(line.replace("\n", ""))
    else:
        dark_img_fnames = dark_files.split(",")


    # All flat images
    if flat_files[0] == "@":
        infile = flat_files[1:]
        
        if not os.path.exists(infile):
            print "  Not able to locate file %s" %infile
        
        flat_img_fnames = []
        with open(infile, "r") as fd:
            for line in fd.readlines():
                if len(line) > 1:
                    flat_img_fnames.append(line.replace("\n", ""))
    else:
        flat_img_fnames = flat_files.split(",")


    # Read bias images and generate an image cube
    print "  Reading Bias images"
    img_header = csusb.fitshead(bias_img_fnames[0])
    nx, ny = img_header['NAXIS1'], img_header['NAXIS2']
    nimgs = len(bias_img_fnames)
    bias_img_cube = np.zeros([nimgs,ny,nx], dtype = np.float32)
    for i in range(nimgs):
        bias_img_cube[i,:,:] = csusb.fitsread(bias_img_fnames[i])

    # Read dark images and generate an image cube
    print "  Reading Darks"
    img_header = csusb.fitshead(dark_img_fnames[0])
    nx, ny = img_header['NAXIS1'], img_header['NAXIS2']
    nimgs = len(dark_img_fnames)
    dark_img_cube = np.zeros([nimgs,ny,nx], dtype = np.float32)
    for i in range(nimgs):
        dark_img_cube[i,:,:] = csusb.fitsread(dark_img_fnames[i])

    # Read flat images and generate an image cube
    print "  Reading Flat Fields"
    img_header = csusb.fitshead(flat_img_fnames[0])
    nx, ny = img_header['NAXIS1'], img_header['NAXIS2']
    nimgs = len(flat_img_fnames)
    flat_img_cube = np.zeros([nimgs,ny,nx], dtype = np.float32)
    for i in range(nimgs):
        flat_img_cube[i,:,:] = csusb.fitsread(flat_img_fnames[i])
        
    # Generate master bias image
    print "  Generating master bias image"
    master_bias_image = csusb.masterbias(bias_img_cube)

    # Generate master dark image
    print "  Generating master dark image"
    master_dark_image = csusb.masterdark(dark_img_cube)

    # Create normalized flat field
    print "  Generating normalized flat field image"
    master_flat_image = csusb.masterflat(flat_img_cube, master_dark_image)
        
    nimgs = len(sci_img_fnames)
    for i in range(nimgs):
        sci_file = sci_img_fnames[i]
        
        print " Science image : ", sci_file

        sci_image, header = csusb.fitsread(sci_file, header = True)

        # Reduced the science frames
        sci_red_image = csusb.imreduce(sci_image, master_dark_image, master_flat_image)

        # Write the reduced and average FITS image
        if output != "":
            red_file = output + "_final.fits"
        else:
            red_file = os.path.splitext(sci_file)[0] + '_final.fits'
        
        if os.path.exists(red_file):
            os.remove(red_file)
        
        csusb.fitswrite(sci_red_image, red_file, header = header)

        print "  Reduced science image : ", red_file

    return


if __name__ == "__main__":
    usage = "Usage: python %prog [options] sci_image bias_image dark_image flat_image"
    description = "Description. Utility to reduce raw science CSUSB telescope images."
    parser = OptionParser(usage = usage, version = "%prog 0.1", description = description)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default = False,
                      help = "print result messages to stdout"
                      )
    parser.add_option("-q", "--quiet",
                    action="store_false", dest="verbose", default = True,
                    help = "don't print result messages to stdout"
                    )
    parser.add_option("-t", "--threshold", dest = "threshold",
                    action='store', metavar="THRESHOLD", help = "Threshold for normalized flatfields (default is 0.8)",
                    default = 0.8
                    )
    parser.add_option("-o", "--output", dest = "output",
                    action="store", metavar="OUTPUT", help = "Output file name",
                    default = ""
                    )                    
                                        
    (options, args) = parser.parse_args()  
    if len(args) != 4:
        parser.error("REDUCE: Incorrect number of arguments")
        
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output
 
    process(args[0], args[1], args[2], args[3], options.threshold, options.output)    

    # Reset verbosity
    if not options.verbose:
        sys.stdout = old_stdout