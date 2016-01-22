#!/usr/bin/env python

"""
    --------------------------------------------------------------------------
    Routine to perform fast aperture photometry on CSUSB science frames.
    
    Usage: python fastphot.py [options] image coords
        
                                
    Author:
        Navtej Saini

    Organization:
        Caltech, Pasadena, CA, USA

    Version:
        21 January 2016     0.1     Initial implementation 
    --------------------------------------------------------------------------        
"""

import os, sys, dateutil
import numpy as np, warnings
from StringIO import StringIO
from optparse import OptionParser

try:
    import matplotlib.pylab as plt
except ImportError:
    plot_flag = False
else:
    try:
        import seaborn
    except ImportError:
        pass
    plot_flag = True


import csusb


def plotter(phot_data, ts, outfile):
    """
    Plot light curve. 
    
    Parameters
    ----------
    phot_data : numpy array
        Photometry array
        
    nframes : int
        Number of image cube frames
        
    exptime : float
        Kinetic or accumulation time
        
    outfile : string
        Name of the out png image
        
    Returns
    -------
    None
    """   
    params = {'backend': 'ps',
	      'font.size': 10,
              'axes.labelweight': 'medium',
	      'figure.dpi' : 300,
              'savefig.dpi': 300,
              'savefig.jpeg_quality': 100
              }
    plt.rcParams.update(params)
	   
    plt.figure(figsize=(6,4))
    plt.title("Normalized Light Curve : %s" %phot_data[0]['DATETIME'].split('T')[0])
    plt.xlabel("Time (secs)")
    plt.ylabel("Normalized Flux")
    plt.plot(ts, phot_data['FLUX_ADU']/np.mean(phot_data['FLUX_ADU']), "r-")    
    plt.savefig(outfile, dpi = 300, bbox_inches = "tight")
    
    return

                
def timedelta(dt1, dt2):
    """
    Determine datetime difference between two datetimes.
    
    Parameters
    ----------
    dt1, dt : string
        Datetime in iso string format
        
    Returns
    -------
    deltasec : int
        Time difference in seconds
    """            
    # Convert string iso format to python datetime object
    dt1 = dateutil.parser.parse(dt1)
    dt2 = dateutil.parser.parse(dt2)
                                
    # Time delta between two in seconds
    deltatime = dt2 - dt1                          
               
    return deltatime.seconds
                                                                                               
                                                                                                                                                                                                                                                               
def process(infile, coords, method, inner_radius, outer_radius, cen_method, window_size, output):
    """
    Entry point function to process science image.
    
    Parameters
    ----------
    infile : string
        Science image or list of science images
        
    coords : string
        Input text file with coordinates of stars
        
    method : string
        FWHM of the stelar psf in pixels
        
    inner_radius : float
        Sky background sigma
        
    outer_radius : int
        Inner sky annulus radius in pixels
        
                
    Returns
    -------
    None 
    """
    print "FASTPHOT: CHIMERA Fast Aperture Photometry Routine"
    
    inner_radius = float(inner_radius)
    outer_radius = float(outer_radius)
    
    # Check if input is a string of FITS images or a text file with file names
    if infile[0] == "@":
        infile = infile[1:]
        
        if not os.path.exists(infile):
            print "REGISTER: Not able to locate file %s" %infile
        
        images = []
        with open(infile, "r") as fd:
            for line in fd.readlines():
                if len(line) > 1:
                    images.append(line.replace("\n", ""))
    else:
        images = infile.split(",")

    # Number of images
    nimgs = len(images)

    dtype = [("DATETIME", "S25"),("XCEN", "f4"),("YCEN", "f4"),("MSKY", "f8"),("NSKY", "f8"),("AREA", "f8"),("FLUX_ADU", "f8"),("FLUX_ELEC", "f8"),("FERR", "f8"),("MAG", "f8")]
    phot_data = np.zeros([nimgs], dtype = dtype)
    for i in range(nimgs):
        sci_file = images[i]
        print "  Processing science image %s" %sci_file

        # Read FITS image and star coordinate
        image = csusb.fitsread(sci_file)
        pos = np.loadtxt(coords, ndmin = 2)

        # Instantiate an Aperphot object
        ap = csusb.Aperphot(sci_file, coords)
        
        # Set fwhmpsf, sigma, annulus and dannulus
        ap.method = method
        ap.inner_radius = inner_radius
        ap.outer_radius = outer_radius
        
        # Determine nominal aperture radius for photometry
        if i == 0:
            nom_aper = ap.cog(window_size, cen_method)
         
        print "  Nominal aperture radius : %4.1f pixels" %nom_aper
           
        # Perform aperture photometry on all the frames
        objpos = csusb.recenter(image, pos, window_size, cen_method)
        aperphot_data = ap.phot(image, objpos, nom_aper)
        pos = np.copy(objpos)
            
        phot_data[i]['DATETIME'] = ap.utcstart
        phot_data[i]['XCEN'] = aperphot_data["xcenter_raw"]
        phot_data[i]['YCEN'] = aperphot_data["ycenter_raw"]
        phot_data[i]['MSKY'] = aperphot_data["msky"]
        phot_data[i]['NSKY'] = aperphot_data["nsky"]
        phot_data[i]['AREA'] = aperphot_data["area"]
        phot_data[i]['FLUX_ADU'] = aperphot_data["flux"]
        phot_data[i]['FLUX_ELEC'] = phot_data[i]['FLUX_ADU'] * ap.epadu
        phot_data[i]['MAG'] = ap.zmag - 2.5 * np.log10(phot_data[i]['FLUX_ELEC']/ap.exptime)
            
        # Calculate error in flux - using the formula
        # err = sqrt(flux * gain + npix * (1 + (npix/nsky)) * (flux_sky * gain + R**2))
        phot_data[i]['FERR'] = np.sqrt(phot_data[i]['FLUX_ELEC'] + phot_data[i]['AREA'] * (1 + phot_data[i]['AREA']/phot_data[i]['NSKY']) * (phot_data[i]['MSKY'] * ap.epadu + ap.readnoise**2))
                        
    # Save photometry data in numpy binary format
    print "  Saving photometry data as numpy binary"
    if output != "":
        npy_outfile = output + ".npy"
    else:
        npy_outfile = sci_file.replace(".fits", ".phot.npy")
    if os.path.exists(npy_outfile):
        os.remove(npy_outfile)
            
    np.save(npy_outfile, phot_data)

    # Plot first pass light curve
    if plot_flag:
        print "  Plotting normalized light curve"
        if output != "":
            plt_outfile = output + ".png"
        else:
            plt_outfile = sci_file.replace(".fits", ".lc.png")
        ts = np.zeros(nimgs, dtype = np.int32)
        for i in range(nimgs):
            ts[i] = timedelta(phot_data[0]['DATETIME'], phot_data[i]['DATETIME'])    
        plotter(phot_data, ts, plt_outfile)
        
    return


if __name__ == "__main__":
    usage = "Usage: python %prog [options] sci_image coords"
    description = "Description. Utility to perform fast aperture photometry in CHIMERA science images."
    parser = OptionParser(usage = usage, version = "%prog 0.1", description = description)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default = False,
                      help = "print result messages to stdout"
                      )
    parser.add_option("-q", "--quiet",
                    action="store_false", dest="verbose", default = True,
                    help = "don't print result messages to stdout"
                    )
    parser.add_option("-m", "--method", dest = "method",
                    action="store", metavar="METHOD", help = "Method to use for determining overlap between aperture and pixels (default is exact)",
                    default = "exact"
                    )
    parser.add_option("-i", "--inner_radius", dest = "inner_radius",
                    action="store", metavar="INNER_RADIUS", help = "Inner radius of sky annlus in pixels (default is 20)",
                    default = 20
                    )
    parser.add_option("-d", "--outer_radius", dest = "outer_radius",
                    action="store", metavar="OUTER_RADIUS", help = "Radius of sky annulus in pixels (default is 40)",
                    default = 40
                    )
    parser.add_option("-c", "--cen_method", dest = "cen_method",
                    action="store", metavar="CEN_METHOD", help = "Centroid method (default is 2dg)",
                    default = "2dg"
                    )
    parser.add_option("-w", "--window_size", dest = "window_size",
                    action="store", metavar="WINDOW_SIZE", help = "Window size for centroid (default is 45)",
                    default = 45
                    )
    parser.add_option("-o", "--output", dest = "output",
                    action="store", metavar="OUTPUT", help = "Output file name",
                    default = ""
                    )                                                                                                       
                                        
    (options, args) = parser.parse_args()  
    if len(args) != 2:
        parser.error("FASTPHOT: Incorrect number of arguments")
        
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output
 
    # Switch off warnings
    warnings.filterwarnings('ignore')
    
    process(args[0], args[1], options.method, options.inner_radius, options.outer_radius, options.cen_method, options.window_size, options.output)    

    # Reset verbosity
    if not options.verbose:
        sys.stdout = old_stdout
