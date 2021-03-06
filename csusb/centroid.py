from __future__ import division


import numpy as np
from photutils.morphology import centroid_com, centroid_1dg, centroid_2dg


def recenter(image, pos, window_size = 15, method = "2dg", threshold = 10.):
    """
    Recenter each star in each frame of the image cube before performing
    aperture photometry to take care of slight misalignments between frames
    because of atmospheric turbulence and tracking/pointing errors.
        
    Parameters
    ----------
    image : numpy array
        2D image
        
    pos : list
        List of (x,y) tuples for star positions
            
    window_size : int
        Window size in which to fit the gaussian to the star to calculate
        new center
        
    method : string
        Method used to find center of the star. Options are 1d Gaussian fit,
        2d gaussian fit or com (center of mass)    

    threshold : int
	Threshold of star movement in pixels. If shift between old and new
	coordinates is > threshold, old value is assigned to new coordinates.	
                    
    Returns
    -------
    xcen, ycen : float
        Source x and y centers               
    """        
    pos = np.asarray(pos)
    
    ny, nx = image.shape        
    window_size = int(window_size)
    threshold = float(threshold)
    nstars = pos.shape[0]
        
    star_pos = np.zeros([nstars,2], dtype = np.float32)
    for i in range(nstars):
        x, y = pos[i][0], pos[i][1]
            
        xmin, xmax = int(x) - int(window_size/2), int(x) + int(window_size/2) + 1
        ymin, ymax = int(y) - int(window_size/2), int(y) + int(window_size/2) + 1
        
        if xmin < 0:     
            xmin = 0
                
        if ymin < 0:
            ymin = 0
                
        if xmax > nx:
            xmax = nx
                
        if ymax > ny:
            ymax = ny                                                                         
                
        if method == "1dg":
            xcen, ycen = centroid_1dg(image[ymin:ymax,xmin:xmax])
        elif method == "2dg":
            xcen, ycen = centroid_2dg(image[ymin:ymax,xmin:xmax])
        elif method == "com":
            xcen, ycen = centroid_com(image[ymin:ymax,xmin:xmax])
        
        if (np.abs(xmin + xcen - x)) > threshold or (np.abs(ymin + ycen - y)) > threshold:
            star_pos[i,0] = x
            star_pos[i,1] = y
        else:
            star_pos[i,0] = xmin + xcen
            star_pos[i,1] = ymin + ycen 

    return star_pos
