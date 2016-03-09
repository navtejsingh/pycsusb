
import csusb
import os, pyfits
import numpy as np
import ConfigParser

    
class Aperphot:
    def __init__(self, sci_file, coords):
        self.sci_file = sci_file
        self.coords = coords
        
        # Set header keyword parameters
        self.setkeywords()
        

    def readconfig(self):
        """
        Read configuration file.
        
        Parameters
        ----------
        
        Returns
        -------        
        """
        Config = ConfigParser.ConfigParser()
        Config.read(os.path.join(csusb.pkg_dir, "csusb.cfg"))
        
        return Config
        

    def setkeywords(self):
        """
        Set FITS image header keyword parameters.
        
        Parameters
        ----------
        
        Returns
        -------
        None
        """
        header = pyfits.getheader(self.sci_file, ignore_missing_end = True)
        self.nx = header["NAXIS1"]
        self.ny = header["NAXIS2"]
        self.exptime = header["EXPOSURE"]
        self.utcstart = header["DATE-OBS"]

        # Read detector and phot parameters from config file
        settings = self.readconfig() 
        self.zmag = float(settings.get('Phot', 'Zmag'))
        self.epadu = float(settings.get('Detector', 'Gain'))
        self.readnoise = float(settings.get('Detector', 'ReadNoise'))
        
    def cog(self, window_size, method, tolerance = 0.01):
        """
        Curve of growth to determine nominal aperture for photometry using 
        astropy photutils.
        
        Parameters
        ----------
        tolerance : float
            Magnitude difference tolerance between different apertures
        
        Returns
        -------
        aperture : float
            Nominal aperture radius for photmetry
        """
        # Aperture values in pixels
        apertures = np.linspace(2,25,24)
        naper = apertures.shape[0]
        
        # Read input image and star position
        image = csusb.fitsread(self.sci_file)
        pos = np.loadtxt(self.coords, ndmin = 2)
        nstars = pos.shape[0]
        
        # Iterate through the frames and determine nominal aperture
        mags_arr = np.zeros(len(apertures))
        objpos = csusb.recenter(image, pos, window_size, method)
        for i in range(naper):
            flux = self.phot(image, objpos[0,:], aper = apertures[i])
            mags_arr[i] = -2.5 * np.log10(flux['flux'])
        mags_diff = np.diff(mags_arr)
        idx = np.where(np.abs(mags_diff) < 0.01)
        if len(idx[0]) != 0:
            nom_aper = apertures[idx[0][0]]
        else:
            nom_aper = 16.0
            
        return nom_aper
        
                
        
    def phot(self, image, objpos, aper):
        """
        Aperture photometry using Astropy's photutils.
        
        Parameters
        ----------
        image : numpy array
            2D image array
            
        objpos : list of tuple
            Object poistions as list of tuples
            
        aper : float
            Aperture radius in pixels
         
        Returns 
        -------
        phot_table : astropy table
             Output table with stellar photometry   
        """
        try:
            from astropy.table import hstack
            from photutils import aperture_photometry, CircularAnnulus, CircularAperture
        except ImportError:
            pass
    
        apertures = CircularAperture(objpos, r = aper) 
        annulus_apertures = CircularAnnulus(objpos, r_in = self.inner_radius, r_out = self.outer_radius)
        
        rawflux_table = aperture_photometry(image, apertures = apertures, method = self.method)
        bkgflux_table = aperture_photometry(image, apertures = annulus_apertures, method = self.method)
        phot_table = hstack([rawflux_table, bkgflux_table], table_names = ["raw", "bkg"])
        
        bkg = phot_table["aperture_sum_bkg"] / annulus_apertures.area()
        phot_table["msky"] = bkg
        phot_table["area"] = apertures.area()
        phot_table["nsky"] = annulus_apertures.area()
                
        bkg_sum = bkg * apertures.area()
        final_sum = phot_table["aperture_sum_raw"] - bkg_sum
        phot_table["flux"] = final_sum
        
        return phot_table