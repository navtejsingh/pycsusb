import ez_setup
ez_setup.use_setuptools()

from setuptools import find_packages, setup

setup(
      name="CSUSB",
      version = "0.1",
      packages=find_packages(),
      
      # PROPER uses numpy and scipy
      install_requires = ['numpy>=1.8', 'scipy>=0.14', 'pyfits>=3.0', 'photutils>=0.2'],
      
      package_data = {
        # If any package contains *.txt, *.rst or *.fits files, include them:
        '': ['*.txt', '*.rst', '*.fits', "*.cfg"]
      },
      
      # Metadata for upload to PyPI
      author="Navtej Singh",
      author_email = "navtej@astro.caltech.edu",
      description="CSUSB telescope data processing pipeline",
      license = "BSD",
      platforms=["OSX", "Linux", "Unix"],
      url="",
) 
