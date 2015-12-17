# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 22:03:50 2015

# Author        Jonas Hartmann @ Gilmour Group @ EMBL Heidelberg

# Version       0.9 (BETA)
                Please report bugs to jonas.hartmann@embl.de or open an issue
                on github at https://github.com/WhoIsJack/tissueRipper

# Date          13.12.2015

# Descript      Example code to illustrate usage of the tissueRipper module.
                Only functional if the appropriate example files (see below)
                are present in the working directory!

# Usage         Required example files:
                    - TR_EXAMPLE.tif
                    - TR_EXAMPLE_segmentation.tif
                    - TR_EXAMPLE_green.tif
                    - TR_EXAMPLE_red.tif
                These files can be found at github (/WhoIsJack/tissueRipper).
                They roughly simulate a tissue consisting of 5 cells, with a
                full cell segmentation stack, a stack of a membrane marker and 
                a stack of some internal spots (e.g. from single molecule 
                fluorescence in-situ hybridisation, smFISH). 
                
# REQUIRES      Python 2.7, NumPy 1.9 (or similar), SciPy 0.15 (or similar),
                scikit-image 0.11.2 (or similar), tifffile 0.3.1 (or similar)
                
# License       None. You are free to use, change and redistribute this code 
                as you like. Mention of the original author is encouraged, 
                but not required.

# Warranty      This code is supplied "as is", without warranties or 
                conditions of any kind. The author takes no responsibility 
                for any inconvenience, damage or loss of data that may be 
                incurred as a result of using this code.
"""

#------------------------------------------------------------------------------

# PREPARATIONS AND IMPORTS

from __future__ import division, print_function
import numpy as np
import tissueRipper as tR


#------------------------------------------------------------------------------

# EXAMPLE 1: MULTI-CHANNEL INPUT AND ALTERNATIVE WAYS OF RUNNING THE RIPPER

# User specifications
input_filename = 'TR_EXAMPLE.tif'
output_filename = 'TR_EXAMPLE_RIPPED.tif'
scale_factor = 1.5 

# Import data using input function
# Note: Data can be imported in any other way (e.g. if skimage.io is unavailable)
data_stack = tR.input_from_tif(input_filename,channel_num=3)
print(np.shape(data_stack))

# Expand the tissue
data_out = tR.rip_tissue(data_stack,scale_factor)

# Alternatively: Do the two steps individually
#cell_labels, translocation = TR.translocation_from_centroids(data_stack,scale_factor)   # Get centroid translocations
#data_out = TR.translocate_pixels(data_stack,scale_factor,cell_labels,translocation)     # Translocate pixels

# Save the data using output function
# Note: Data can be exported in any other way (e.g. if tifffile.imsave is unavailable)
tR.output_to_tif(output_filename,data_out)


#------------------------------------------------------------------------------

# EXAMPLE 2: MULTIPLE SINGLE-CHANNEL INPUTS AND ALTERNATIVE WAYS OF SAVING THE OUTPUT

# User specifications
in_main_filename = 'TR_EXAMPLE_segmentation.tif'
in_other_filenames = ['TR_EXAMPLE_green.tif','TR_EXAMPLE_RED.tif']
output_filename = 'TR_EXAMPLE_RIPPED2.tif'
scale_factor = 1.3

# Import data using input function
data_stack = tR.input_from_tif(in_main_filename,additional_channels=in_other_filenames)

# Expand the tissue
data_out = tR.rip_tissue(data_stack,scale_factor,be_verbose=True)

# Save the data using output function
tR.output_to_tif(output_filename,data_out)

# Alternatively: Take the channels apart again and save them separately
#outfile_names = ['TR_EXAMPLE_segmentation_ripped.tif', 'TR_EXAMPLE_green_ripped.tif', 'TR_EXAMPLE_red_ripped.tif']
#from tifffile import imsave
#for index, filename in enumerate(outfile_names):
#    imsave(filename, data_out[index,:,:,:], bigtiff=True)


#------------------------------------------------------------------------------



