# -*- coding: utf-8 -*-
"""
Created on Wed Jul 01 15:11:57 2015

# Author        Jonas Hartmann @ Gilmour Group @ EMBL Heidelberg

# Version       0.9 (BETA)
                Please report bugs to jonas.hartmann@embl.de or open an issue
                on github at https://github.com/WhoIsJack/tissueRipper

# Date          13.12.2015

# Descript      "Ripped Tissue Visualization"
                Trick for visualizing segmentations of 3D imaging data of
                biological tissues; shifts apart the segmented cells by a 
                given scale factor, preserving overall cell shapes and tissue
                structure, yet showing the cells individually.
                Requires a 3D matrix of labeled entities as input (i.e. a
                segmented confocal microscopy stack), returns an appropriately
                resized matrix with shifted entities. Other channels (such as
                the original unsegmented intensity data) can also be shifted
                accordingly.
                
# Usage         A usage example can be found in tR_example.py on github:
                https://github.com/WhoIsJack/tissueRipper
                This example can also be used to test if everything works as
                intended before you try running your own data.
                
                The key function in this module is "rip_tissue", which takes an 
                input numpy array and a scale factor (which indicates by how 
                much the entire matrix should be expanded) and returns the 
                expanded ("ripped") matrix. See help(rip_tissue) for more usage
                information.
                
                Alternatively, the functions "translocation_from_centroids" and 
                "translocate_pixels" can be called separately, for example if
                further modifications on their output and input are desired.
                
                For import of tif stacks to the required numpy array format
                and for export of the result into tif stacks, the additional
                functions "input_from_tif" and "output_to_tif" are available.
                They require the additional modules scikit-image and tifffile,
                respectively (see # REQUIRES).
                
# REQUIRES      Python 2.7, NumPy 1.9 (or similar), SciPy 0.15 (or similar),
                scikit-image 0.11.2* (or similar), tifffile 0.3.1* (or similar)
                  *only required when importing/saving from/to tif stacks!
                
# License       None. You are free to use, change and redistribute this code 
                as you like. Mention of the original author is encouraged, 
                but not required.

# Warranty      This code is supplied "as is", without warranties or 
                conditions of any kind. The author takes no responsibility 
                for any inconvenience, damage or loss of data that may be 
                incurred as a result of using this code.
                
# Notes         (1) It's possible that the cell pixel translocation loop could
                be vectorized to speed up the program dramatically. I might
                look into this at some point, but I can't think of a simple  
                solution right off the bat, so I'll leave it for now.
"""

#------------------------------------------------------------------------------

# PREPARATIONS AND IMPORTS

from __future__ import division, print_function
import numpy as np
from warnings import warn


#------------------------------------------------------------------------------

# FUNCTION TO CALCULATE TRANSLOCATION BASED ON CENTROIDS

def translocation_from_centroids(seg_stk, scale_factor, return_centroids=False):
    """Help on translocation_from_centroids:
    
    translocation_from_centroids(seg_stk, scale_factor, return_centroids=False)
    
        Calculate by how many pixels the centroids of each segmented cell in 
        the seg_stk array must be shifted in x, y and z in order to expand the
        tissue by scale_factor.
        
        Parameters
        ----------
        seg_stk: numpy array (3D or 4D)
            Array to be expanded. If 3D, the array should contain a labeled
            segmentation only. If 4D, dimension 0 should be the channels, with
            the first channel containing the segmentation.      
        scale_factor: int or float
            Factor by which the tissue will be expanded. For example, an input 
            array of shape (100,100,100) with a scale_factor of 1.3 will result
            in an output array of shape (130,130,130). Note that the output
            array is a factor 2 larger than the input; higher scale factors
            will therefore rapidly increase output array size, which may lead
            to memory issues!
        return_centroids: bool (default: False)
            If true, the centroid coordinates, both in the old and the new 
            (shifted) array will be returned. This can be used to make further
            modifications as needed before translocating the pixels.
            
        Returns
        -------
        cell_labels: numpy array (1D)
            Enumeration of the labels of all labeled entities in the seg_stk
            array. Useful to iterate over all cells.
        translocation: numpy array (Nx3)
            Array containing the number of pixels each of the N cell centroids 
            needs to be shifted in x, y and z.
        old_centroids: numpy array (Nx3)
            Only returned if return_centroids==True. Pixel coordinates of the
            centroids for each of the N cells in the unexpanded condition.
        new_centroids: numpy array (Nx3)
            Only returned if return_centroids==True. Pixel coordinates of the
            centroids for each of the N cells in the new, expanded condition.
            
        Requires
        --------
        NumPy 1.9 or similar
        SciPy 0.15 or similar
    """
        
    # Remove channel dimension for 4D arrays; work on segmentation only
    if len(seg_stk.shape) == 4:
        seg_stk = seg_stk[0,:,:,:]
    
    # Create list of cell labels [assume 0 to be background]
    cell_labels = np.unique(seg_stk)[1:]    
    
    # Grab centroids
    from scipy.ndimage import measurements
    old_centroids = np.array(measurements.center_of_mass(seg_stk,labels=seg_stk,index=cell_labels))
    
    # Calculate new centroid positions using scale_factor
    new_centroids = old_centroids * scale_factor

    # Calculate centroid translocation; all pixels belonging to each centroid
    #   will be translocated by the appropriate amount
    translocation = (new_centroids - old_centroids).astype(np.int32)

    # Return results
    if return_centroids:
        return cell_labels, translocation, old_centroids, new_centroids
    else:
        return cell_labels, translocation


#------------------------------------------------------------------------------

# FUNCTION TO TRANSLOCATE ALL CELL PIXELS

def translocate_pixels(seg_stk, scale_factor, cell_labels, translocation, be_verbose=False):
    """Help on translocate_pixels:
    
    translocate_pixels(seg_stk, scale_factor, cell_labels, translocation, be_verbose=False)
    
        Returns the expanded array ripped_seg where all labeled entities have
        been shifted apart by scale_factor. The shift is done based on pixel
        translocation values which are calculated using the cell's centroids
        (see function translocation_from_centroids).
        
        Parameters
        ----------
        seg_stk: numpy array (3D or 4D)
            Array to be expanded. If 3D, the array should contain a labeled
            segmentation only. If 4D, dimension 0 should be the channels, with
            the first channel containing the segmentation.      
        scale_factor: int or float
            Factor by which the tissue will be expanded. For example, an input 
            array of shape (100,100,100) with a scale_factor of 1.3 will result
            in an output array of shape (130,130,130). Note that the output
            array is a factor 2 larger than the input; higher scale factors
            will therefore rapidly increase output array size, which may lead
            to memory issues!
        cell_labels: numpy array (1D)
            Enumeration of the labels of all labeled entities in the seg_stk
            array. Used to iterate over all cells.
        translocation: numpy array (Nx3)
            Array containing the number of pixels each of the N centroids needs
            to be shifted in x, y and z to achieve an expansion by scale_factor.
        be_verbose: bool (default: False)
            If true, more progress information about the running algorithm will
            be printed.
            
        Returns
        -------
        ripped_seg: numpy array (3D or 4D)
            Expanded ("ripped") input array with the same dimensions and 
            channel order.
        
        Requires
        --------
        NumPy 1.9 or similar
    """
    
    # Warn in case of high bitdepth
    if not (type(seg_stk.flatten()[0]) == np.uint8 or type(seg_stk.flatten()[0]) == np.int8):
        warn('It is recommended to use 8-bit arrays (np.uint8) to minimize memory consumption!')
    
    # Add 1-size dimension to seg_stk for 3D arrays
    if len(seg_stk.shape) == 3:
        seg_stk = np.expand_dims(seg_stk,0)
        
    # Initialize output segmentation array
    stk_dim = np.shape(seg_stk)
    ripped_seg = np.zeros((stk_dim[0],np.int(stk_dim[1]*scale_factor),np.int(stk_dim[2]*scale_factor),np.int(stk_dim[3]*scale_factor)),dtype=type(seg_stk.flatten()[0]))
    
    # Be verbose about progress
    if be_verbose:
        print("\nStarting loop...")
    
    # For each cell...
    for centroid_id,cell_id in enumerate(cell_labels):
    
        # Find all pixel positions in old stack
        old_pxl_indices = np.array(np.where(seg_stk[0,:,:,:]==cell_id))
        
        # Calculate new pixel positions using (broadcasted) centroid translocation
        cell_transloc = translocation[centroid_id].repeat(np.shape(old_pxl_indices)[1]).reshape(3,np.shape(old_pxl_indices)[1])
        new_pxl_indices = old_pxl_indices + cell_transloc
        
        # Assign actual values to the new pixels
        for channel in range(np.shape(seg_stk)[0]):      # For each channel...
            ripped_seg[channel,new_pxl_indices[0],new_pxl_indices[1],new_pxl_indices[2]] = seg_stk[channel,old_pxl_indices[0],old_pxl_indices[1],old_pxl_indices[2]]
        
        # Be verbose about progress
        if be_verbose:
            print("  Done with loop", centroid_id+1, "of", len(cell_labels))
    
    # Be verbose about progress
    if be_verbose:    
        print("Done with loop!")
    
    # Return results
    return ripped_seg
    

#------------------------------------------------------------------------------

# ASSEMBLED TISSUE "RIPPING" FUNCTION

def rip_tissue(seg_stk, scale_factor, be_verbose=False):
    """Help on rip_tissue:
    
    rip_tissue(seg_stk, scale_factor, be_verbose=False)
    
        Expand a segmented tissue by scale_factor, shifting apart the cells but
        preserving cell size, shape and overall tissue organization.
        
        Parameters
        ----------
        seg_stk: numpy array (3D or 4D)
            Array to be expanded. If 3D, the array should contain a labeled
            segmentation only. If 4D, dimension 0 should be the channels, with
            the first channel containing the segmentation.      
        scale_factor: int or float
            Factor by which the tissue will be expanded. For example, an input 
            array of shape (100,100,100) with a scale_factor of 1.3 will result
            in an output array of shape (130,130,130). Note that the output
            array is a factor 2 larger than the input; higher scale factors
            will therefore rapidly increase output array size, which may lead
            to memory issues!
        be_verbose: bool (default: False)
            If true, more progress information about the running algorithm will
            be printed.
            
        Returns
        -------
        ripped_seg: numpy array (3D or 4D)
            Expanded ("ripped") input array with the same dimensions and 
            channel order.
        
        Requires
        --------
        NumPy 1.9 or similar
        SciPy 0.15 or similar
    """
        
    # Be verbose about progress
    if be_verbose:
        print('\nCalculating translocations based on centroids')
        
    # Calculate translocation based on centroids
    cell_labels, translocation = translocation_from_centroids(seg_stk, scale_factor)
    
    # Be verbose about progress
    if be_verbose:
        print('\nTranslocating pixels')
        
    # Translocate all cell pixels
    ripped_seg = translocate_pixels(seg_stk, scale_factor, cell_labels, translocation, be_verbose=be_verbose)
    
    # Be verbose about progress
    if be_verbose:
        print('\nAll done; returning results')
        
    # Return results
    return ripped_seg


#------------------------------------------------------------------------------

# OPTIONAL INPUT LOADING FUNCTION FOR TIF FILES

def input_from_tif(filename, channel_num=1, additional_channels=None):
    """Help on input_from_tif:
    
    input_from_tif(filename, channel_num=1, additional_channels=None)
    
        Import a 3-dimensional tif file and convert it to a numpy array. The 
        stack may contain multiple channels reduced to 3D in the order CZXY. In
        this case, a 4D numpy array will be produced, with channels in the 0th
        dimension.
        Optionally, also import several additional 3-dimensional tif files and
        concatenate them with the first array into a larger 4D array. These 
        additional channels must be pure 3D (one channel only).
        
        Parameters
        ----------
        filename: string
            Name of the tif file to be imported. Must end on '.tif' or '.tiff'. 
            File must be in the working directory or filename must contain the
            full path. For multi-channel stacks reduced to 3D, the stack order
            must be CZXY. Timecourses are not permitted.
            In the context of tissueRipper, this should be the segmentation on
            its own (if it is 3D) or the segmentation plus other channels (if 
            it is reduced 4D).
        channel_num: int (default: 1)
            Number of channels reduced to 3D in the file to be imported. By
            default, a pure 3D stack (channel_num = 1) is assumed.          
        additional_channels: list of strings or None (default: None)
            List of filenames for additional 3-dimensional tif files to be
            imported and concatenated to the first imported array.
            Concatenation is done in dimension 0 (that is, as channels).
            In the context of tissue ripper, these could be further channels
            to be expanded ("ripped") in addition to the segmentation channel.
        
        Returns
        -------
        main_img: numpy array (3D or 4D)
            Imported image(s).
        
        Requires
        --------
        NumPy 1.9 or similar
        scikit-image 0.11.2 or similar
    """
    
    # Prep
    from skimage import io
    
    # SUBFUNCTION FOR IMPORT
    def import_tif(filename):
        
        # Check that an acceptable filename has been specified    
        if not (filename[-4:] == '.tif' or filename[-5:] == '.tiff'):
            raise NameError("INVALID FILE EXTENSION! Please specify .tif or .tiff in the filename!")
        # Check that an existing file has been specified
        from os.path import isfile
        if not isfile(filename):
            raise IOError("INPUT TIF FILE NOT FOUND! ("+filename+")")
        
        # If all is well, proceed to load the file using skimage.io and convert to numpy  
        imported_img = io.imread(filename)
        imported_img = np.array(imported_img)
        
        # Then, check for channel_num and split the array into 4D by channels if necessary
        if channel_num > 1:
            
            # Initialize 4D array of appropriate size
            img4d = np.zeros((channel_num,imported_img.shape[0]/channel_num,imported_img.shape[1],imported_img.shape[2]),dtype=type(imported_img.flatten()[0]))
            
            # Funky looping construct to transform   
            for i,j in enumerate(range(0,len(imported_img),channel_num)):
                for k in range(channel_num):
                    img4d[k,i,:,:] = imported_img[j+k,:,:]
            
            # Assign back in order to return the result
            imported_img = img4d
        
        # Return result
        return imported_img
    
    # Import the main file
    main_img = import_tif(filename)
    
    # If there are additional channels, import and add them
    if additional_channels:
        
        # If main_img is only one channel so far ("pure 3D"), add a fourth dimension for concatenation
        if len(main_img.shape) == 3:
            main_img = np.expand_dims(main_img,0)
        
        # For each additional channel filename, add the file
        for additional_filename in additional_channels:
            
            # Import data
            add_img = import_tif(additional_filename)
            
            # Concatenate, assume additional image to be 3D
            main_img = np.concatenate([main_img,np.expand_dims(add_img,0)])
    
    # Return the final stack
    return main_img
    

#------------------------------------------------------------------------------

# OPTIONAL OUTPUT WRITING FUNCTION FOR TIF FILES

def output_to_tif(filename, ripped_seg):
    """Help on output_to_tif:
    
    output_to_tif(filename, ripped_seg)
    
        Write a multidimensional numpy (image) array to a tif file.
        
        Parameters
        ----------
        filename: string
            Name of the file to be written. Must end on '.tif' or '.tiff'.
        ripped_seg: numpy array (3D or 4D)
            Multidimensional numpy image array. In the context of tissueRipper,
            this would be an extended ("ripped") segmentation, either 3D for 
            the segmentation only or 4D when other channels are included.
        
        Returns
        -------
        nothing: writes a 3D tif file
            If 4D data was given, it will be reduced to 3D in the resulting tif
            stack, following the order CZXY. 
            WARNING: If a file with an identical path/filename is present, it
            will be overwritten without prompt!
        
        Requires
        --------
        NumPy 1.9 or similar
        tifffile 0.3.1 or similar
    """

    # Prep
    from tifffile import imsave

    # Check that an acceptable filename has been specified    
    if not (filename[-4:] == '.tif' or filename[-5:] == '.tiff'):
        raise NameError("INVALID FILE EXTENSION! Please specify .tif or .tiff in the filename!")
    
    # For 4D, reduce to 3D with order CZXY
    if len(ripped_seg.shape) == 4:
        
        # Initialize the 3D version
        img3d = np.zeros((ripped_seg.shape[0]*ripped_seg.shape[1],ripped_seg.shape[2],ripped_seg.shape[3]),dtype=type(ripped_seg.flatten()[0]))
        
        # Loop to create the image
        for i,j in enumerate(range(0,img3d.shape[0],ripped_seg.shape[0])):
            for k in range(ripped_seg.shape[0]):
                img3d[j+k,:,:] = ripped_seg[k,i,:,:]
        
        # Assign back to save
        ripped_seg = img3d
    
    # Save
    print('\nSaving to ' + filename + '...')
    imsave(filename, ripped_seg, bigtiff=True)
    print('Done!')
    
    # Done
    return


#------------------------------------------------------------------------------



