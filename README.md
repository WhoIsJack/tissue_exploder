# tissueRipper
A simple algorithm to create a nice "expanded" visualization of 3D single-cell segmentations of tissues imaged with confocal fluorescence microscopy.

## Description
- tissueRipper is	an algorithm that transforms 3D single-cell segmentations of confocal fluorescence microscopy data in a very particular way for the purpose of visualization. 
- All it really does is shift the segmented cells apart such that cell size and shape as well as overall tissue architecture are preserved but all cells can be seen individually. 
- tissueRipper is	written in Python 2.7 and requires the packages NumPy and SciPy (and optionally scikit-image and tifffile for import from and export to tif files).
- Note: tissueRipper is	not	a useful intermediate step in a quantitative analysis pipeline; it is well suited for showing off the power of 3D single-cell segmentation, but it does not itself improve segmentations or help  with the analysis of segmented data, so it is purely a visualization aid.

## Files
- tissueRipper.py           -- Contains all the code needed to actually run the algorithm. You can import this as a module and incorporate its functions into your image analysis pipeline.
- tR_example.py             -- Example code showing how to use the functions in tissueRipper.py.
- Example Images.zip        -- Contains example images that can be used with tR_example.py to test if everything is working. These images are (obviously) not real data, just a quickly made-up construction for testing and illustration.
- tissueRipper miniDoc.pdf  -- A very brief and general explanation of what the algorithm is, how it works, and how it is meant to be used in the context of image analysis pipelines. Does not include information on how to do things on the coding level; see tR_example.py for this purpose.
- tissueRipper-0.9.zip      -- Source distribution for installation. Unpack and run "python setup.py install" in the unpacked tissueRipper-0.9 directory to install the module.
- tissueRipper-0.9.win-amd64.exe -- Windows x64 installer making installation even easier if you happen to be working on this operating system.

## Installation
You have several options:
- Download tissueRipper.py and put it in the working directory you want to use it in. Not a real installation but works fine.
- Download tissueRipper-0.9.zip, unpack it and run "python setup.py install" in the unpacked directory.
- For windows x64 only -- Download tissueRipper-0.9.win-amd64.exe, run it and follow the installation wizard.
