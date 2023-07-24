import numpy as np
from skimage.measure import regionprops_table

def explode(segmentation, scale_factor=1.7, other_channels=None):
    """Explode a segmented tissue by scale_factor, shifting apart the cells but
    preserving cell size, shape and overall tissue organization.

    Parameters
    ----------
    segmentation : numpy array of shape ([z], y, x)
        A 2D or 3D segmentation image, where all pixels/voxels of a given cell
        share the same integer label.
    scale_factor : float, optional
        Factor by which the centroids of cells will be shifted apart.
        The default value, which is often quite good, is 1.7.
        Note that high values will lead to very large output images.
    other_channels : numpy array of shape (channels, [z], y, x), optional
        Other channels (e.g. fluorescence signal) that will be shifted apart
        in the same manner as the segmentation. Note that the 1st dimension
        must be the channels (so if only a single other channel is passed,
        its 1st dimension must be of size 1). The default is None.

    Returns
    -------
    new_segmentation : numpy array (2D or 3D)
        Exploded version of `segmentation`.
    new_other_channels : numpy array (3D or 4D)
        Exploded version of `other_channels`.
        Only returned if `other_channels` is not `None`.

    Notes
    -----
    This is a simplified and much faster version of the old implementation. It
    uses indexing arrays to vectorize any loops over individual cells. Note 
    that expanding this vectorization to time courses would be non-trivial, as 
    it isn't guaranteed that all centroids are present at all time points, so 
    the centroids arrays may no longer be dense.
    """

    # Find centroids
    old_centroids = regionprops_table(segmentation, properties=('centroid',))
    old_centroids = np.array([old_centroids['centroid-'+str(d)]
                              for d in range(len(old_centroids))])

    # Calculate xy translocation per centroid
    new_centroids = old_centroids * scale_factor
    translocation = (new_centroids - old_centroids).astype(int)

    # Initialize output arrays
    new_shape = tuple((np.array(segmentation.shape)*scale_factor+1).astype(int))
    new_segmentation = np.zeros(new_shape, dtype=segmentation.dtype)
    if other_channels is not None:
        new_other_channels = np.zeros((other_channels.shape[0], *new_shape),
                                      dtype=segmentation.dtype)

    # Construct indexing arrays
    # Note that `old_pointers` just indexes every pixel/voxel in `segmentation`
    # and `new_pointers` is the same but shifted by `translocation`, making use
    # of the fact that `segmentation-1` can index directly into `translocation`
    old_pointers = np.meshgrid(*[np.arange(s) for s in segmentation.shape],
                               indexing='ij')
    new_pointers = [old_pointers[d] + translocation[d, segmentation.astype(int)-1]
                    for d in range(segmentation.ndim)]

    # Exclude background from transfer
    old_pointers = tuple([op[segmentation>0] for op in old_pointers])
    new_pointers = tuple([np[segmentation>0] for np in new_pointers])

    # Translocate segmentation labels
    new_segmentation[new_pointers] = segmentation[old_pointers]

    # Translocate other channels
    if other_channels is not None:
        for ch in range(other_channels.shape[0]):
            new_other_channels[ch][new_pointers] = other_channels[ch][old_pointers]
            new_other_channels[ch][new_segmentation==0] = 0

    # Return results
    if other_channels is not None:
        return new_segmentation, new_other_channels
    else:
        return new_segmentation
