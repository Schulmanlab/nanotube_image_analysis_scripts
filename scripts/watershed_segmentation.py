import numpy as np

from skimage.transform import (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
from skimage.feature import canny
from skimage import data
from skimage import io 
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import path

from scipy import ndimage as ndi
from skimage.filters import sobel, threshold_yen
from skimage.morphology import watershed
from skimage.measure import label, regionprops

#plt.ion()

image = io.imread(str(sys.argv[1]))

#thresh = threshold_yen(image_unthresholded)
#image = image_unthresholded > thresh

#plt.hist(image, bins=np.arange(0, 256))
#plt.savefig("grey_histo.pdf")
#plt.show()
#plt.close()
elevation_map = sobel(image)
#io.imsave("elevation_map_test.jpeg", elevation_map)

markers = np.zeros_like(image)
markers[image < 80] = 1
markers[image > 100] = 2 

#io.imsave("markers_test.jpeg", markers)

segmentation = watershed(elevation_map, markers)

segmentation = ndi.binary_fill_holes(segmentation - 1)

print segmentation

labelled_tubes, num_tubes = ndi.label(segmentation)

print num_tubes
print len(labelled_tubes)

'''for region in regionprops(label_image):
	if region.area/2.2>=3:
		print region.area/2.2
		print "eccentricity "+ str(region.eccentricity)
'''
plt.imshow(segmentation)
plt.show()



