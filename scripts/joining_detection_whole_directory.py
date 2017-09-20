import numpy as np

from skimage.transform import (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
from skimage.feature import canny
from skimage import data
from skimage import io 
from skimage.morphology import closing, disk
from skimage.morphology import skeletonize
from skimage.measure import label, regionprops
from skimage.filters import threshold_otsu, threshold_local, rank
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import path


from scipy import ndimage as ndi

def line_length(line):
	p0, p1 = line
	a = np.array((p0[0],p0[1]))
	b = np.array((p1[0],p1[1]))
	dist = np.linalg.norm(a-b)
	#print dist
	return dist

def generate_rectangles(lines, width):
	#for each line, define a rectange that is the length of the line and has a width of 10 (5 on either side)
	rectangles = []
	for line in lines:
		p0, p1 = line
		a = np.array((p0[0],p0[1],0))
		b = np.array((p1[0],p1[1],0))
		z_vec = np.array((0,0,1))
		slope_vector = a - b
		perpendicular_to_slope = np.cross(slope_vector,z_vec)
		normalized_perpendicular_to_slope = perpendicular_to_slope / np.linalg.norm(perpendicular_to_slope)
		#first 2 corners
		c1 = a + width*normalized_perpendicular_to_slope
		c2 = a - width*normalized_perpendicular_to_slope
		#second 2 corners 
		c3 = b - width*normalized_perpendicular_to_slope
		c4 = b + width*normalized_perpendicular_to_slope



		rectangle = np.array([[c1[0],c1[1]],[c2[0],c2[1]],[c3[0],c3[1]],[c4[0],c4[1]]])
		rectangle_path = path.Path(rectangle) #converting vertices into path so that contains_point can be used
		#print rectangle_path

		rectangles.append(rectangle_path)
	return rectangles

#def generate_rectangle_from_bounding_box


# Line finding using the Probabilistic Hough Transform
lengths_cy3_joined = []
lengths_647_joined = []
lengths_cy3_unjoined = []
lengths_647_unjoined = []

i=0
while i < len(os.listdir(os.getcwd())):
	file_list = os.listdir(os.getcwd())
	cy3_file = file_list[i]
	print "cy3 filename is "+str(cy3_file)
	atto488_file = file_list[i+1]
	atto647_file = file_list[i+2]
	image_unthresholded = io.imread(cy3_file)
	image_647_unthresholded = io.imread(atto647_file)
	#thresh = threshold_otsu(image_unthresholded)
	#image = image_unthresholded>thresh

	#block_size = 35
	#image = threshold_local(image_unthresholded, block_size, offset=10)

	radius = 3
	selem = disk(radius)

	#thresholding both files
	image = rank.otsu(image_unthresholded, selem)
	image_647 = rank.otsu(image_647_unthresholded, selem)


	#perfoming edge detection and morphological filling
	edges_open = canny(image, 2, 1, 25)
	selem = disk(1.5)
	edges = closing(edges_open, selem)
	fill_tubes = ndi.binary_fill_holes(edges)

	edges_open_647 = canny(image_647, 2, 1, 25)
	selem = disk(1.5)
	edges_647 = closing(edges_open_647, selem)
	fill_tubes_647 = ndi.binary_fill_holes(edges_647)

	#label image 
	label_image = label(fill_tubes)
	label_image_647 = label(fill_tubes_647)

	regions_joined_647 = []
	print "printing 647 components that are joined"
	for region_647 in regionprops(label_image_647):
		is_joined = 0
		if region_647.area/2.2 >= 3 and region_647.eccentricity >=.5:
			for region in regionprops(label_image):
				if region.area <.3 or region.eccentricity < .5:
					continue
				region_647_coords = region_647.coords.tolist()
				region_coords = region.coords.tolist()
				if any(i in region_647_coords for i in region_coords):
					is_joined = 1
					#print "we have joining!"
					break
		if is_joined == 1:
			lengths_647_joined.append(region_647.area/2.2)
			regions_joined_647.append(region_647)
#print regions_joined_647
#print regionprops(label_image_647)

	print "printing cy3 components that are joined"
	for region_647 in regionprops(label_image_647):
		is_joined = 0
		if region_647.area/2.2 >= 3 and region_647.eccentricity >=.5:
			for region in regionprops(label_image):
				if region.area <.3 or region.eccentricity < .5:
					continue
				region_647_coords = region_647.coords.tolist()
				region_coords = region.coords.tolist()
				if any(i in region_647_coords for i in region_coords):
					is_joined = 1
					joined_cy3_region = region
					break
		if is_joined == 1:
			lengths_cy3_joined.append(joined_cy3_region.area/2.2)
			#print region.coords

	print "printing 647 components that are not joined"
	for region_647 in regionprops(label_image_647):
		if region_647 not in regions_joined_647 and region_647.area/2.2>= 3 and region_647.eccentricity>=.5:
			lengths_647_unjoined.append(region_647.area/2.2)

	i+=3

print "printing joined 647 lengths"
for length in lengths_647_joined:
	print length

print "printing joined cy3 lengths"
for length in lengths_cy3_joined:
	print length

print "printing unjoined 647 lengths"
for length in lengths_647_unjoined:
	print length
	#fill_tubes_skeleton = skeletonize(fill_tubes)

