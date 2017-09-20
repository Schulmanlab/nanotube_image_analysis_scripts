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
from skimage import img_as_uint


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
cy3_file_list = os.listdir('cy3')
atto647_file_list = os.listdir('atto647')
atto488_file_list = os.listdir('atto488')
for i in range(len(cy3_file_list)):
	cy3_file = cy3_file_list[i]
	atto488_file = atto488_file_list[i]
	atto647_file = atto647_file_list[i]
	print "cy3 filename is "+str(cy3_file)
	image_unthresholded = io.imread('cy3/'+cy3_file)
	image_647_unthresholded = io.imread('atto647/'+atto647_file)
	#thresh = threshold_otsu(image_unthresholded)
	#image = image_unthresholded>thresh

	block_size = 15
	#image = threshold_local(image_unthresholded, block_size, offset=10)
	#image_647 = threshold_local(image_647_unthresholded, block_size, offset=10)

	radius = 5
	selem = disk(radius)

	#thresholding both files (getting rid of this because it should not be necessary!)
	#image = rank.otsu(image_unthresholded, selem)
	#image_647 = rank.otsu(image_647_unthresholded, selem)

	image = image_unthresholded
	image_647 = image_647_unthresholded


	#perfoming edge detection and morphological filling
	edges_open = canny(image, 3, 1, 25) #originally 2,1,25
	#edges_open = canny(image, 2) #originally 2,1,25
	selem = disk(5)
	edges = closing(edges_open, selem)
	fill_tubes = ndi.binary_fill_holes(edges)
	io.imsave(cy3_file+"fill_tubes.png", img_as_uint(fill_tubes), cmap=cm.gray)

	edges_open_647 = canny(image_647, 2, 1, 25)
	selem = disk(2)
	edges_647 = closing(edges_open_647, selem)
	fill_tubes_647 = ndi.binary_fill_holes(edges_647)
	io.imsave(atto647_file+"fill_tubes.png", img_as_uint(fill_tubes_647), cmap=cm.gray)

	#label image 
	label_image = label(fill_tubes)
	label_image_647 = label(fill_tubes_647)

	regions_joined_647 = []
	regions_joined_cy3 = []
	print "printing 647 components that are joined"
	print len(regionprops(label_image_647))
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
			regions_joined_cy3.append(joined_cy3_region)
			#print region.coords

	print "printing 647 components that are not joined"
	for region_647 in regionprops(label_image_647):
		if region_647 not in regions_joined_647 and region_647.area/2.2>= 3 and region_647.eccentricity>=.5:
			lengths_647_unjoined.append(region_647.area/2.2)

	print "printing cy3 components that are not joined"
	for region in regionprops(label_image):
		if region not in regions_joined_cy3 and region.area/2.2>= 3 and region.eccentricity>=.5:
			lengths_cy3_unjoined.append(region.area/2.2)

	i+=3

print "printing joined 647 lengths"
f1=open('joined647.dat','w+')
for length in lengths_647_joined:
	if length<=200:
		print >>f1, length
f1.close()

print "printing joined cy3 lengths"
f1=open('joinedcy3.dat','w+')
for length in lengths_cy3_joined:
	if length<=200:
		print >>f1, length
f1.close()

print "printing unjoined 647 lengths"
f1=open('unjoined647.dat','w+')
for length in lengths_647_unjoined:
	if length<=200:
		print >>f1, length
f1.close()

print "printing unjoined cy3 lengths"
f1=open('unjoinedcy3.dat','w+')
for length in lengths_cy3_unjoined:
	if length<=200:
		print >>f1, length
f1.close()
	#fill_tubes_skeleton = skeletonize(fill_tubes)

