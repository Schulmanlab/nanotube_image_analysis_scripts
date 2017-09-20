import numpy as np

from skimage.transform import (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
from skimage.feature import canny
from skimage import data
from skimage import io 
from skimage.morphology import closing, disk
from skimage.morphology import skeletonize
from skimage.measure import label, regionprops
from skimage.filters import threshold_otsu, threshold_local, rank, threshold_yen, threshold_adaptive
import sys

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
image_unthresholded = io.imread(str(sys.argv[1]))
#image_647_unthresholded = io.imread(str(sys.argv[2]))
#thresh = threshold_otsu(image_unthresholded)
#image = image_unthresholded>thresh

#block_size = 35
#image = threshold_local(image_unthresholded, block_size, offset=10)

radius = 3
selem = disk(radius)

'''#image=image_unthresholded
image = rank.otsu(image_unthresholded, selem)
#image_647 = rank.otsu(image_647_unthresholded, selem)
'''

thresh = threshold_yen(image_unthresholded)
image = image_unthresholded > thresh

#image = rank.yen(image_unthresholded, selem)
'''#edges_open = canny(image, 2, 1, 25)
edges_open = canny(image, 2, 1, 25)
selem = disk(1.5)
edges = closing(edges_open, selem)
fill_tubes = ndi.binary_fill_holes(edges)'''



#label image 
fill_tubes = image
print fill_tubes
label_image = label(fill_tubes)
print len(regionprops(label_image))
for region in regionprops(label_image):
	if region.area/2.2>=3:
		print region.area/2.2
		print "eccentricity "+ str(region.eccentricity)

'''print "printing A components not joined yet"
for region_647 in regionprops(label_image_647):
	is_joined = 0
	if region_647.area/2.2 >= 3:
		for region in regionprops(label_image):
			region_647_coords = region_647.coords.tolist()
			region_coords = region.coords.tolist()
			#print not any(i in region_647_coords for i in region_coords)
			print is_joined
			if not any(i in region_647_coords for i in region_coords):
				is_joined = 1
	if is_joined == 1:
		print region_647.area/2.2'''

'''regions_joined_647 = []

print "printing A components that are joined"
for region_647 in regionprops(label_image_647):
	is_joined = 0
	if region_647.area/2.2 >= 3:
		for region in regionprops(label_image):
			region_647_coords = region_647.coords.tolist()
			region_coords = region.coords.tolist()
			if any(i in region_647_coords for i in region_coords):
				is_joined = 1
				#print "we have joining!"
				break
	if is_joined == 1:
		print region_647.area/2.2
		regions_joined_647.append(region_647)
#print regions_joined_647
#print regionprops(label_image_647)

print "printing B components that are joined"
for region_647 in regionprops(label_image_647):
	is_joined = 0
	if region_647.area/2.2 >= 3:
		for region in regionprops(label_image):
			region_647_coords = region_647.coords.tolist()
			region_coords = region.coords.tolist()
			if any(i in region_647_coords for i in region_coords):
				is_joined = 1
				joined_B_region = region
				break
	if is_joined == 1:
		print joined_B_region.area/2.2
		#print region.coords

print "printing A components that are not joined"
for region_647 in regionprops(label_image_647):
	if region_647 not in regions_joined_647 and region_647.area/2.2>= 3:
		print region_647.area/2.2
#fill_tubes_skeleton = skeletonize(fill_tubes)'''

lines = probabilistic_hough_line(fill_tubes, threshold=10, line_length=3,
                                 line_gap=1)
#removing all lines that begin or end within a rectangle defined by a longer line with a width of 10
culled_lines = []
sorted_lines = sorted(lines, key = line_length)
sorted_rectangles = generate_rectangles(sorted_lines, 10.0)
#print sorted_rectangles

#for each line, determine whether either point lies within a rectangle defined by a longer line, if its 
#inside the rectangle, discard this line
'''for line in sorted_lines: 
	
	p0, p1 = line
	a = np.array((p0[0],p0[1],0))
	b = np.array((p1[0],p1[1],0))
	for rectangle in sorted_rectangles:
		idx_line = sorted_lines.index(line)
		idx_rectangle = sorted_rectangles.index(rectangle)
		if idx_line == idx_rectangle:
			continue
		keep_line = True
		print rectangle
		print a
		print b 
		print rectangle.contains_point(a)
		print rectangle.contains_point(b)

		if rectangle.contains_point(a) or rectangle.contains_point(b):
			print "discarding line"
			print rectangle
			print a
			print b 
			print rectangle.contains_point(a)
			print rectangle.contains_point(b)
			keep_line = False
			break
	print keep_line
	if keep_line:
		culled_lines.append(line)

print len(lines)
print len(culled_lines)'''

#lines = culled_lines


# Generating figure 2
'''fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(image, cmap=cm.gray)
ax[0].set_title('Input image')

ax[1].imshow(fill_tubes, cmap=cm.gray)
ax[1].set_title('Canny edges')

ax[2].imshow(edges * 0)
for line in lines:
    p0, p1 = line
    ax[2].plot((p0[0], p1[0]), (p0[1], p1[1]))
ax[2].set_xlim((0, image.shape[1]))
ax[2].set_ylim((image.shape[0], 0))
ax[2].set_title('Probabilistic Hough')

for a in ax:
    a.set_axis_off()
    a.set_adjustable('box-forced')

plt.tight_layout()'''
plt.imshow(image)
plt.show()