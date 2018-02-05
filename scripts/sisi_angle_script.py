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
import math
from matplotlib import cm
from matplotlib import path
from skimage import img_as_uint
from scipy import ndimage
from scipy.spatial import distance
from scipy import ndimage as ndi
from numpy import unravel_index
import Tkinter, tkFileDialog

#modifying the joining detection script to measure the angle of Sisi's nanotubes relative to the x-axis of her images 

#constants
#constants
tube_width = 5.0
length_cutoff = 3.0 
eccentricity_cutoff = 0.5
end_to_end_distance_cutoff = 10.0

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(abs(dotproduct(v1, v2) )/ (length(v1) * length(v2)))

def line_length(line):
	p0, p1 = line
	a = np.array((p0[0],p0[1]))
	b = np.array((p1[0],p1[1]))
	dist = np.linalg.norm(a-b)
	#print dist
	return dist

def make_endpoints_mask(filled_binary_image):
	#function to determine the endpoints of a nanotube identified via edge detection/morphological filling
	#need to find all endpoint candidates and find the pair separated by the longest path

	#first skeletonize the filled binary image (must be a binary int image)
	filled_binary_image = filled_binary_image.astype(int)
	skeleton = skeletonize(filled_binary_image)
	skeleton = skeleton.astype(int)

	#now we make a kernel to compute the endpoints of the skeletonized image
	kernel = np.uint8([[1,  1, 1], [1, 10, 1], [1,  1, 1]])

	#now we convolve the kernel with the skeletonized image 
	convolved_skeleton = ndimage.convolve(skeleton, kernel, mode='constant', cval = 1)

	#now produce an output mask with only pixels with value 11, these are the endpoints
	endpoint_mask = np.zeros_like(convolved_skeleton)
	endpoint_mask[np.where(convolved_skeleton == 11)] = 1

	return endpoint_mask

def endpoints(region_coords, endpoint_mask):
	#using a previously genereated endpoint mask to find the endpoints for a particular tube
	#this will return a pair of tubles with the x,y coordinates of the two endpoints 
    endpoints_labelled = label(endpoint_mask)
    potential_endpoints = []
    for endpoint in regionprops(endpoints_labelled):
    	if any(i in region_coords for i in endpoint.coords.tolist()):
    		potential_endpoints.append(endpoint.centroid)
    
    #now we will find the pair of potential endpoints with the maximal separation distance, those are the true endpoints
    if len(potential_endpoints) <= 1:
    	return None 

    pairwise_distances = distance.cdist(potential_endpoints, potential_endpoints, 'euclidean')
    indices_of_max_distance = unravel_index(pairwise_distances.argmax(), pairwise_distances.shape)

    endpoint1 = potential_endpoints[indices_of_max_distance[0]]
    endpoint2 = potential_endpoints[indices_of_max_distance[1]]
    #print endpoint1
    #print endpoint2
    endpoints = [endpoint1, endpoint2]
    return endpoints

def are_joined(endpoint1, endpoint2):
	#given two endpoints calculate the distance between them and return True or False for whether they meet the joining criteria
	cutoff = 5.0 
	distance = distance(endpoint1,endpoint2)
	if distance <= cutoff: 
		return True 

	else:
		return False 


def calc_distance(endpoint1, endpoint2):
	#simple distance calculation
	distance_squared = (endpoint1[0]-endpoint2[0]) * (endpoint1[0]-endpoint2[0]) + (endpoint1[1]-endpoint2[1]) * (endpoint1[1]-endpoint2[1])
	distance = math.sqrt(distance_squared)

	return distance

	

# Line finding using the Probabilistic Hough Transform
tube_lengths = []
tube_angles = []


i=0
#cy3_file_list = os.listdir('6_nt')
root = Tkinter.Tk()
root.withdraw()

file_paths = tkFileDialog.askopenfilenames()
cy3_file_list = list(file_paths)

for i in range(len(cy3_file_list)):
	cy3_file = cy3_file_list[i]

	print "cy3 filename is "+str(cy3_file)
	image_unthresholded = io.imread(cy3_file)

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


	#perfoming edge detection and morphological filling
	edges_open = canny(image, 2, 1, 50) #originally 2,1,25 last param can go up to 500 for improved performance, must lower for poorer images
	#edges_open = canny(image, 2) #originally 2,1,25
	selem = disk(3)#originally 5
	edges = closing(edges_open, selem)
	fill_tubes = ndi.binary_fill_holes(edges)
	io.imsave(cy3_file+"fill_tubes.png", img_as_uint(fill_tubes), cmap=cm.gray)
	cy3_endpoint_mask = make_endpoints_mask(fill_tubes)



	#label image 
	label_image = label(fill_tubes)


	print "detecting nanotube angles...."
	print len(regionprops(label_image))
	for region in regionprops(label_image):
		if region.area/tube_width >= length_cutoff and region.eccentricity >= eccentricity_cutoff:
			region_coords = region.coords.tolist()
			region_endpoints = endpoints(region_coords, cy3_endpoint_mask)
			if region_endpoints == None:
				continue
			endpoint_to_endpoint_vector  = np.subtract(region_endpoints[0], region_endpoints[1])
			x_axis_vector = np.array([0, 1])
			angle_with_x_axis = angle(endpoint_to_endpoint_vector, x_axis_vector)
			angle_with_x_axis *= 180.0/math.pi
			print 'angle with x axis is: ', angle_with_x_axis
			tube_angles.append(angle_with_x_axis)

				



print "printing angles"
f1=open('angles.dat','w+')
for angle in tube_angles:
	print >>f1, angle
f1.close()


