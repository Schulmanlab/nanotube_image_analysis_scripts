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
from skimage.external import tifffile


#modifying the joining detection script to measure the angle of Sisi's nanotubes relative to the x-axis of her images 
#modifying this script further to measure a time series of angles for individual tubes and performa an autocorrelation 
#analysis to determine the relaxation time 
#need a method to assign nanotubes a label based on their endpoint positions, referencing against a record of nanotubes
#detected from all previous images in the time series 


#constants
tube_width = 7.0
length_cutoff = 7.0 
eccentricity_cutoff = 0.5


def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
	#going to store a time series of endpoint to endpoint vectors instead of angles 
	#this will make calculating the MSAD for a given delta t easier
	#we will not need to worry about finding the minimum angle between two time points (just use the dot product)
	print v1
	print v2

	return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

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

def find_seed_attachment_point(endpoints):
	#hack, going to hard code the approximate attachment point and find the endpoint that is closer to that one
	#this could be done more elegantly by going through the image stack and finding the point that does not change
	approx_seed_attachment_point = np.array([33,33])

	endpoint1_distance = calc_distance(endpoints[0], approx_seed_attachment_point)
	endpoint2_distance = calc_distance(endpoints[1], approx_seed_attachment_point)

	if endpoint1_distance < endpoint2_distance:
		return endpoints[0], endpoints[1]
	else:
		return endpoints[1], endpoints[0]


def process_images():
	tube_lengths = []
	tube_angles = []


	i=1

	cy3_image_stack = tifffile.imread("0_P1_tube_2.tif")


	for image in cy3_image_stack:
		total_images = len(cy3_image_stack)
		current_frame = i
		print "processing frame " +str(i) + " of "+str(total_images)


		#perfoming edge detection and morphological filling
		edges_open = canny(image, 2, 1, 50) #originally 2,1,25
		selem = disk(3)#originally 5
		edges = closing(edges_open, selem)
		fill_tubes = ndi.binary_fill_holes(edges)
		io.imsave(str(i)+"_fill_tubes.png", img_as_uint(fill_tubes), cmap=cm.gray)
		cy3_endpoint_mask = make_endpoints_mask(fill_tubes)


		#label image 
		label_image = label(fill_tubes)


		print "detecting nanotube angles...."
		print len(regionprops(label_image))

		if len(regionprops(label_image)) == 0:
			#image contains no nanotubes
			print 'angle cannot be determined: assigning uniform distribution'
			tube_angles.append('uniform')


		#we need to determine if a nanotube is present and if it is elongated enough to make a determination of 
		#the angle in the x-y plane

		for region in regionprops(label_image):
			if region.area/tube_width >= length_cutoff and region.eccentricity >= eccentricity_cutoff:
				region_coords = region.coords.tolist()
				region_endpoints = endpoints(region_coords, cy3_endpoint_mask)
				if region_endpoints == None:
					continue

				#before calculating the endpoint to endpoint vector we need to make sure that the root point is the same
				#the second endpoint in the subtraction should always be the seed attachment point 
				#seed_attachment_point = find_seed_attachment_point(region_endpoints, previous_region_endpoints)
				seed_attachment_point, diffusing_end = find_seed_attachment_point(region_endpoints)
				endpoint_to_endpoint_vector  = np.subtract(diffusing_end, seed_attachment_point)
				print endpoint_to_endpoint_vector
				x_axis_vector = np.array([0.0, 1.0])
				angle_with_x_axis = angle(endpoint_to_endpoint_vector, x_axis_vector)
				angle_with_x_axis *= 180.0/math.pi
				print 'angle with x axis is: ', angle_with_x_axis
				tube_angles.append(angle_with_x_axis)

				#previous_region_endpoints = region_endpoints
		i+=1

	return tube_angles

				



print "printing angles"
f1=open('angles.dat','w+')
tube_angles = process_images()
for angle in tube_angles:
	print >>f1, angle
f1.close()


