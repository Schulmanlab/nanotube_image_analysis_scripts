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

#modifying the joining detection script to measure the angle of Sisi's nanotubes relative to the x-axis of her images 

#constants
#constants
tube_width = 5.0
length_cutoff = 3.0 
eccentricity_cutoff = 0.5
end_to_end_distance_cutoff = 10.0

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
cy3_file_list = os.listdir('6_nt')

for i in range(len(cy3_file_list)):
	cy3_file = cy3_file_list[i]

	print "cy3 filename is "+str(cy3_file)
	image_unthresholded = io.imread('cy3/'+cy3_file)

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
	edges_open = canny(image, 3, 1, 25) #originally 2,1,25
	#edges_open = canny(image, 2) #originally 2,1,25
	selem = disk(5)
	edges = closing(edges_open, selem)
	fill_tubes = ndi.binary_fill_holes(edges)
	io.imsave(cy3_file+"fill_tubes.png", img_as_uint(fill_tubes), cmap=cm.gray)
	cy3_endpoint_mask = make_endpoints_mask(fill_tubes)


'''
	#label image 
	label_image = label(fill_tubes)
	label_image_647 = label(fill_tubes_647)

	regions_joined_647 = []
	regions_joined_cy3 = []
	print "detecting joining"
	print len(regionprops(label_image_647))
	for region_647 in regionprops(label_image_647):
		is_joined = 0
		if region_647.area/tube_width >= length_cutoff and region_647.eccentricity >= eccentricity_cutoff:
			for region in regionprops(label_image):
				if region.area/tube_width < length_cutoff or region.eccentricity < eccentricity_cutoff:
					continue
				region_647_coords = region_647.coords.tolist()
				region_coords = region.coords.tolist()

				region_647_endpoints = endpoints(region_647_coords, atto647_endpoint_mask)
				region_endpoints = endpoints(region_coords, cy3_endpoint_mask)

				if region_647_endpoints == None or region_endpoints == None:
					continue 
				#print region_647_endpoints
				#print region_endpoints

				#now calculate all pairwsie distances between the two sets of endpoints
				pairwise_distances = distance.cdist(region_647_endpoints, region_endpoints, 'euclidean')
				minimum_distance = pairwise_distances.min()
				
				if minimum_distance < end_to_end_distance_cutoff:
					print minimum_distance
					is_joined = 1
					#print "we have joining!"
					break
		if is_joined == 1:
			lengths_647_joined.append(region_647.area/tube_width)
			regions_joined_647.append(region_647)
			lengths_cy3_joined.append(region.area/tube_width)
			regions_joined_cy3.append(region)


	print "printing 647 components that are not joined"
	for region_647 in regionprops(label_image_647):
		if region_647 not in regions_joined_647 and region_647.area/tube_width>= length_cutoff and region_647.eccentricity>= eccentricity_cutoff:
			lengths_647_unjoined.append(region_647.area/tube_width)

	print "printing cy3 components that are not joined"
	for region in regionprops(label_image):
		if region not in regions_joined_cy3 and region.area/tube_width>= length_cutoff and region.eccentricity>= eccentricity_cutoff:
			lengths_cy3_unjoined.append(region.area/tube_width)

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

	'''