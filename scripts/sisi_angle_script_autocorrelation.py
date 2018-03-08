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
import scipy.stats
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
#from pymc3 import *
import pymc
from StringIO import StringIO

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

	#v1 = np.array([int(v1[0]), int(v1[1])])
	#v2 = np.array([int(v2[0]), int(v2[1])])
	#return sum((a*b) for a, b in zip(v1, v2))
	return np.dot(v1, v2)

def length(v):
	#print "dot product result is: " + str(math.sqrt(dotproduct(v, v)))
	#return math.sqrt(dotproduct(v, v))
	return np.linalg.norm(v)

def angle(v1, v2):
	#going to store a time series of endpoint to endpoint vectors instead of angles 
	#this will make calculating the MSAD for a given delta t easier
	#we will not need to worry about finding the minimum angle between two time points (just use the dot product)
	#returns angle is radians not degrees 
	print 'dot product is: ', dotproduct(v1,v2)
	print v1
	print v2
	print length(v1)
	print length(v2)
	print dotproduct(v1, v2) / (length(v1) * length(v2))
	if dotproduct(v1, v2) / (length(v1) * length(v2)) >= .9999 and dotproduct(v1, v2) / (length(v1) * length(v2)) <= 1.0001:
		print "hit condition"
		return 0.0
	else:
		return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def signed_angle_with_x_axis(vector):
	#need a signed angle (to distinguish +45 from -45 from the x axis) this is purely needed for calculating correlation coefficients
	#will need to compare the y coordinate of the diffusing end to the y coordinate of the seed attachment point
	#to determine whether the angle is positive or negative 

	x_axis_vector = np.array([0.0, 1.0]) #looks like y but remember matrix indices are row, column

	#find unsigned angle
	unsigned_angle = math.acos(dotproduct(vector, x_axis_vector) / (length(vector) * length(x_axis_vector)))

	#find the sign using the y component of the vector 
	if vector[0] >= 0:
		return unsigned_angle

	else:
		return -1.0 * unsigned_angle 

def vector_from_x_axis_angle(angle):
	#note this will only work if the mean angle from the x axis is very close to zero... I think I fixed this  
	#convert angle to radians
	#print 'angle is ', angle 
	angle *= math.pi/180.0
	tan_angle = math.tan(angle)
	#tan = opposite over adjacent = y/x = vector[0]/vector[1]
	#if angle is negative the y component should be negative
	if angle <= 0.0:
		vector0 = -1.0
		vector1 = 1.0/tan_angle

	else:
		vector0 = 1.0
		vector1 = 1.0/tan_angle

	vector = np.array([vector0, vector1])

	return vector 

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

def calc_angle_variance(tube_vectors, mean_tube_vector):
	#calculate the variance in the distribution of angles
	#this needs to be done by using the dot product for each point 
	squared_angle_differences = []
	for tube_vector in tube_vectors:
		angle_difference = angle(tube_vector, mean_tube_vector)
		squared_angle_differences.append(angle_difference * angle_difference)

	n, min_max, mean, var, skew, kurt = scipy.stats.describe( squared_angle_differences )

	return mean 



def find_seed_attachment_point(endpoints):
	#hack, going to hard code the approximate attachment point and find the endpoint that is closer to that one
	#this could be done more elegantly by going through the image stack and finding the point that does not change
	approx_seed_attachment_point = np.array([31,23])

	endpoint1_distance = calc_distance(endpoints[0], approx_seed_attachment_point)
	endpoint2_distance = calc_distance(endpoints[1], approx_seed_attachment_point)

	if endpoint1_distance < endpoint2_distance:
		return endpoints[0], endpoints[1]
	else:
		return endpoints[1], endpoints[0]

def msad(tube_vectors, delta_t):
	#for a given delta_t find the distribution of mean square angular deviations
	#images are spaced 2 seconds apart
	#for 'uniform' entries we will select an angle from a uniform distribution

	#first let's go through the list and replace 'uniform' entries with a random vector
	final_tube_vectors = []
	for tube_vector in tube_vectors:
		print tube_vector
		if tube_vector == 'uniform':
			new_vector = np.array([np.random.random()*1.0 - 2.0, np.random.random()*1.0 - 2.0])
		else: 
			new_vector = tube_vector
		final_tube_vectors.append(new_vector)

	print final_tube_vectors
	angular_displacements = []
	dot_products = []
	for i in range(len(tube_vectors)-delta_t):
		angular_displacement = angle(final_tube_vectors[i], final_tube_vectors[i+delta_t])
		angular_displacements.append(angular_displacement * angular_displacement)
		dot_product = np.dot(final_tube_vectors[i]/np.linalg.norm(final_tube_vectors[i]), final_tube_vectors[i+delta_t]/np.linalg.norm(final_tube_vectors[i+delta_t]) )
		dot_products.append(dot_product)
	print angular_displacements
	n, min_max, mean, var, skew, kurt = scipy.stats.describe( angular_displacements )
	mean_angular_displacement = mean
	n, min_max, mean, var, skew, kurt = scipy.stats.describe( dot_products )
	mean_dot_product = mean
	return mean_angular_displacement, mean_dot_product

def mean_msad(tube_vectors, delta_t):
	#just repeating the msad calculation n times and taking the average 
	msad_list = []
	mean_dot_product_list = []
	if delta_t == 0.0:
		return 0.0, 1.0, 0.0
	for i in range(10):
		msad_, dot_product = msad(tube_vectors, delta_t)
		msad_list.append(msad_)
		mean_dot_product_list.append(dot_product)
	n, min_max, mean, var, skew, kurt = scipy.stats.describe( msad_list)
	mean_msad = mean
	msad_variance = var
	n, min_max, mean, var, skew, kurt = scipy.stats.describe( mean_dot_product_list)
	mean_dot_product = mean

	return mean_msad, mean_dot_product, msad_variance

def pearson_r(tube_vectors, tube_angles, delta_t):
	#calculating pearson correlation coefficient for a given delta t
	#this will have to be repeated many times to account for uncertainty in some angular measurments 
	#need to implement positive/negative angles with a given axis here 

	#first we need to calculate the mean angle with the x-axis for the signal and the shifted signal 
	#correction: I am going to assume that the mean and variance do not change with time and are the same for the 
	#original and shifted signals
	final_tube_angles = []
	final_tube_vectors = []
	for i in range(len(tube_angles)):
		print tube_angles[i]
		if tube_angles[i] == 'uniform':
			new_angle = np.random.random()*360.0 - 180.0 
			new_vector = vector_from_x_axis_angle(new_angle)
		else: 
			new_angle = tube_angles[i]
			new_vector = tube_vectors[i]
		final_tube_angles.append(new_angle)
		final_tube_vectors.append(new_vector)
	print final_tube_angles 

	n, min_max, mean, var, skew, kurt = scipy.stats.describe(final_tube_angles)
	mean_angle_with_x_axis = mean

	#now we need to convert the mean angle into a vector
	mean_vector = vector_from_x_axis_angle(mean_angle_with_x_axis)
	#this is not a valid way to calc variance here b/c we are using the dot product to find the difference 
	#between each point and the mean 

	var_angle_with_x_axis = calc_angle_variance(final_tube_vectors, mean_vector)
	print "mean angle with x_axis is: ", mean_angle_with_x_axis
	print "mean vector is: ", mean_vector

	#now we calculate < x_i - x_mean > using the dot product to find the difference between the given angle and the mean
	#we do this for both the original and shifted signals 

	correlation_numerator = []
	for i in range(len(tube_vectors)-delta_t):
		angle_mean_difference = angle(final_tube_vectors[i], mean_vector)
		#print "mean_vector is: ", mean_vector
		#print "current vector is: ", final_tube_vectors[i]
		print "angle difference is: ", angle_mean_difference
		shifted_angle_mean_difference = angle(final_tube_vectors[i+delta_t], mean_vector)
		#print "shifted vector is: ", final_tube_vectors[i+delta_t]
		print "shifted angle difference is: ", shifted_angle_mean_difference
		correlation_numerator.append(angle_mean_difference * shifted_angle_mean_difference)

	n, min_max, mean, var, skew, kurt = scipy.stats.describe( correlation_numerator )
	correlation_numerator_mean = mean
	correlation_denominator = var_angle_with_x_axis

	pearson_correlation = correlation_numerator_mean/correlation_denominator


	return pearson_correlation 











'''def plot_mean_msad_vs_delta_t():
	#generate plot

def estimate_diffusion_coefficient():
	#estimate diffusion coefficient'''



def process_images():
	tube_lengths = []
	tube_angles = []
	tube_vectors = []


	i=1

	cy3_image_stack = tifffile.imread("0_2um_tube.tif")


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
			tube_vectors.append('uniform')
			continue


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
				signed_x_axis_angle = signed_angle_with_x_axis(endpoint_to_endpoint_vector)
				signed_x_axis_angle *= 180.0/math.pi
				print 'angle with x axis is: ', angle_with_x_axis
				print 'signed angle with x axis is: ', signed_x_axis_angle
				tube_angles.append(signed_x_axis_angle)
				tube_vectors.append(endpoint_to_endpoint_vector)

				#previous_region_endpoints = region_endpoints

			else:
				tube_angles.append('uniform')
				tube_vectors.append('uniform')
				break


		i+=1
	print tube_vectors
	print len(tube_vectors)
	return tube_angles, tube_vectors

				

tube_angles, tube_vectors = process_images()
mean_angle_displacement, mean_dot_product = msad(tube_vectors, 0)
print "mean angle displacement: ", mean_angle_displacement
print "mean dot product: ", mean_dot_product
random_test_tube_angles = []
random_test_tube_vectors = []
#for i in range(100000):
#	random_test_tube_angles.append('uniform')
#	random_test_tube_vectors.append('uniform')

#pearson_correlation = pearson_r(random_test_tube_vectors, random_test_tube_angles, 10)
#pearson_correlation = pearson_r(tube_vectors, tube_angles, 10)
#print "pearson correlation is: ", pearson_correlation
msad_list = []
delta_t_list = range(0,10)
delta_t_plot_list = []
dot_product_list = []
msad_variance_list = []
for delta_t in delta_t_list:
	msad_, dot_product, msad_variance = mean_msad(tube_vectors,delta_t)
	msad_list.append(msad_)
	dot_product_list.append(dot_product)
	delta_t_plot_list.append(delta_t*2.0)
	msad_variance_list.append(msad_variance)

plt.plot(delta_t_plot_list, msad_list)
plt.xlabel('n * delta_t (seconds) ')
plt.ylabel('MSAD (radians^2)')
#plt.title('Hill model brute force optimization')
plt.savefig('msad_plot.pdf')
plt.close()

print dot_product_list
print "msad list is: ",msad_list
print "msad variance list is: ",msad_variance_list
print "delta_t list is: ",delta_t_plot_list
plt.plot(delta_t_plot_list, dot_product_list)
plt.xlabel('n * delta_t (seconds) ')
plt.ylabel('< dot product between vec(t), vec(t+delta_t) >')
#plt.title('Hill model brute force optimization')
plt.savefig('autocorrelation_plot.pdf')
plt.close()

x_data = delta_t_plot_list
y_data = msad_list
y_variance = sum(msad_variance_list)/float(len(msad_variance_list))
print "mean msad_variance is: ", y_variance

data = dict(x = x_data, y = y_data)

'''with Model() as model: # model specifications in PyMC3 are wrapped in a with-statement
    # Define priors
    #sigma = HalfCauchy('sigma', beta=10, testval=1.)
    sigma = y_variance 
    intercept = 0
    x_coeff = Normal('x', .25, sd=.3)

    # Define likelihood
    likelihood = Normal('y', mu=intercept + x_coeff * x_data,
                        sd=sigma, observed=y_data)

    # Inference!
    trace = sample(3000, cores=2) # draw 3000 posterior samples using
'''



alpha = pymc.Uniform('alpha', lower=0.1, upper=1)

x = pymc.Normal('x', mu=0, tau=1, value=x_data, observed=True)

@pymc.deterministic(plot=False)
def linear_regress(x=x, alpha=alpha):
    return x*alpha

y = pymc.Normal('output', mu=linear_regress, tau = 1.0/y_variance, value=y_data, observed=True)

model = pymc.Model([x, y, alpha])
mcmc = pymc.MCMC(model)
mcmc.sample(iter=100000, burn=10000, thin=10)

alpha.summary()

result = StringIO()
 
sys.stdout = result

alpha.summary()

result_string = result.getvalue()

print "printing MCMC summary"
f1=open('bayesian_regression.dat','w+')
#tube_angles = process_images()
print >>f1, result_string
f1.close()






'''print "printing angles"
f1=open('angles.dat','w+')
#tube_angles = process_images()
for angle in tube_angles:
	print >>f1, angle
f1.close()
'''

