from math import sqrt
from skimage import data
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.color import rgb2gray
from skimage import io
import sys, os
import numpy
import matplotlib.pyplot as plt
plt.switch_backend('agg')

question = str(raw_input('Would you like to use the current directory? [y/n]: '))
if question == 'y':
	pathname = os.getcwd()
else: pathname = str(raw_input('Enter path name: '))

dilution = float(raw_input('Enter the additional dilution factor: '))

#Initiate variables
numS = 0
avgS = 0
n = 0

pathnames = os.listdir(pathname)
output = []

for num in range(len(pathnames)):

	if os.path.isdir(pathnames[num]):
		continue

	image_gray = io.imread(pathnames[num])
	image = io.imread(str(pathnames[num]))
	#Comment out if images saved as color images
	#image_gray = rgb2gray(image)

	#Calculate seeds using threshold of 0.075
	blobs_log = blob_log(image_gray, max_sigma=1, num_sigma=10, threshold=.075, overlap=.1)
	numS = len(blobs_log)

	#Average number of seeds over thresholds 0.07 to 0.08
	'''for i in range(10):
		numS += len(blob_log(image_gray, max_sigma=1, num_sigma=10, threshold=.07 + float(i)/float(1000), overlap=.1))
	numS /= 11'''

	#Add filename and number of seeds to output; increment number of files by 1
	print pathnames[num], numS
	output.append([pathnames[num], numS])
	n += 1

	# Compute radii in the 3rd column.
	blobs_log[:,2] = blobs_log[:,2] * sqrt(2)

	blobs_list = [blobs_log]
	colors = ['yellow']
	titles = ['Laplacian of Gaussian']
	sequence = zip(blobs_list, colors, titles)

	fig, ax = plt.subplots(figsize=(9, 3))
	ax.set_aspect('equal')

	for idx, (blobs, color, title) in enumerate(sequence):
		ax.imshow(image, interpolation='nearest')

	for blob in blobs:
		y, x, r = blob
		c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
		ax.add_patch(c)

	name, ext = os.path.splitext(pathnames[num])
	plt.savefig(name + 'plot.pdf')

column1 = [row[1] for row in output]

#Calculate outliers using the IQR
'''med = numpy.median(column1)
q1 = numpy.percentile(column1, 25)
q3 = numpy.percentile(column1, 75)
iqr = q3 - q1
leftbound = q1 - 1.25*iqr
rightbound = q3 + 1.25*iqr
for i in range(len(column1)):
	if column1[i] < leftbound or column1[i] > rightbound:
		print output[i], "is an outlier"
	else: avgS += column1[i]'''

#Calculate outliers 2 SDs away from the mean
for i in range(len(column1)):
	if abs(column1[i] - numpy.mean(column1)) > 2 * numpy.std(column1):
		print output[i], "is an outlier"
		n -= 1
	else: avgS += column1[i]

avg = str(avgS/float(n))
fov = float(avg) * dilution
print "Average number of seeds: " + avg
print "Seeds per fov for 0.3uL in 19.7uL: " + str(fov)