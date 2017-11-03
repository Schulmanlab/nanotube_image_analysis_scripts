import numpy as np
import math 
from itertools import combinations
from scipy.stats import ttest_ind
import random

with open("hill_errors_boot_list.dat") as f:
	content = f.readlines()
hill_errors = [float(x.strip()) for x in content]

with open("bernie_errors_boot_list.dat") as f:
	content = f.readlines()
bernie_errors = [float(x.strip()) for x in content]

with open("constant_errors_boot_list.dat") as f:
	content = f.readlines()
constant_errors = [float(x.strip()) for x in content]

data = {
	'hill' : hill_errors,
	'bernie' : bernie_errors,
	'constant' : constant_errors,
}
for list1, list2 in combinations(data.keys(), 2):
	t, p = ttest_ind(data[list1], data[list2])
	print list1, list2, p

success = 0 
for i in range (10000):
	hill_error = random.choice(hill_errors)
	bernie_error = random.choice(bernie_errors)
	constant_error = random.choice(constant_errors)

	if constant_error < bernie_error < hill_error:
		success += 1

prob = float(success)/10000.0
print "probability of constant < bernie < hill is: "