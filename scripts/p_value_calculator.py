import numpy as np
import math 
from itertools import combinations
from scipy.stats import ttest_ind

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