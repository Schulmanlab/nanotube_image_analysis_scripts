# -*- coding: utf-8 -*-
from __future__ import unicode_literals 
#script to generate data files needed for plotting
import numpy as np
import joining_ode_class
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from tabulate import tabulate
import random
import scipy as sp
import scipy.stats
import math
from mpi4py import MPI
comm = MPI.COMM_WORLD

def confidence_interval(data, confidence=0.90):
	s = np.array(data)
	n, min_max, mean, var, skew, kurt = scipy.stats.describe( s )
	normalized_var = var/math.sqrt(0.5)
	std = math.sqrt(normalized_var)
	interval = scipy.stats.norm.interval(0.90, loc = mean, scale = std )
	return interval 


#here starts the main script to brute force try many constant joining values and evaluate the error with expt
#*abdul sees 3e6 per molar per second in his work*
#the kjoins predicted by bernie/hill are ~20 and ~2, so fairly close to 1, makes sense that roughly the same range of param would work

jode_constant = joining_ode_class.ODE_joining('constant')
jode_hill = joining_ode_class.ODE_joining('hill')
jode_bernie = joining_ode_class.ODE_joining('bernie')
plot_dir_name = "bin10/discontinuity_testing/"
n_bootstrap = 2 
n_bins = 10
max_tube_length = 10.0

def optimal_parameter_hill(hill_model):
	param_list = np.linspace(1e10, 1.5e11, 400)
	kjoin_list = [param_list[i] for i in range(len(param_list))]
	#original: kjoin_list = [1e10+(i*1e9) for i in range(200)]
	squared_error_list = []
	#kjoin_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14]
	#kjoin_list = [10000*(i*10) for i in range(10)]

	optimal_kjoin_table = []
	for kjoin in kjoin_list:
		squared_error_list.append( hill_model.vahid_error(float(kjoin)) )

	best_kjoin = kjoin_list[squared_error_list.index(min(squared_error_list))]
	best_error = min(squared_error_list)
	#optimal_kjoin_table.append(["Hill", best_kjoin, best_error])
	return best_kjoin, best_error

def optimal_parameter_bernie(bernie_model):
	#original: kjoin_list = [1e7+(i*1e6) for i in range(200)]
	param_list = np.linspace(1e7, 1.5e8, 400)
	kjoin_list = [param_list[i] for i in range(len(param_list))]
	squared_error_list = []
	#kjoin_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14]
	#kjoin_list = [10000*(i*10) for i in range(10)]
	for kjoin in kjoin_list:
		squared_error_list.append( bernie_model.vahid_error(float(kjoin)) )

	best_kjoin = kjoin_list[squared_error_list.index(min(squared_error_list))]
	best_error = min(squared_error_list)

	return best_kjoin, best_error

def optimal_parameter_constant(constant_model):
	#original: kjoin_list = [1e6+(i*1e5) for i in range(100)]
	param_list = np.linspace(1e6, .8e7, 400)
	kjoin_list = [param_list[i] for i in range(len(param_list))]
	squared_error_list = []
	#kjoin_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14]
	#kjoin_list = [10000*(i*10) for i in range(10)]
	for kjoin in kjoin_list:
		squared_error_list.append( constant_model.vahid_error(float(kjoin)) )

	best_kjoin = kjoin_list[squared_error_list.index(min(squared_error_list))]
	best_error = min(squared_error_list)

	return best_kjoin, best_error


'''param_list = np.linspace(1e10, 1.5e11, 400)
kjoin_list = [param_list[i] for i in range(len(param_list))]
#original: kjoin_list = [1e10+(i*1e9) for i in range(200)]
squared_error_list = []
#kjoin_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14]
#kjoin_list = [10000*(i*10) for i in range(10)]

optimal_kjoin_table = []
for kjoin in kjoin_list:
	squared_error_list.append( jode_hill.vahid_error(float(kjoin)) )'''

if comm.rank == 0:
	best_kjoin_full_data_hill, best_error_full_data_hill = optimal_parameter_hill(jode_hill)

best_kjoin_list = []
error_boot_list = []
for i in range(n_bootstrap):
	hill_model = joining_ode_class.ODE_joining(joining_model = 'hill', n_bins = n_bins, max_tube_length = max_tube_length, bootstrap = True)
	best_kjoin, best_error = optimal_parameter_hill(hill_model)
	error_boot = hill_model.vahid_error(best_kjoin_full_data_hill)
	error_boot_list.append(error_boot)
	#print 'best parameter: '+str(best_kjoin)
	#print 'best error: '+str(best_error)
	#print best_kjoin, best_error, 
	#if best_kjoin >= 1e11:
	best_kjoin_list.append(best_kjoin)

#here we will send out the best_kjoin_list for this process to the head process (is it ok for a process to send to itself? we will find out)
if comm.rank != 0:
	comm.send(best_kjoin_list, dest=0, tag = 11)
	comm.send(error_boot_list, dest=0, tag = 12)
#here is where we want to append the other best_kjoin_lists from the all processes running
if comm.rank == 0: 
	for i in range(comm.size - 1):
		slave_kjoin_list = comm.recv(source = i+1 , tag = 11)
		slave_error_list = comm.recv(source = i+1 , tag = 12)
		for slave_kjoin in slave_kjoin_list:
			best_kjoin_list.append( slave_kjoin)
		for slave_error in slave_error_list:
			error_boot_list.append( slave_error)

	interval_parameter = confidence_interval(best_kjoin_list)
	interval_error = confidence_interval(error_boot_list)

	#print best_kjoin_list
	print "best Hill parameter: ", best_kjoin_full_data_hill
	print "best Hill error: ", best_error_full_data_hill
	print ".90 CI Hill parameter: ", interval_parameter
	print ".90 CI Hill error: ", interval_error


#optimal_kjoin_table.append(["Hill", best_kjoin_full_data_hill, best_error_full_data_hill, interval_parameter, interval_error])









param_list = np.linspace(1e7, 1.5e8, 400)
kjoin_list = [param_list[i] for i in range(len(param_list))]
#original: kjoin_list = [1e7+(i*1e6) for i in range(200)]
squared_error_list = []
#kjoin_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14]
#kjoin_list = [10000*(i*10) for i in range(10)]
for kjoin in kjoin_list:
	squared_error_list.append( jode_bernie.vahid_error(float(kjoin)) )

best_kjoin_full_data_bernie = kjoin_list[squared_error_list.index(min(squared_error_list))]
best_error_full_data_bernie = min(squared_error_list)
#optimal_kjoin_table.append(["Bernie", best_kjoin, best_error])

plt.plot(kjoin_list, squared_error_list)
plt.xlabel('kjoin prefactor')
plt.ylabel('squared error')
plt.title('Bernie model brute force optimization')
plt.savefig(plot_dir_name + 'Bernie_model_brute_force_optimization.pdf')
plt.close()


best_kjoin_list = []
error_boot_list = []
for i in range(n_bootstrap):
	bernie_model = joining_ode_class.ODE_joining(joining_model = 'bernie', n_bins = n_bins, max_tube_length = max_tube_length, bootstrap = True)
	best_kjoin, best_error = optimal_parameter_bernie(bernie_model)
	error_boot = bernie_model.vahid_error(best_kjoin_full_data_bernie)
	error_boot_list.append(error_boot)
	#print best_kjoin, best_error, bernie_model.random_seed
	
	best_kjoin_list.append(best_kjoin)

interval_parameter = confidence_interval(best_kjoin_list)
interval_error = confidence_interval(error_boot_list)


optimal_kjoin_table.append(["Bernie", best_kjoin_full_data_bernie, best_error_full_data_bernie, interval_parameter, interval_error])




#optimal_kjoin_table = []
param_list = np.linspace(1e6, .8e7, 400)
kjoin_list = [param_list[i] for i in range(len(param_list))]
kjoin_list = [1e6+(i*1e5) for i in range(100)]
squared_error_list = []
#kjoin_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14]
#kjoin_list = [10000*(i*10) for i in range(10)]

for kjoin in kjoin_list:
	squared_error_list.append( jode_constant.vahid_error(float(kjoin)) )

best_kjoin_full_data_constant = kjoin_list[squared_error_list.index(min(squared_error_list))]
error_specific = squared_error_list [kjoin_list.index(5.8e6)]
#print 'error for 5.8e6: '+str(error_specific)
best_error_full_data_constant = min(squared_error_list)
#optimal_kjoin_table.append(["constant", best_kjoin, best_error])

plt.plot(kjoin_list, squared_error_list)
plt.xlabel('kjoin prefactor')
plt.ylabel('squared error')
plt.title('Constant kjoin model brute force optimization')
#plt.ylim([.35,.45])
plt.savefig(plot_dir_name + 'Constant_kjoin_model_brute_force_optimization.pdf')
plt.close()


best_kjoin_list = []
error_boot_list = []
for i in range(n_bootstrap):
	constant_model = joining_ode_class.ODE_joining(joining_model = 'constant', n_bins = n_bins, max_tube_length = max_tube_length, bootstrap = True)
	best_kjoin, best_error = optimal_parameter_constant(constant_model)
	error_boot = constant_model.vahid_error(best_kjoin_full_data_constant)
	error_boot_list.append(error_boot)
	#print best_kjoin, best_error, constant_model.random_seed
	#print 'best parameter: '+str(best_kjoin)
	#print 'best error: '+str(best_error)
	#print 'error associated with best full data parameter: '+str(error_boot)
	
	best_kjoin_list.append(best_kjoin)

interval_parameter = confidence_interval(best_kjoin_list)
interval_error = confidence_interval(error_boot_list)


optimal_kjoin_table.append(["constant", best_kjoin_full_data_constant, best_error_full_data_constant, interval_parameter, interval_error])

#print tabulate(optimal_kjoin_table, headers = ["model", "fitted parameter", "squared error", ".90 parameter CI", ".90 error CI"])




