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


#xfmt = matplotlib.ticker.ScalarFormatter()
#xfmt.set_powerlimits((-3,3))

f_cost_fxns, axarr_cost_fxns = plt.subplots(3,3, sharex = False, sharey = False, figsize = (9,9))


#here starts the main script to brute force try many constant joining values and evaluate the error with expt
#*abdul sees 3e6 per molar per second in his work*
#the kjoins predicted by bernie/hill are ~20 and ~2, so fairly close to 1, makes sense that roughly the same range of param would work

#this script will perform the same brute force search as the original brute_force_parameter_search but it will use three different cost functions
#and will output a 3x3 grid of subplots for inclusion in our paper

jode_constant = joining_ode_class.ODE_joining('constant')
jode_hill = joining_ode_class.ODE_joining('hill')
jode_bernie = joining_ode_class.ODE_joining('bernie')
plot_dir_name = "bin10/optimal_parameter/"
n_bootstrap = 1 
n_bins = 10
max_tube_length = 10.0

models = ['constant','hill','bernie']
cost_fxns = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

for model in models:
	for cost_fxn in cost_fxns:
		jode = joining_ode_class.ODE_joining(joining_model = model, weight_joining_percentage = cost_fxn[0], weight_a_cdf = cost_fxn[1], weight_c_cdf = cost_fxn[2])
		if model == 'hill':
			param_list = np.linspace(1e10, 1.5e11, 100)
		if model == 'bernie':
			param_list = np.linspace(1e7, 1.5e8, 100)
		if model == 'constant':
			param_list = np.linspace(1e6, .8e7, 100)

		kjoin_list = [param_list[i] for i in range(len(param_list))]
		#original: kjoin_list = [1e10+(i*1e9) for i in range(200)]
		squared_error_list = []
		#kjoin_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14]
		#kjoin_list = [10000*(i*10) for i in range(10)]

		
		for kjoin in kjoin_list:
			squared_error_list.append( jode.vahid_error(float(kjoin)) )


		axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].plot(kjoin_list, squared_error_list)
		if model == 'constant':
			if cost_fxn[2] == 1.0:
				axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].set_ylim([0.3, 0.5])
			else:
				axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].set_ylim([0.0, 0.2])

		if model == 'hill':
			if cost_fxn[0] == 1.0:
				axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].set_ylim([0.0, 0.2])
			if cost_fxn[1] == 1.0:
				axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].set_ylim([0.3, 0.5])
			if cost_fxn[2] == 1.0:
				axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].set_ylim([1.4, 1.6])

		if model == 'bernie':
			if cost_fxn[0] == 1.0:
				axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].set_ylim([0.0, 0.2])
			if cost_fxn[1] == 1.0:
				axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].set_ylim([0.0, 0.2])
			if cost_fxn[2] == 1.0:
				axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].set_ylim([0.7, 0.9])

		axarr_cost_fxns[models.index(model),cost_fxns.index(cost_fxn)].get_xaxis().get_major_formatter().set_powerlimits((0, 0))


f_cost_fxns.set_tight_layout(True)
f_cost_fxns.savefig(plot_dir_name + 'cost_fxn_subplots.pdf')







'''
param_list = np.linspace(1e7, 3e8, 400)
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







#optimal_kjoin_table = []
param_list = np.linspace(1e6, .8e7, 400)
kjoin_list = [param_list[i] for i in range(len(param_list))]
#kjoin_list = [1e6+(i*1e5) for i in range(100)]
squared_error_list = []
#kjoin_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14]
#kjoin_list = [10000*(i*10) for i in range(10)]

for kjoin in kjoin_list:
	squared_error_list.append( jode_constant.vahid_error(float(kjoin)) )

best_kjoin_full_data_constant = kjoin_list[squared_error_list.index(min(squared_error_list))]
#error_specific = squared_error_list [kjoin_list.index(5.8e6)]
#print 'error for 5.8e6: '+str(error_specific)
best_error_full_data_constant = min(squared_error_list)
#optimal_kjoin_table.append(["constant", best_kjoin, best_error])

plt.plot(kjoin_list, squared_error_list)
plt.xlabel('kjoin prefactor')
plt.ylabel('squared error')
plt.title('Constant kjoin model brute force optimization')
plt.ylim([.0,.2])
plt.savefig(plot_dir_name + 'Constant_kjoin_model_brute_force_optimization.pdf')
plt.close()'''

