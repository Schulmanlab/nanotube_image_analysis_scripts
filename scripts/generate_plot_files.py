# -*- coding: utf-8 -*-
from __future__ import unicode_literals 
#script to generate data files needed for plotting
import numpy as np
import joining_ode_class
import matplotlib.pyplot as plt

jode_constant = joining_ode_class.ODE_joining('constant')
jode_hill = joining_ode_class.ODE_joining('hill')
jode_bernie = joining_ode_class.ODE_joining('bernie')
plot_dir_name = "bin10/bootstrapping/"

params = [.4e7, .6e11, 1e8] #optimized parameters
models = [jode_constant, jode_hill, jode_bernie]

########################################################################################################
#plotting simulated distributions with time for each variable of interest here
for i in range(len(models)):
	model = models[i]
	param = params[i]
	line_types = ['red', 'green', 'blue', 'black']
	t, w = model.perform_integration(param)

	t, A_conc_data, B_conc_data, C_conc_data, B_joined_conc_data = model.unpack_concentrations_from_ode_format( w, t)

	#simulated C tube distributions with time
	for timepoint in range(len(t)):
		C_conc_timepoint = C_conc_data[timepoint]
		lengths = [float((bin+1))*model.bin_length for bin in range(len(C_conc_timepoint))]
		line_type = line_types[timepoint]
		plt.plot(lengths, C_conc_timepoint, line_type, label = str(t[timepoint]) + 'seconds')
	plt.legend(loc = 'best')
	plt.xlabel('tube length (µm)')
	plt.ylabel('concentration (M)')
	plt.title('C tube simulated distributions')
	plt.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_C_tubes.pdf')
	plt.close()

	#simulated unjoined B tube distributions with time
	for timepoint in range(len(t)):
		B_conc_timepoint = B_conc_data[timepoint]
		lengths = [float((bin+1))*model.bin_length for bin in range(len(B_conc_timepoint))]
		line_type = line_types[timepoint]
		plt.plot(lengths, B_conc_timepoint, line_type, label = str(t[timepoint]) + 'seconds')
	plt.legend(loc = 'best')
	plt.xlabel('tube length (µm)')
	plt.ylabel('concentration (M)')
	plt.title('unjoined B tube simulated distributions')
	plt.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_unjoined_B_tubes.pdf')
	plt.close()

	#simulated joined B tube distributions with time
	for timepoint in range(len(t)):
		B_joined_conc_timepoint = B_joined_conc_data[timepoint]
		lengths = [float((bin+1))*model.bin_length for bin in range(len(B_joined_conc_timepoint))]
		line_type = line_types[timepoint]
		plt.plot(lengths, B_joined_conc_timepoint, line_type, label = str(t[timepoint]) + 'seconds')
	plt.legend(loc = 'best')
	plt.xlabel('tube length (µm)') 
	plt.ylabel('concentration (M)')
	plt.title('joined B tube simulated distributions')
	plt.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_joined_B_tubes.pdf')
	plt.close()
########################################################################################################

########################################################################################################
#plotting simulated vs experimental length distributions for each variable of interest at each timepoint
#TODO: normalize both distributions for fair comparison
for i in range(len(models)):
	model = models[i]
	param = params[i]
	#line_types = ['red', 'green', 'blue', 'black']
	t, w = model.perform_integration(param)

	#simulated data
	t, sim_A_conc_data, sim_B_conc_data, sim_C_conc_data, sim_B_joined_conc_data = model.unpack_concentrations_from_ode_format( w, t)

	#experimental data
	timepoint_prefixes = ['0hr', '2hr', '4hr']
	expt_B_data = []
	expt_C_data = []
	expt_B_joined_data = []

	for timepoint_prefix in timepoint_prefixes:
		A_unjoined, B_unjoined, A_joined, B_joined, C_joined = model.read_ABC_concentrations(timepoint_prefix)
		expt_B_data.append(B_unjoined)
		expt_C_data.append(C_joined)
		expt_B_joined_data.append(B_joined)

	timepoint_titles = ['.5hr', '2.5hr', '4.5hr']

	for j in range(len(timepoint_titles)):
		expt_B_data_timepoint = [expt_B_data[j][k] / sum(expt_B_data[j]) for k in range(len(expt_B_data[j]))] #normalizing
		#expt_B_data_timepoint = expt_B_data[j]
		expt_C_data_timepoint = [expt_C_data[j][k] / sum(expt_C_data[j]) for k in range(len(expt_C_data[j]))] 
		#expt_C_data_timepoint = expt_C_data[j]
		expt_B_joined_data_timepoint = [expt_B_joined_data[j][k] / sum(expt_B_joined_data[j]) for k in range(len(expt_B_joined_data[j]))]
		#expt_B_joined_data_timepoint = expt_B_joined_data[j]

		sim_B_data_timepoint = [sim_B_conc_data[j+1][k] / sum(sim_B_conc_data[j+1]) for k in range(len(sim_B_conc_data[j+1])) ]#plus one is needed to account for the 0hr timepoint in the simulated data
		sim_C_data_timepoint = [sim_C_conc_data[j+1][k] / sum(sim_C_conc_data[j+1]) for k in range(len(sim_C_conc_data[j+1])) ]
		sim_B_joined_data_timepoint = [sim_B_joined_conc_data[j+1][k] / sum(sim_B_joined_conc_data[j+1]) for k in range(len(sim_B_joined_conc_data[j+1])) ]
		#sim_C_data_timepoint = sim_C_conc_data[j+1]
		#sim_B_joined_data_timepoint = sim_B_joined_conc_data[j+1]

		lengths = [float((bin+1))*model.bin_length for bin in range(len(expt_B_data_timepoint))]

		plt.plot(lengths, sim_B_data_timepoint, 'red', label = 'simulation')
		plt.plot(lengths, expt_B_data_timepoint, 'blue', label = 'experiment')
		plt.legend(loc = 'best')
		plt.xlabel('tube length (µm)') 
		plt.ylabel('fraction')
		plt.title('unjoined B length distribution at '+timepoint_titles[j])
		plt.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_' + timepoint_titles[j]+ '_expt_v_sim_unjoined_B_tubes.pdf')
		plt.close()

		plt.plot(lengths, sim_C_data_timepoint, 'red', label = 'simulation')
		plt.plot(lengths, expt_C_data_timepoint, 'blue', label = 'experiment')
		plt.legend(loc = 'best')
		plt.xlabel('tube length (µm)') 
		plt.ylabel('fraction')
		plt.title('C length distribution at '+timepoint_titles[j])
		plt.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_' + timepoint_titles[j]+ '_expt_v_sim_C_tubes.pdf')
		plt.close()

		plt.plot(lengths, sim_B_joined_data_timepoint, 'red', label = 'simulation')
		plt.plot(lengths, expt_B_joined_data_timepoint, 'blue', label = 'experiment')
		plt.legend(loc = 'best')
		plt.xlabel('tube length (µm)') 
		plt.ylabel('fraction')
		plt.title('joined B length distribution at '+timepoint_titles[j])
		plt.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_' + timepoint_titles[j]+ '_expt_v_sim_joined_B_tubes.pdf')
		plt.close()

########################################################################################################


########################################################################################################
#plot kjoin (without prefactor) vs n for fixed m for bernie and hill 
m = 10.0
lengths = []
kjoins = []
for i in range(jode_hill.n_bins):
	kjoin = jode_hill.kjoin( (i+1)*jode_hill.bin_length, m )
	kjoins.append(kjoin)
	lengths.append((i+1)*jode_hill.bin_length)
plt.plot(lengths, kjoins)
plt.xlabel('tube length (µm)')
plt.ylabel('kjoin (without prefactor)')
plt.title('kjoin (without prefactor) vs n for fixed m = 10 µm')
plt.savefig(plot_dir_name+str(jode_hill.model_name) + '_' + str(model.n_bins) + '_10_kjoin_vs_n.pdf')
plt.close()

m = 10.0
lengths = []
kjoins = []
for i in range(jode_bernie.n_bins):
	kjoin = jode_bernie.kjoin( (i+1)*jode_hill.bin_length, m )
	kjoins.append(kjoin)
	lengths.append((i+1)*jode_bernie.bin_length)
plt.plot(lengths, kjoins)
plt.xlabel('tube length (µm)')
plt.ylabel('kjoin (without prefactor)')
plt.title('kjoin (without prefactor) vs n for fixed m = 10 µm')
plt.savefig(plot_dir_name+str(jode_bernie.model_name) + '_' + str(model.n_bins) + '_10_kjoin_vs_n.pdf')
plt.close()



				

'''B_conc_timepoint = B_conc_data[timepoint]
		B_joined_timepoint = B_joined_conc_data[timepoint]

	lengths = [float((bin+1))*jode_constant.bin_length for bin in range(len(C_conc_timepoint))]
	plt.plot(lengths, C_conc_timepoint)
	plt.savefig('constant_' + str(t[timepoint]) + '_C_tubes.pdf')
	plt.close()

	f1=open('constant_' + str(t[timepoint]) + '_C_tubes.dat','w+')
	for bin in range(len(C_conc_timepoint)):
		print >>f1, float((bin+1))*jode_constant.bin_length, C_conc_timepoint[bin]
	f1.close()

	

	f1=open('constant_' + str(t[timepoint]) + '_B_tubes.dat','w+')
	for bin in range(len(B_conc_timepoint)):
		print >>f1, float((bin+1))*jode_constant.bin_length, B_conc_timepoint[bin]
	f1.close()

	f1=open('constant_' + str(t[timepoint]) + '_B_tubes_joined.dat','w+')
	for bin in range(len(B_joined_timepoint)):
		print >>f1, float((bin+1))*jode_constant.bin_length, B_joined_timepoint[bin]
	f1.close()

#dumping everything for 'hill' joining model....
t, w_constant = jode_hill.perform_integration(2e7)

t, A_conc_data, B_conc_data, C_conc_data, B_joined_conc_data = jode_hill.unpack_concentrations_from_ode_format( w_constant, t)

for timepoint in range(len(t)):
	C_conc_timepoint = C_conc_data[timepoint]
	B_conc_timepoint = B_conc_data[timepoint]
	B_joined_timepoint = B_joined_conc_data[timepoint]


	f1=open('hill_' + str(t[timepoint]) + '_C_tubes.dat','w+')
	for bin in range(len(C_conc_timepoint)):
		print >>f1, float((bin+1))*jode_hill.bin_length, C_conc_timepoint[bin]
	f1.close()

	f1=open('hill_' + str(t[timepoint]) + '_B_tubes.dat','w+')
	for bin in range(len(B_conc_timepoint)):
		print >>f1, float((bin+1))*jode_hill.bin_length, B_conc_timepoint[bin]
	f1.close()

	f1=open('hill_' + str(t[timepoint]) + '_B_tubes_joined.dat','w+')
	for bin in range(len(B_joined_timepoint)):
		print >>f1, float((bin+1))*jode_hill.bin_length, B_joined_timepoint[bin]
	f1.close()

#dumping everything for 'bernie' joining model....
t, w_constant = jode_hill.perform_integration(2e7)

t, A_conc_data, B_conc_data, C_conc_data, B_joined_conc_data = jode_bernie.unpack_concentrations_from_ode_format( w_constant, t)

for timepoint in range(len(t)):
	C_conc_timepoint = C_conc_data[timepoint]
	B_conc_timepoint = B_conc_data[timepoint]
	B_joined_timepoint = B_joined_conc_data[timepoint]


	f1=open('bernie_' + str(t[timepoint]) + '_C_tubes.dat','w+')
	for bin in range(len(C_conc_timepoint)):
		print >>f1, float((bin+1))*jode_bernie.bin_length, C_conc_timepoint[bin]
	f1.close()

	f1=open('bernie_' + str(t[timepoint]) + '_B_tubes.dat','w+')
	for bin in range(len(B_conc_timepoint)):
		print >>f1, float((bin+1))*jode_bernie.bin_length, B_conc_timepoint[bin]
	f1.close()

	f1=open('bernie_' + str(t[timepoint]) + '_B_tubes_joined.dat','w+')
	for bin in range(len(B_joined_timepoint)):
		print >>f1, float((bin+1))*jode_bernie.bin_length, B_joined_timepoint[bin]
	f1.close()'''