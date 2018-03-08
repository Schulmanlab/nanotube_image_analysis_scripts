# -*- coding: utf-8 -*-
from __future__ import unicode_literals 
#script to generate data files needed for plotting
import numpy as np
import joining_ode_class
import matplotlib.pyplot as plt
import scipy.stats
from scipy.stats import gaussian_kde
from scipy.interpolate import UnivariateSpline

jode_constant = joining_ode_class.ODE_joining('constant')
jode_hill = joining_ode_class.ODE_joining('hill')
jode_bernie = joining_ode_class.ODE_joining('bernie')
plot_dir_name = "bin10/optimal_parameter/"

params = [3.86e6, 2.98e10, 7.175e7] #optimized parameters
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
plt.savefig(plot_dir_name+"scaling_range_"+str(jode_hill.model_name) + '_' + str(jode_hill.n_bins) + '_10_kjoin_vs_n.pdf')
plt.close()

m = 10.0
lengths = []
kjoins = []
for i in range(jode_bernie.n_bins):
	kjoin = jode_bernie.kjoin( (i+1)*jode_bernie.bin_length, m )
	kjoins.append(kjoin)
	lengths.append((i+1)*jode_bernie.bin_length)
plt.plot(lengths, kjoins)
plt.xlabel('tube length (µm)')
plt.ylabel('kjoin (without prefactor)')
plt.title('kjoin (without prefactor) vs n for fixed m = 10 µm')
plt.savefig(plot_dir_name+"scaling_range_"+str(jode_bernie.model_name) + '_' + str(jode_bernie.n_bins) + '_10_kjoin_vs_n.pdf')
plt.close()

########################################################################################################
#plotting mean B tube length vs time with standard error of the mean calculated from bootstrapping
expt_B_unjoined_data = []
expt_B_joined_data = []
timepoint_prefixes = ['0hr', '2hr', '4hr']
for timepoint_prefix in timepoint_prefixes:
		A_unjoined_raw_data, B_unjoined_raw_data, C_joined_raw_data, A_joined_raw_data, B_joined_raw_data = jode_constant.read_ABC_lengths(timepoint_prefix)
		expt_B_unjoined_data.append(B_unjoined_raw_data)
		expt_B_joined_data.append(B_joined_raw_data)

timepoints = [.5, 2.5, 4.5]
length_means = []
length_standard_error_of_mean = []
for i in range(len(timepoints)):
	#first calculate mean 
	total_b_tubes = expt_B_joined_data[i] + expt_B_unjoined_data[i]
	length_means.append(np.mean(total_b_tubes))

	#now we calculate the standard error of the mean using bootstrapping
	
	length_means_boot = []
	for j in range(100):
		jode_constant_bootstrap = joining_ode_class.ODE_joining('constant', bootstrap = True)
		A_unjoined_raw_data, B_unjoined_raw_data, C_joined_raw_data, A_joined_raw_data, B_joined_raw_data = jode_constant_bootstrap.read_ABC_lengths(timepoint_prefixes[i])
		total_b_tubes_boot = B_unjoined_raw_data + B_joined_raw_data
		length_means_boot.append(np.mean(total_b_tubes_boot))
		#print np.mean(total_b_tubes_boot)
	standard_error_of_mean = scipy.stats.sem(length_means_boot)
	length_standard_error_of_mean.append(standard_error_of_mean)

print length_means
print length_standard_error_of_mean

plt.errorbar(timepoints, length_means, yerr = length_standard_error_of_mean)
plt.xlabel('time (hours)')
plt.ylabel('mean nanotube length µm')
plt.title('minority species nanotube length does not change with time')
plt.xlim(0.0, 5.0)
plt.ylim(4.0, 6.0)
plt.savefig(plot_dir_name+"nt_length_vs_time.pdf")
plt.close()



###### new ######
########################################################################################################
#plotting experimental excess A tube length distribution at each time point, overlaid with error bars


timepoint_prefixes = ['0hr', '2hr', '4hr']
timepoint_labels = ['.5 hrs', '2.5 hrs', '4.5 hrs']
colors = ['red','blue','green']


for timepoint_prefix in timepoint_prefixes:
	A_unjoined, B_unjoined, A_joined, B_joined, C_joined = model.read_ABC_concentrations(timepoint_prefix)

	lengths = [float((bin+1))*jode_constant.bin_length for bin in range(len(A_unjoined))]
	expt_A_data_timepoint = [A_unjoined[k] / sum(A_unjoined) for k in range(len(A_unjoined))] #normalizing
	#expt_A_data_timepoint = [A_unjoined[k] for k in range(len(A_unjoined))] #unnormalized
	density = gaussian_kde(expt_A_data_timepoint)
	f = UnivariateSpline(lengths, expt_A_data_timepoint, s=.05)
	#plt.plot(lengths, density(lengths), width = jode_constant.bin_length, alpha = .5, color = colors[timepoint_prefixes.index(timepoint_prefix)], label = timepoint_labels[timepoint_prefixes.index(timepoint_prefix)])
	plt.errorbar(lengths, f(lengths), yerr = [.013,.05,.048,.047,.022,.014,.017,.018,.014,.0098,0,0,0,0,0,0,0,0,0,0] ,color = colors[timepoint_prefixes.index(timepoint_prefix)], label = timepoint_labels[timepoint_prefixes.index(timepoint_prefix)])


plt.xlabel('length µm')
plt.ylabel('density')
#plt.title('minority species nanotube length does not change with time')
plt.xlim(0.0, 10.0)
plt.ylim(0.0,.30)
plt.legend(loc = 'upper right', fontsize = 16)
#plt.ylim(4.0, 6.0)
plt.savefig(plot_dir_name+"excess_distribution_does_not_change.pdf")
plt.close()

########################################################################################################
#plotting simulated vs experimental length distributions for each variable of interest at each timepoint
#TODO: normalize both distributions for fair comparison
#using subplot to make a massive 3x3 figure for each model and each timepoint 

f_unjoined_B, axarr_unjoined_B = plt.subplots(3,3, sharex = True, sharey = True, figsize = (9,9))
f_C, axarr_C = plt.subplots(3,3, sharex = True, sharey = True, figsize = (9,9))
f_joined_B, axarr_joined_B = plt.subplots(3,3, sharex = True, sharey = True, figsize = (9,9))

#f_unjoined_B.text(0.5, 0.04, 'length (µm)', ha='center')
#f_unjoined_B.text(0.04, 0.5, 'fraction', va='center', rotation='vertical')

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


		#axarr_unjoined_B[i,j].plot(lengths, sim_B_data_timepoint, 'red', label = 'simulation')
		#axarr_unjoined_B[i,j].plot(lengths, expt_B_data_timepoint, 'blue', label = 'experiment')
		axarr_unjoined_B[i,j].bar(lengths, sim_B_data_timepoint, width = model.bin_length, alpha = .5, color = 'red', label = 'simulation')
		axarr_unjoined_B[i,j].bar(lengths, expt_B_data_timepoint, width = model.bin_length, alpha = .5, color = 'blue', label = 'experiment')
		axarr_unjoined_B[i,j].set_xlim([0,2*model.n_bins])
		if i == j == 0:
			axarr_unjoined_B[i,j].legend(loc = 'upper right', fontsize = 8)
			axarr_unjoined_B[i,j].set_title("t = .5 hrs")
			axarr_unjoined_B[i,j].set_ylabel("constant kjoin")


		if i == 0 and j == 1:
			axarr_unjoined_B[i,j].set_title("t = 2.5 hrs")

		if i == 0 and j == 2:
			axarr_unjoined_B[i,j].set_title("t = 4.5 hrs")

		if j == 0 and i == 1:
			axarr_unjoined_B[i,j].set_ylabel("Hill model")

		if j == 0 and i == 2:
			axarr_unjoined_B[i,j].set_ylabel("Bernie model")



		#axarr_unjoined_B[i,j].xlabel('tube length (µm)') 
		#axarr_unjoined_B[i,j].ylabel('fraction')
		#axarr_unjoined_B[i,j].set_title(model.model_name+' unjoined B length distribution at '+timepoint_titles[j])
		#axarr_unjoined_B.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_' + timepoint_titles[j]+ '_expt_v_sim_unjoined_B_tubes.pdf')
		#axarr_unjoined_B.close()

		#axarr_C[i,j].bar(lengths, sim_C_data_timepoint, 'red', label = 'simulation')
		#axarr_C[i,j].bar(lengths, expt_C_data_timepoint, 'blue', label = 'experiment')
		axarr_C[i,j].bar(lengths, sim_C_data_timepoint, width = model.bin_length, alpha = .5, color = 'red', label = 'simulation')
		axarr_C[i,j].bar(lengths, expt_C_data_timepoint, width = model.bin_length, alpha = .5, color = 'blue', label = 'experiment')
		axarr_C[i,j].set_xlim([0,2*model.n_bins])
		if i == j == 0:
			axarr_C[i,j].legend(loc = 'upper right', fontsize = 8)
			axarr_C[i,j].set_title("t = .5 hrs")
			axarr_C[i,j].set_ylabel("constant kjoin")


		if i == 0 and j == 1:
			axarr_C[i,j].set_title("t = 2.5 hrs")

		if i == 0 and j == 2:
			axarr_C[i,j].set_title("t = 4.5 hrs")

		if j == 0 and i == 1:
			axarr_C[i,j].set_ylabel("Hill model")

		if j == 0 and i == 2:
			axarr_C[i,j].set_ylabel("Bernie model")
		
		#axarr_C[i,j].xlabel('tube length (µm)') 
		#axarr_C[i,j].ylabel('fraction')
		#axarr_C[i,j].set_title(model.model_name+' C length distribution at '+timepoint_titles[j])
		#axarr_C.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_' + timepoint_titles[j]+ '_expt_v_sim_C_tubes.pdf')
		#axarr_C.close()

		#axarr_joined_B[i,j].plot(lengths, sim_B_joined_data_timepoint, 'red', label = 'simulation')
		#axarr_joined_B[i,j].plot(lengths, expt_B_joined_data_timepoint, 'blue', label = 'experiment')
		axarr_joined_B[i,j].bar(lengths, sim_B_joined_data_timepoint, width = model.bin_length, alpha = .5, color = 'red', label = 'simulation')
		axarr_joined_B[i,j].bar(lengths, expt_B_joined_data_timepoint, width = model.bin_length, alpha = .5, color = 'blue', label = 'experiment')
		axarr_joined_B[i,j].set_xlim([0,2*model.n_bins])
		if i == j == 0:
			axarr_joined_B[i,j].legend(loc = 'upper right', fontsize = 8)
			axarr_joined_B[i,j].set_title("t = .5 hrs")
			axarr_joined_B[i,j].set_ylabel("constant kjoin")


		if i == 0 and j == 1:
			axarr_joined_B[i,j].set_title("t = 2.5 hrs")

		if i == 0 and j == 2:
			axarr_joined_B[i,j].set_title("t = 4.5 hrs")

		if j == 0 and i == 1:
			axarr_joined_B[i,j].set_ylabel("Hill model")

		if j == 0 and i == 2:
			axarr_joined_B[i,j].set_ylabel("Bernie model")
		
		#axarr_joined_B[i,j].xlabel('tube length (µm)') 
		#axarr_joined_B[i,j].ylabel('fraction')
		#axarr_joined_B[i,j].set_title(model.model_name+' joined B length distribution at '+timepoint_titles[j])
		#axarr_joined_B.savefig(plot_dir_name+str(model.model_name) + '_' + str(model.n_bins) + '_' + timepoint_titles[j]+ '_expt_v_sim_joined_B_tubes.pdf')
		#axarr_joined_B.close()

		

f_unjoined_B.set_tight_layout(True)
f_unjoined_B.savefig(plot_dir_name + '' + str(model.n_bins) + '_unjoined_B_subplot.pdf')

f_joined_B.set_tight_layout(True)
f_joined_B.savefig(plot_dir_name + '' + str(model.n_bins) + '_joined_B_subplot.pdf')

f_C.set_tight_layout(True)
f_C.savefig(plot_dir_name + '' + str(model.n_bins) + '_C_subplot.pdf')

plt.close()

#f_unjoined_B.close()
#f_joined_B.close()
#f_C.close()
########################################################################################################






				

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