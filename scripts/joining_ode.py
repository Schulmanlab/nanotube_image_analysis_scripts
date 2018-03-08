import numpy as np
from scipy.integrate import odeint

#python script to setup and integrate ODEs to model nanotube joining
#we need ODEs for A,B, and C tubes: one ODE per bin
def read_initial_concentrations():
	#read initial concentrations of A and B tubes from files 
	Aconc_initial = []
	Bconc_initial = []
	
	A_data_filename = '0hr_unjoined647.dat'
	B_data_filename = '0hr_unjoinedcy3.dat'

	with open(A_data_filename) as f:
		content = f.readlines()
	A_raw_data = [float(x.strip())*.17 for x in content]

	with open(B_data_filename) as f:
		content = f.readlines()
	B_raw_data = [float(x.strip())*.17 for x in content]

	n_bins = 10 #n_bins applies to the number of bins for A and B tubes, C tubes will have twice as many bins because they can be twice the length
	max_tube_length = 25.0
	bin_length = max_tube_length/float(n_bins)

	#print A_raw_data
	#print B_raw_data

	#now we bin the data for A/B tubes
	A_count_binned = [0 for i in range( n_bins*2 ) ]
	B_count_binned = [0 for i in range( n_bins*2 ) ]
	bin_edges = [bin_length * (i) for i in range( n_bins )]

	print bin_edges

	#because the lengths of A/B/C bin counts must match, A/B will be half populated with zeroes
	#calculating the bins here over the relevant range and then inserting into the actual list which is 2x the length
	A_count_binned_pre, edges = np.histogram( A_raw_data, bin_edges)
	B_count_binned_pre, edges = np.histogram( B_raw_data, bin_edges)


	for i in range(len(bin_edges)-1):
		A_count_binned[i] = A_count_binned_pre[i]
		B_count_binned[i] = B_count_binned_pre[i]

	
	return A_count_binned, B_count_binned

	#for tube in A_raw_data:
	#	for i in range(len(A_count_binned)):
	#		if tube > bin_upper_limits


def vectorfield(w, t, p):
	#for A tubes: dAi/dt = -1 * sum_over_m[kjoin(i,m)[Ai][Bm]]
	#for B tubes: dBi/dt = -1 * sum_over_m[kjoin(i,m)[Bi][Am]]
	#for C tubes: dCi/dt = 
	
	n_bins = 10 #n_bins applies to the number of bins for A and B tubes, C tubes will have twice as many bins because they can be twice the length
	max_tube_length = 25.0
	bin_length = max_tube_length/n_bins

	Aconc = []
	Bconc = []
	Cconc = []

	#w list is structured [A0, B0, C0, A1, B1, B2]
	#populate concentration arrays
	i = 0 
	while i<len(w):
		Aconc.append(w[i])
		Bconc.append(w[i+1])
		Cconc.append(w[i+2])
		i+=3

	dAdt_list = [0 for i in range( n_bins * 2) ]
	dBdt_list = [0 for i in range( n_bins * 2) ]
	dCdt_list = [0 for i in range( n_bins * 2) ]

	#f list is structured [dA0dt, dB0dt, dC0dt, dA1dt, dB1dt, dC1dt]
	#first we will populate the indivitual dAdt, dBdt, and dCdt lists according to the joining odes
	#now we will populate the f list of derivatives according to the joining odes

	#for A tubes: dAi/dt = -1 * sum_over_m[kjoin(i,m)[Ai][Bm]]
	#for B tubes: dBi/dt = -1 * sum_over_m[kjoin(i,m)[Bi][Am]]
	#for C tubes: dCi/dt = 
	
	for i in range(n_bins):
		dAdt = 0 
		dBdt = 0
		
		for j in range(n_bins):
			dAdt -= 1.0 * kjoin( (i+1)*bin_length, (j+1)*bin_length ) * Aconc[i] * Bconc[j]
			dBdt -= 1.0 * kjoin( (i+1)*bin_length, (j+1)*bin_length ) * Bconc[i] * Aconc[i]
			dCdt = 1.0 * kjoin( (i+1)*bin_length, (j+1)*bin_length ) * Aconc[i] * Bconc[j] + 1.0 * kjoin( (i+1)*bin_length, (j+1)*bin_length ) * Bconc[i] * Aconc[i]
			dCdt_list[i+j] += dCdt

		dAdt_list[i] = dAdt
		dBdt_list[i] = dBdt

	#now we will populate the f list of derivatives according to the joining odes
	f = []

	for i in range(n_bins):
		f.append(dAdt_list[i])
		f.append(dBdt_list[i])
		f.append(dCdt_list[i])







def kjoin(n, m):
	#calculate kjoin according to type of joining model being considered
	return 1

read_initial_concentrations()




