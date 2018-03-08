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

	#print bin_edges

	#because the lengths of A/B/C bin counts must match, A/B will be half populated with zeroes
	#calculating the bins here over the relevant range and then inserting into the actual list which is 2x the length
	A_count_binned_pre, edges = np.histogram( A_raw_data, bin_edges)
	B_count_binned_pre, edges = np.histogram( B_raw_data, bin_edges)


	for i in range(len(bin_edges)-1):
		A_count_binned[i] = A_count_binned_pre[i]
		B_count_binned[i] = B_count_binned_pre[i]

	#need to convert these to concentrations
	#subsequent concentrations at later time points will be normalized
	#first convert count to moles, then divide by 6 microliters, then divide by two for half on coverslip half on slide
	for i in range(len(bin_edges)-1):
		#tried doing this conversion but missing the .087, why is this here? Using Deepak's calculation of N = 42[conc in pm]/per FOV
		#A_count_binned[i] = A_count_binned[i] * (1.0/6.022e23) * (1.0/6e-6) * (1.0/(.087*.087)) * 2.0
		#B_count_binned[i] = B_count_binned[i] * (1.0/6.022e23) * (1.0/6e-6) * 2.0
		A_count_binned[i] = A_count_binned[i] * (1.0/42.0) * (1.0/25.0) * 1e-12 * 5.0#1/25 per fov
		B_count_binned[i] = B_count_binned[i] * (1.0/42.0) * (1.0/25.0) * 1e-12 * 5.0

	#print sum(A_count_binned), sum(B_count_binned)

	return A_count_binned, B_count_binned

	




def vectorfield(w, t, kjoin):
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
			dAdt -= 1.0 * kjoin * Aconc[i] * Bconc[j]
			dBdt -= 1.0 * kjoin * Bconc[i] * Aconc[i]
			dCdt = 1.0 * kjoin * Aconc[i] * Bconc[j] + 1.0 * kjoin * Bconc[i] * Aconc[i]
			dCdt_list[i+j] += dCdt

		dAdt_list[i] = dAdt
		dBdt_list[i] = dBdt

	#now we will populate the f list of derivatives according to the joining odes
	f = []

	for i in range(n_bins * 2):
		f.append(dAdt_list[i])
		f.append(dBdt_list[i])
		f.append(dCdt_list[i])

	return f







'''def kjoin(n, m):
	#calculate kjoin according to type of joining model being considered
	return 1
	'''

def perform_integration(kjoin):
	#perform the actual integration for a given value of kjoin

	n_bins = 10 #n_bins applies to the number of bins for A and B tubes, C tubes will have twice as many bins because they can be twice the length
	max_tube_length = 25.0
	bin_length = max_tube_length/n_bins

	
	A_count_binned, B_count_binned = read_initial_concentrations()
	C_count_binned = [0 for i in range( n_bins*2 ) ]

	#packaging up initial concentrations

	w = []
	#30 minutes = 1800 seconds
	#2hrs = 120 minutes = 7200 seconds + 30 minutes
	#4hrs = 240 minutes = 14400 seconds + 30 minutes 
	t = [0.0,1800, 9000.0, 16200.0]

	for i in range(n_bins*2):
		w.append(A_count_binned[i])
		w.append(B_count_binned[i])
		w.append(C_count_binned[i])

	#ODE solver parameters
	abserr = 1.0e-8
	relerr = 1.0e-6
	

	# Call the ODE solver.
	wsol = odeint(vectorfield, w, t, args=(kjoin,), atol=abserr, rtol=relerr)

	#for t1, w1 in zip(t, wsol):
		#print t1, w1[0]

	return t, wsol

def evaluate_error(kjoin):
	#evaluate the error between the simulated distributions and the experimental ones
	#in this version we will only evaluate the error in the percentage of joined B tubes



	experimental_joining_percentage = [0.0, .30, .472, .51]
	simulated_joining_percentage = []

	#calculating simulated joining percentage
	#first perform the integration for the given kjoin
	t, wsol = perform_integration(kjoin)
	for w in wsol:
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
		unjoined_B_tubes = sum(Bconc)
		joined_B_tubes = sum(Cconc)
		joining_fraction = joined_B_tubes/(joined_B_tubes+unjoined_B_tubes)
		simulated_joining_percentage.append(joining_fraction)

	squared_error = 0 
	for i in range(len(experimental_joining_percentage)):
		squared_error += (experimental_joining_percentage[i]-simulated_joining_percentage[i]) * (experimental_joining_percentage[i]-simulated_joining_percentage[i])

	#print simulated_joining_percentage, experimental_joining_percentage
	
	return squared_error
	

		

#here starts the main script to brute force try many constant joining values and evaluate the error with expt
#*abdul sees 3e6 per molar per second in his work*

kjoin_list = [1e6+(i*1e6) for i in range(100)]
for i in range(100):
	kjoin = kjoin_list[i]
	squared_error = evaluate_error(float(kjoin))
	print kjoin, squared_error

#evaluate_error(1000000.0)

#cost fucntion... 







