import numpy as np
from scipy.integrate import odeint
import math 


class ODE_joining:
	def __init__(self, joining_model = 'constant', n_bins = 83, max_tube_length = 10.0, bootstrap = False):
		self.joining_model = joining_model
		self.n_bins = n_bins
		self.max_tube_length = max_tube_length
		self.bin_length = max_tube_length/float(n_bins)
		self.model_name = joining_model
		self.bootstrap = bootstrap
		self.n_slides = 25.0
		self.dilution_factor = 5.0 
		np.random.seed()
		self.random_seed =  np.random.randint(0,2**32 - 1) #generating a random seed here, each time read_ABC_concentrations is called the same random seed will be provided so identical calls should return the same values
		#print 'constructor random seed is: ' + str(self.random_seed)

	#python script to setup and integrate ODEs to model nanotube joining
	#we need ODEs for A,B, and C tubes: one ODE per bin

	def read_initial_concentrations(self):
		#read initial concentrations of A and B tubes from files
		#although there is some joining at the "0hr" timepoint we are going to asssume everything is unjoined 
		#for the purposes of running the ODE, therefore we will sum the joined and unjoined tubes for both A and B 
		time_prefix = '0hr'
		A_unjoined, B_unjoined, A_joined, B_joined, C_joined = self.read_ABC_concentrations(time_prefix)
		A_total = [sum(x) for x in zip(A_unjoined, A_joined)]
		B_total = [sum(x) for x in zip(B_unjoined, B_joined)]
		return A_total, B_total

	def read_ABC_concentrations(self, time_prefix):
		np.random.seed(self.random_seed) #resetting the state of the random number generator so we get the same concs on
										#sucesive calls of read_ABC_concentrations when bootstrapping is active

		#read in the A/B/C distributions for specified experimental timepoint 	
		A_unjoined_filename = time_prefix + '_unjoined647.dat'
		B_unjoined_filename = time_prefix + '_unjoinedcy3.dat'
		A_joined_filename = time_prefix + '_joinedcy3.dat'
		B_joined_filename = time_prefix + '_joinedcy3.dat'
		C_joined_filename = time_prefix + '_Cjoined.dat'
		

		with open(A_unjoined_filename) as f:
			content = f.readlines()
		A_unjoined_raw_data = [float(x.strip())*.17 for x in content]

		with open(B_unjoined_filename) as f:
			content = f.readlines()
		B_unjoined_raw_data = [float(x.strip())*.17 for x in content]

		with open(A_joined_filename) as f:
			content = f.readlines()
		A_joined_raw_data = [float(x.strip())*.17 for x in content]

		with open(B_joined_filename) as f:
			content = f.readlines()
		B_joined_raw_data = [float(x.strip())*.17 for x in content]

		with open(C_joined_filename) as f:
			content = f.readlines()
		C_joined_raw_data = [float(x.strip())*.17 for x in content]

		#n_bins = 10 #n_bins applies to the number of bins for A and B tubes, C tubes will have twice as many bins because they can be twice the length
		#max_tube_length = 25.0
		#self.bin_length =self.max_tube_length/float(n_bins)

		#print A_raw_data
		#print B_raw_data

		#now we bin the data for A/B tubes
		A_unjoined_count_binned = [0 for i in range(self.n_bins*2 ) ]
		B_unjoined_count_binned = [0 for i in range(self.n_bins*2 ) ]
		A_joined_count_binned = [0 for i in range(self.n_bins*2 ) ]
		B_joined_count_binned = [0 for i in range(self.n_bins*2 ) ]
		C_joined_count_binned = [0 for i in range(self.n_bins*2 ) ]

		bin_edges = [self.bin_length * (i) for i in range(self.n_bins + 1 )]
		bin_edges_C = [self.bin_length * (i) for i in range(self.n_bins*2 + 1 )]



		#here if we wish to bootstrap we will include only half of the raw data in the binning procedure
		#the half of the data that is selected will be randomized
		#modify this to take exactly 1/2 of the data every time rather than doing it randomly
		#this will ensure that the size of the sample is the same for each bootstrapping run 
		if self.bootstrap == True:
			raw_data = [A_unjoined_raw_data, B_unjoined_raw_data, C_joined_raw_data, A_joined_raw_data, B_joined_raw_data]
			#print "length before: "
			#print len(A_unjoined_raw_data)
			for data_type in raw_data:
				values_to_exclude = np.random.choice(data_type, len(data_type)/2, False)
				for value in values_to_exclude:
					data_type.remove(value)
			#print "length after: "
			#print len(A_unjoined_raw_data)



			#print B_joined_raw_data

		#because the lengths of A/B/C bin counts must match, A/B will be half populated with zeroes
		#calculating the bins here over the relevant range and then inserting into the actual list which is 2x the length
		A_unjoined_count_binned_pre, edges = np.histogram( A_unjoined_raw_data, bin_edges)
		B_unjoined_count_binned_pre, edges = np.histogram( B_unjoined_raw_data, bin_edges)
		A_joined_count_binned_pre, edges = np.histogram( A_joined_raw_data, bin_edges)
		B_joined_count_binned_pre, edges = np.histogram( B_joined_raw_data, bin_edges)
		C_joined_count_binned_pre, edges = np.histogram( C_joined_raw_data, bin_edges_C)


		for i in range(len(bin_edges) - 1):
			A_unjoined_count_binned[i] = A_unjoined_count_binned_pre[i]
			B_unjoined_count_binned[i] = B_unjoined_count_binned_pre[i]
			A_joined_count_binned[i] = A_joined_count_binned_pre[i]
			B_joined_count_binned[i] = B_joined_count_binned_pre[i]
		for i in range(len(bin_edges_C) - 1):
			C_joined_count_binned[i] = C_joined_count_binned_pre[i]


		#need to convert these to concentrations
		#subsequent concentrations at later time points will be normalized
		#first convert count to moles, then divide by 6 microliters, then divide by two for half on coverslip half on slide
		#divide by half the number of slides here if you are bootstrapping
		if self.bootstrap == True:
			effective_slides = self.n_slides/2.0
		else:
			effective_slides = self.n_slides 

		for i in range(len(bin_edges_C) - 1):
			#tried doing this conversion but missing the .087, why is this here? Using Deepak's calculation of N = 42[conc in pm]/per FOV
			#A_count_binned[i] = A_count_binned[i] * (1.0/6.022e23) * (1.0/6e-6) * (1.0/(.087*.087)) * 2.0
			#B_count_binned[i] = B_count_binned[i] * (1.0/6.022e23) * (1.0/6e-6) * 2.0
			
			A_unjoined_count_binned[i] = A_unjoined_count_binned[i] * (1.0/42.0) * (1.0/effective_slides) * 1e-12 * self.dilution_factor #1/25 per fov
			B_unjoined_count_binned[i] = B_unjoined_count_binned[i] * (1.0/42.0) * (1.0/effective_slides) * 1e-12 * self.dilution_factor
			A_joined_count_binned[i] = A_joined_count_binned[i] * (1.0/42.0) * (1.0/effective_slides) * 1e-12 * self.dilution_factor #1/25 per fov
			B_joined_count_binned[i] = B_joined_count_binned[i] * (1.0/42.0) * (1.0/effective_slides) * 1e-12 * self.dilution_factor
			C_joined_count_binned[i] = C_joined_count_binned[i] * (1.0/42.0) * (1.0/effective_slides) * 1e-12 * self.dilution_factor
			
		#print sum(A_count_binned), sum(B_count_binned)
		#print 'C_joined_count_binned'
		#print C_joined_count_binned
		#print B_joined_count_binned
		return A_unjoined_count_binned, B_unjoined_count_binned, A_joined_count_binned, B_joined_count_binned, C_joined_count_binned
		#this is not actually a count, it is a concentration...

	def calc_joined_B_dist_from_unjoined_B_dist(self, unjoined_B_dist, initial_B_dist):
		#given an unjoined B distribution from the ode solver, use the initial B distribution to back out the joined B distribution
		#the total concentration of B tubes should be conserved so this should be valid...
		joined_B_dist = [(x[0] - x[1]) for x in zip(initial_B_dist, unjoined_B_dist)]
		return joined_B_dist

	def unpack_concentrations_from_ode_format(self, w, t):
		#method to unpack the binned concentration data from the default output format used by the ode solver
		A_conc_data = []
		B_conc_data = []
		C_conc_data = []
		B_joined_conc_data = []

		initial_A, initial_B = self.read_initial_concentrations()

		for timepoint in range(len(t)):
			conc_data = w[timepoint]
			A_conc_data_timepoint = []
			B_conc_data_timepoint = []
			C_conc_data_timepoint = []

			i = 0 
			while i<len(conc_data):
				A_conc_data_timepoint.append(conc_data[i])
				B_conc_data_timepoint.append(conc_data[i+1])
				C_conc_data_timepoint.append(conc_data[i+2])
				i+=3
			A_conc_data.append(A_conc_data_timepoint)
			B_conc_data.append(B_conc_data_timepoint)
			C_conc_data.append(C_conc_data_timepoint)
			B_joined_conc_data.append(self.calc_joined_B_dist_from_unjoined_B_dist(B_conc_data_timepoint, initial_B))


		return t, A_conc_data, B_conc_data, C_conc_data, B_joined_conc_data

	def vectorfield(self, w, t, param):
		#for A tubes: dAi/dt = -1 * sum_over_m[kjoin(i,m)[Ai][Bm]]
		#for B tubes: dBi/dt = -1 * sum_over_m[kjoin(i,m)[Bi][Am]]
		#for C tubes: dCi/dt = 
		
		#n_bins = 10 #n_bins applies to the number of bins for A and B tubes, C tubes will have twice as many bins because they can be twice the length
		#max_tube_length = 25.0
		#self.bin_length =self.max_tube_length/n_bins

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

		dAdt_list = [0 for i in range( self.n_bins * 2) ]
		dBdt_list = [0 for i in range( self.n_bins * 2) ]
		dCdt_list = [0 for i in range( self.n_bins * 2) ]

		#f list is structured [dA0dt, dB0dt, dC0dt, dA1dt, dB1dt, dC1dt]
		#first we will populate the indivitual dAdt, dBdt, and dCdt lists according to the joining odes
		#now we will populate the f list of derivatives according to the joining odes

		#for A tubes: dAi/dt = -1 * sum_over_m[kjoin(i,m)[Ai][Bm]]
		#for B tubes: dBi/dt = -1 * sum_over_m[kjoin(i,m)[Bi][Am]]
		#for C tubes: dCi/dt = 
		
		for i in range(self.n_bins):
			dAdt = 0 
			dBdt = 0
			
			for j in range(self.n_bins):
				dAdt -= 1.0 * param * self.kjoin( (i+1)*self.bin_length, (j+1)*self.bin_length ) * Aconc[i] * Bconc[j]
				dBdt -= 1.0 * param * self.kjoin( (i+1)*self.bin_length, (j+1)*self.bin_length ) * Bconc[i] * Aconc[j]
				if i == j:
					dCdt = 1.0 * param * self.kjoin( (i+1)*self.bin_length, (j+1)*self.bin_length ) * Aconc[i] * Bconc[j] 
				else:
					dCdt = 1.0 * param * self.kjoin( (i+1)*self.bin_length, (j+1)*self.bin_length ) * Aconc[i] * Bconc[j] + 1.0 * param * self.kjoin( (i+1)*self.bin_length, (j+1)*self.bin_length ) * Bconc[i] * Aconc[j]

				dCdt_list[i+j] += dCdt

			dAdt_list[i] = dAdt
			dBdt_list[i] = dBdt

		#now we will populate the f list of derivatives according to the joining odes
		f = []

		for i in range(self.n_bins * 2):
			f.append(dAdt_list[i])
			f.append(dBdt_list[i])
			f.append(dCdt_list[i])

		return f







	def kjoin(self, n, m):
		#calculate kjoin according to type of joining model being considered
		#for bernie/hill models there is a scaling factor to be considered, this is related to the ratio of the tube length to diameter
		#going to set scaling to match Abdul's translational diffusion coefficient calculation
		#D = (kT/N*friction)*ln(N)
		#friction = 6*pi*viscosity*(diameter/2)
		#lengths must be divided by tube diameter!
		d = 0.015 #microns, lengths are in microns

		if self.joining_model == 'constant':
			return 1

		elif self.joining_model == 'hill':
			kjoin_hill = (n/d * math.log(m/d) + m/d * math.log(n/d)) / (n/d * m/d * (n/d+m/d))
			return kjoin_hill

		elif self.joining_model == 'bernie':
			kjoin_bernie = (n/d * math.log(m/d) + m/d * math.log(n/d)) / (n/d * m/d)
			return kjoin_bernie

		else:
			return None 

		

	def perform_integration(self, param):
		#perform the actual integration for a given value of kjoin


		
		A_count_binned, B_count_binned = self.read_initial_concentrations()
		C_count_binned = [0 for i in range(self.n_bins*2 ) ]

		#packaging up initial concentrations

		w = []
		#30 minutes = 1800 seconds
		#2hrs = 120 minutes = 7200 seconds + 30 minutes
		#4hrs = 240 minutes = 14400 seconds + 30 minutes 
		t = [0.0,1800.0, 9000.0, 16200.0]

		for i in range(self.n_bins*2):
			w.append(A_count_binned[i])
			w.append(B_count_binned[i])
			w.append(C_count_binned[i])

		#ODE solver parameters
		abserr = 1.0e-8
		relerr = 1.0e-6
		

		# Call the ODE solver.
		wsol = odeint(self.vectorfield, w, t, args=(param,), atol=abserr, rtol=relerr)

		#for t1, w1 in zip(t, wsol):
			#print t1, w1[0]

		return t, wsol

	def evaluate_error(self, param):
		#evaluate the error between the simulated distributions and the experimental ones
		#param has different meanings depending on the model
		#for the constant joining rate model, param is kjoin
		#for the hill and bernie models, param is the prefactor that the respective length dependence is multiplied by



		experimental_joining_percentage = [0.0, .30, .472, .51]
		simulated_joining_percentage = []

		#calculating simulated joining percentage
		#first perform the integration for the given kjoin
		t, wsol = self.perform_integration(param)
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
		
	def vahid_error(self, param):
		#evaluate the squared error based on vahid's cost function defined in his latest report
		#evaluate the error between the simulated distributions and the experimental ones
		#param has different meanings depending on the model
		#for the constant joining rate model, param is kjoin
		#for the hill and bernie models, param is the prefactor that the respective length dependence is multiplied by 

		#first we need to load in the experimental data for comparison
		#we are comparing at 30 minutes (labelled 0hr), 2 hrs (2.5hr), and 4hrs (4.5hr)
		initial_A, initial_B = self.read_initial_concentrations()
		t0_A_unjoined, t0_B_unjoined, t0_A_joined, t0_B_joined, t0_C_joined = self.read_ABC_concentrations('0hr')
		t2_A_unjoined, t2_B_unjoined, t2_A_joined, t2_B_joined, t2_C_joined = self.read_ABC_concentrations('2hr')
		t4_A_unjoined, t4_B_unjoined, t4_A_joined, t4_B_joined, t4_C_joined = self.read_ABC_concentrations('4hr')

		#calculating simulated joining percentage
		#first perform the integration for the given kjoin

		#original: experimental_joining_percentage = [0.0, .23, .339, .381]
		#experimental_joining_percentage = [0.0, .3, .4, .5]
		#experimental_joining_percentage is now calculated from the data directly rather than hard coded, this will allow this component
		#of the error function to properly handle bootstrapping now
		experimental_joining_percentage = [0.0]
		simulated_joining_percentage = []

		experimental_C_joined = [t0_C_joined, t2_C_joined, t4_C_joined]
		simulated_C_joined = []

		experimental_B_joined = [t0_B_joined, t2_B_joined, t4_B_joined]
		simulated_B_joined = []

		experimental_B_unjoined = [t0_B_unjoined, t2_B_unjoined, t4_B_unjoined]

		for i in range(3): 
			experimental_joining_percentage.append(sum(experimental_B_joined[i])/(sum(experimental_B_joined[i])+sum(experimental_B_unjoined[i])))

		print 'expt joining percentage: ', experimental_joining_percentage

		t, wsol = self.perform_integration(param)
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

			simulated_C_joined.append(Cconc)
			B_joined_conc = self.calc_joined_B_dist_from_unjoined_B_dist(Bconc, initial_B)
			simulated_B_joined.append(B_joined_conc)
			#simulated_B_joined.append(Bconc) #this is wrong :( this is actually the concentration of unjoined B tubes 

		joining_percentage_component = 0 
		Ctubes_component = 0
		Btubes_component = 0	
		squared_error = 0 
		for i in range(len(experimental_joining_percentage)):

			#this is the joining percentage contribution 
			joining_percentage_component += (experimental_joining_percentage[i]-simulated_joining_percentage[i]) * (experimental_joining_percentage[i]-simulated_joining_percentage[i])
			squared_error += 1.0 * (experimental_joining_percentage[i]-simulated_joining_percentage[i]) * (experimental_joining_percentage[i]-simulated_joining_percentage[i])

			
			if i>=1:
				#this is the C tubes cdf component
				Ctubes_component += self.cumulative_distribution_error(experimental_C_joined[i-1], simulated_C_joined[i])
				#squared_error += self.cumulative_distribution_error(experimental_C_joined[i-1], simulated_C_joined[i])

				#this is the joined B tubes cdf component
				Btubes_component += self.cumulative_distribution_error(experimental_B_joined[i-1], simulated_B_joined[i])
				#squared_error += self.cumulative_distribution_error(experimental_B_joined[i-1], simulated_B_joined[i])

		#print simulated_joining_percentage, experimental_joining_percentage
		#print "joining percentage component" + str(joining_percentage_component)
		#print "Ctubes component " + str(Ctubes_component)
		#print "Btubes component " + str(Btubes_component)
		return squared_error
			
	def cumulative_distribution_error(self, binned_data_1, binned_data_2):
		#given two sets of binned data, calculate the total squared error in the cumulative distributions
		#both binned datasets must have the same number of elements
		#print binned_data_1
		#print binned_data_2
		binned_data_1_normalization = float(sum(binned_data_1))
		binned_data_2_normalization = float(sum(binned_data_2))

		squared_error = 0.0
		cumulative_total_1 = 0.0
		cumulative_total_2 = 0.0


		for i in range(len(binned_data_1)):
			cumulative_total_1 += float(binned_data_1[i])
			cumulative_total_2 += float(binned_data_2[i])

			squared_error += ((cumulative_total_1/binned_data_1_normalization) - (cumulative_total_2/binned_data_2_normalization)) * ((cumulative_total_1/binned_data_1_normalization) - (cumulative_total_2/binned_data_2_normalization))

		return squared_error