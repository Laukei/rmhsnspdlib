#
#   A quick and hacky program that fits calculateCountRate's formula to experimental_counts' data
#    using a least-squares fit approach, and minimizing in 3 dimensions.
#
#                   by Rob Heath
#

import numpy as np

input_fluxes = [15.45842,38.8298,97.53604,244.99946,615.41081,1545.84206,3882.9797,9753.60401,24499.94558,61541.08088,154584.20604,388297.96969,975360.40144,2.45E+006,6.15E+006,1.55E+007,3.88E+007,9.75E+007,2.45E+008,6.15E+008,1.55E+009,3.88E+009,9.75E+009,2.45E+010,6.15E+010]
experimental_counts = [0,0,0,0,0,1,0,2,2,2,2,3,6,13,56,218,1253,6358,29384,101638,285965,656111,962296,998279,1.01E+006] #two-photon
#experimental_counts = [18,27,32,22,38,32,62,31,24,47,49,69,77,177,452,1029,2566,6988,18063,40308,84666,207935,481233,866281,995885] #one-photon


def calculateCountRate(eta,D,exponent,input_fluxes):
	f = 1000000 #frequency of the laser
	calculatedCounts = []
	for flux in input_fluxes:
		#flux/f = mu
		#exponent accounts for multiple photon absorption
		#D accounts for dark-count rate
		calculatedCounts.append((f*(1-np.exp(-(flux/f)*eta))**exponent)+D) 
	return calculatedCounts

def returnRsquared(calculated_counts, experimental_counts):
	s_tot = 0
	s_reg = 0
	s_err = 0
	experimental_average = np.average(experimental_counts)
	for c, count in enumerate(calculated_counts):
		#based on: http://en.wikipedia.org/wiki/Coefficient_of_determination (03/04/2013 17:55)
		s_tot += (experimental_counts[c]-experimental_average)**2
		#s_reg += (calculated_counts[c]-experimental_average)**2
		s_err +=(experimental_counts[c]-calculated_counts[c])**2
	R_squared = 1 - (s_err/s_tot)
	return R_squared

def generateGridOfRsquareds(D_range,eta_range,exponent_range,input_fluxes,experimental_counts):
	R_squareds = []
	minima = [0,0,0,0]

	for d, D in enumerate(D_range):
		R_squareds.append([])
		print "On D:\t",d,"of",len(D_range)
		for e, eta in enumerate(eta_range):
			R_squareds[d].append([])
			for x, exponent in enumerate(exponent_range):
				calculated_counts = calculateCountRate(eta,D,exponent,input_fluxes)
				R_squared = returnRsquared(calculated_counts,experimental_counts)
				R_squareds[d][e].append(R_squared)
				if R_squared > minima[0]:
					minima[0], minima[1], minima[2], minima[3] = R_squared, d, e, x
	return R_squareds, minima

def checkForConvergence(R_squareds,minima,D_range,eta_range,exponent_range):
	d, e, x = minima[1], minima[2], minima[3]
	success = [False,False,False]
	try:
		if R_squareds[d][e][x-1] < minima[0] and R_squareds[d][e][x+1] < minima[0]:
			print "Exponent is optimal"
			success[0] = True
		else:
			print "Exponent NOT OPTIMAL: no convergence"
	except IndexError:
		print "Exponent NOT OPTIMAL: range does not permit convergence"

	try:
		if R_squareds[d][e-1][x] < minima[0] and R_squareds[d][e+1][x] < minima[0]:
			print "Eta is optimal"
			success[1] = True
		else:
			print "Eta NOT OPTIMAL: no convergence"
	except IndexError:
		print "Eta NOT OPTIMAL: range does not permit convergence"

	try:
		if R_squareds[d-1][e][x] < minima[0] and R_squareds[d+1][e][x] < minima[0]:
			print "D is optimal"
			success[2] = True
		else:
			print "D NOT OPTIMAL: no convergence"
	except IndexError:
		print "D NOT OPTIMAL: range does not permit convergence"

	if success[0] == success[1] == success[2] == True:
		print "All values optimized."
	else:
		print "Some values sub-optimal."

	print "Eta:\t", eta_range[e]
	print "D:\t", D_range[d]
	print "Exp:\t", exponent_range[x] 
	print "with R^2 value: ", minima[0]

D_range = [2] #from x, to y, spaced as Z
eta_range = np.linspace(0.0001,0.001,101)
exponent_range = np.linspace(1.2,1.8,61)
#D needs locking in place - the weights at the far end excessively outweigh the near end

R_squareds, minima = generateGridOfRsquareds(D_range,eta_range,exponent_range,input_fluxes,experimental_counts)
checkForConvergence(R_squareds,minima,D_range,eta_range,exponent_range)

