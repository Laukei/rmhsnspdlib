#
#   A quick and hacky program that fits calculateCountRate's formula to experimental_counts' data
#    using a least-squares fit approach, and minimizing in 3 dimensions.
#
#                   by Rob Heath
#

import numpy as np

input_fluxes = [15.45842,38.8298,97.53604,244.99946,615.41081,1545.84206,3882.9797,9753.60401,24499.94558,61541.08088,154584.20604,388297.96969,975360.40144,2.45E+006,6.15E+006,1.55E+007,3.88E+007,9.75E+007,2.45E+008,6.15E+008,1.55E+009,3.88E+009,9.75E+009,2.45E+010,6.15E+010]
#experimental_counts = [0,0,0,0,0,1,0,2,2,2,2,3,6,13,56,218,1253,6358,29384,101638,285965,656111,962296,998279,1.01E+006]
experimental_counts = [18,27,32,22,38,32,62,31,24,47,49,69,77,177,452,1029,2566,6988,18063,40308,84666,207935,481233,866281,995885]


def calculateCountRate(eta,D,exponent,input_fluxes):
	f = 1000000 #frequency of the laser
	calculatedCounts = []
	for flux in input_fluxes:
		#flux/f = mu
		#exponent accounts for multiple photon absorption
		#D accounts for dark-count rate
		calculatedCounts.append((f*(1-np.exp(-(flux/f)*eta))**exponent)+D) 
	return calculatedCounts

def returnResiduals(calculated_counts, experimental_counts):
	sum_of_residuals = 0
	for c, count in enumerate(calculated_counts):
		sum_of_residuals+=(experimental_counts[c]-calculated_counts[c])**2
	return sum_of_residuals

def generateGridOfResiduals(D_range,eta_range,exponent_range,input_fluxes,experimental_counts):
	residuals = []
	minima = [float("inf"),0,0,0]

	for d, D in enumerate(D_range):
		residuals.append([])
		print "On D:\t",d,"of",len(D_range)
		for e, eta in enumerate(eta_range):
			residuals[d].append([])
			for x, exponent in enumerate(exponent_range):
				calculated_counts = calculateCountRate(eta,D,exponent,input_fluxes)
				residual = returnResiduals(calculated_counts,experimental_counts)
				residuals[d][e].append(residual)
				if residual < minima[0]:
					minima[0], minima[1], minima[2], minima[3] = residual, d, e, x
	return residuals, minima

def checkForConvergence(residuals,minima,D_range,eta_range,exponent_range):
	d, e, x = minima[1], minima[2], minima[3]
	success = [False,False,False]
	try:
		if residuals[d][e][x-1] > minima[0] and residuals[d][e][x+1] > minima[0]:
			print "Exponent is optimal"
			success[0] = True
		else:
			print "Exponent NOT OPTIMAL: no convergence"
	except IndexError:
		print "Exponent NOT OPTIMAL: range does not permit convergence"

	try:
		if residuals[d][e-1][x] > minima[0] and residuals[d][e+1][x] > minima[0]:
			print "Eta is optimal"
			success[1] = True
		else:
			print "Eta NOT OPTIMAL: no convergence"
	except IndexError:
		print "Eta NOT OPTIMAL: range does not permit convergence"

	try:
		if residuals[d-1][e][x] > minima[0] and residuals[d+1][e][x] > minima[0]:
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
	print "with minima residual value: ", minima[0]

#D_range = [2] #from x, to y, spaced as Z
D_range = np.linspace(10,30,5)
eta_range = np.linspace(0,0.01,1000)
exponent_range = np.linspace(0.8,1.4,61)
#exponent_range = [1]
#D needs locking in place - the weights at the far end excessively outweigh the near end

residuals, minima = generateGridOfResiduals(D_range,eta_range,exponent_range,input_fluxes,experimental_counts)
checkForConvergence(residuals,minima,D_range,eta_range,exponent_range)

