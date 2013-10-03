import numpy as np

eta = 0.5  # eta = 1 = perfect detection
mu = 1 		# number of incident photons
n = 4 		# number of pixels
b = 1		# number of pixels that need to fire for a detection event

eta = float(eta)
mu = float(mu)
n = float(n)
b = float(b)

def p_det(eta,n,i):			#detection probability for n pixels and eta, with i dead pixels
	return (eta*((n-i)/n))

def c_det(eta,n,i):			#the complement of p_det
	return 1 - p_det(eta,n,i)

def output_prefactor(eta,b,n): #calculates the probability of b pixels firing with no failures
	prefactor = 1
	for i in range(0,int(b)):
		prefactor = prefactor * p_det(eta,n,i)
	return prefactor

print 'prefactor:',output_prefactor(eta,b,n)

def recursive_summation(eta,b,n,mu):
	listOfIndividualErrors = []
	for k in range(0,int(b)):
		listOfIndividualErrors.append(c_det(eta,n,k))
	recursiveSum = 0
	for k in range(0,int(mu)-int(b)+1):
		recursiveSum+=combine(listOfIndividualErrors,k)
	return recursiveSum

def p_output(eta,b,n,mu):
	final_result = output_prefactor(eta,b,n) * (recursive_summation(eta,b,n,mu))
	return final_result

def combine(listOfIndividualErrors,levels):
	listOfIndividualErrors = np.array(listOfIndividualErrors)
	currentSet = listOfIndividualErrors
	#print previousSet
	if levels == 0 : #this must be 1 if there are 0 levels of combinations
		return 1
	for level in range(levels-1):
		newSet = np.array([0.0]*len(listOfIndividualErrors))
		for ir, i in enumerate(listOfIndividualErrors):
			for jr, j in enumerate(currentSet[ir:]):
				#print ir,jr,'i',i,'j',j
				newSet[ir] += i*j
		currentSet = newSet
	return sum(currentSet)


#print output_prefactor(eta,b,n)*(combine([0.5,6.25,7.5],0)+combine([0.5,0.625,0.75],1)+combine([0.5,0.625,0.75],2))
print p_output(eta,b,n,mu)
#print recursive_summation(eta,b,n,mu)