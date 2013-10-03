import numpy as np
import matplotlib.pyplot as plt

eta = 0.01  # eta = 1 = perfect detection
mu_max = 5		# number of incident photons
step = 1    # number to increment mu_max by (to avoid incredible computation times...)
n = 4 		# number of pixels
b = 3		# number of pixels that need to fire for a detection event

def p_det(eta,n,i):			#detection probability for n pixels and eta, with i dead pixels
	return (eta*((n-i)/n))

def c_det(eta,n,i):			#the complement of p_det
	return 1 - p_det(eta,n,i)

def output_prefactor(eta,b,n): #calculates the probability of b pixels firing with no failures
	prefactor = 1
	for i in range(0,int(b)):
		prefactor = prefactor * p_det(eta,n,i)
	return prefactor

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
		#print currentSet
	return sum(currentSet)


#print p_output(eta,b,n,mu)

X,Y = [],[]
for mu in range(0,mu_max+1,step):
	eta = float(eta)
	mu = float(mu)
	n = float(n)
	b = float(b)
	X.append(mu)
	Y.append(p_output(eta,b,n,mu))
	print 'mu:',mu,'probability:',Y[-1]
plt.xlabel(r'$\mu$')
plt.ylabel('Probability of output from detector')
plt.title(r'$\eta$: '+str(eta)+', $n$: '+str(int(n))+', $b$: '+str(int(b)))
plt.plot(X,Y) 
plt.show()