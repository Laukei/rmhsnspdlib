import numpy as np
import matplotlib.pyplot as plt

eta = 0.001  # eta = 1 = perfect detection
mu_max = 1000		# number of incident photons
step = 1    # number to increment mu_max by (to avoid incredible computation times...)
n = 1		# number of pixels
b = 1		# number of pixels that need to fire for a detection event


def p_det(eta,n,i):			#detection probability for n pixels and eta, with i dead pixels
	return (eta*((n-i)/n))

def c_det(eta,n,i):			#the complement of p_det
	return 1 - p_det(eta,n,i)

def transition_matrix(eta,b,n):
	transMatrix = np.zeros((int(b)+1,int(b)+1))
	transMatrix[-1][-1] = 1
	for k in range(int(b)):
		transMatrix[k][k]=c_det(eta,n,k)
		transMatrix[k][k+1]=p_det(eta,n,k)
	return transMatrix

def initial_state(b):
	initState = np.zeros((int(b)+1))
	initState[0] = 1
	return initState

def p_output(eta,b,n,mu):
	currentState = initial_state(b)
	for i in range(int(mu)):
		currentState = np.dot(currentState,transition_matrix(eta,b,n))
		#print currentState
	return currentState

print np.linalg.matrix_power(transition_matrix(eta,b,n),1000) ## this is hugely relevant!

X,Y = [],[]
for mu in range(0,mu_max+1,step):
	eta = float(eta)
	mu = float(mu)
	n = float(n)
	b = float(b)
	X.append(mu)
	Y.append(p_output(eta,b,n,mu)[-1])
	print 'mu:',mu,'probability:',Y[-1]
plt.xlabel(r'$\mu$')
plt.ylabel('Probability of output from detector')
plt.title(r'$\eta$: '+str(eta)+', $n$: '+str(int(n))+', $b$: '+str(int(b)))
plt.plot(X,Y) 
plt.show()