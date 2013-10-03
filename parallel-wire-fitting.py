import numpy as np
import matplotlib.pyplot as plt

eta = 0.0009  # eta = 1 = perfect detection
n = 2		# number of pixels
b = 2		# number of pixels that need to fire for a detection event
f = 1000000 # frequency of input pulses

inputFluxes = np.array([15.45842,38.8298,97.53604,244.99946,615.41081,1545.84206,3882.9797,9753.60401,24499.94558,61541.08088,154584.20604,388297.96969,975360.40144,2.45E+006,6.15E+006,1.55E+007,3.88E+007,9.75E+007,2.45E+008,6.15E+008,1.55E+009,3.88E+009,9.75E+009,2.45E+010,6.15E+010])
experimentalCounts = np.array([0,0,0,0,0,1,0,2,2,2,2,3,6,13,56,218,1253,6358,29384,101638,285965,656111,962296,998279,1.01E+006])/f #two-photon

for i, flux in enumerate(inputFluxes):
	inputFluxes[i]=int(flux)/f
print inputFluxes

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
	muthState = np.dot(initial_state(b),np.linalg.matrix_power(transition_matrix(eta,b,n),int(mu)))
	return muthState

#print np.dot(initial_state(float(b)),np.linalg.matrix_power(transition_matrix(float(eta),float(b),float(n)),2))

X,Y = [],[]
#for mu in range(0,mu_max+1,step):
for mu in inputFluxes:
	eta = float(eta)
	mu = float(mu)
	n = float(n)
	b = float(b)
	X.append(mu)
	prob = p_output(eta,b,n,mu)[-1]
	Y.append(prob)
	print 'mu:',mu,'probability:',Y[-1],'\r'

plt.loglog(inputFluxes,experimentalCounts,'ro')
plt.xlabel(r'Photons per pulse ($\mu$)')
plt.ylabel(r'Trigger probability $P_{output}$')
plt.title(r'$\eta$: '+str(eta)+', $n$: '+str(int(n))+', $b$: '+str(int(b)))
plt.grid(b=True)
plt.loglog(X,Y,'b-')
plt.show()