#
#
#    this variation has a shape factor A:
#
#       it generates a bxb matrix and rather than being bi-diagonal
#         is tri-diagonal, allowing jumps from n pixels triggered to n+2 pixels triggered
#           in a single absorption, which simulates a detector's non-discrete regimes.
#
#            This is controlled by parameter A, which is the bias towards one-pixel detection. 1-A is the probability of a two-pixel jump.
#          

import numpy as np
import matplotlib.pyplot as plt
import math

eta = 1e-3  # eta = 1 = perfect detection
n = 2		# number of pixels
b = 2		# number of pixels that need to fire for a detection event
f = 1000000 # frequency of input pulses
mu_max = 10000
step = 1

inputFluxes = np.array([15.45842,38.8298,97.53604,244.99946,615.41081,1545.84206,3882.9797,9753.60401,24499.94558,61541.08088,154584.20604,388297.96969,975360.40144,2.45E+006,6.15E+006,1.55E+007,3.88E+007,9.75E+007,2.45E+008,6.15E+008,1.55E+009,3.88E+009,9.75E+009,2.45E+010,6.15E+010])
experimentalCounts = np.array([0,0,0,0,0,1,0,2,2,2,2,3,6,13,56,218,1253,6358,29384,101638,285965,656111,962296,998279,1.01E+006])/f #two-photon

for i, flux in enumerate(inputFluxes):
	inputFluxes[i]=int(flux)/f
#print inputFluxes

def p_det(eta,n,i):			#detection probability for n pixels and eta, with i dead pixels
	#print (float(n)-float(i)),'/',n,'=',(n-i)/n
	return eta*((float(n)-float(i))/float(n))

def c_det(eta,n,i):			#the complement of p_det
	return 1.0 - p_det(eta,n,i)

def transition_matrix(eta,b,n):
	transMatrix = np.zeros((int(b)+1,int(b)+1))
	transMatrix[-1][-1] = 1.0
	for k in range(int(b)):
		transMatrix[k][k]=c_det(eta,n,k)
		transMatrix[k][k+1]=p_det(eta,n,k)
	return transMatrix

def initial_state(b):
	initState = np.zeros((int(b)+1))
	initState[0] = 1.0
	return initState

def p_output(eta,b,n,mu):
	muthMatrix = np.linalg.matrix_power(transition_matrix(eta,b,n),int(mu))
	muthState = np.dot(initial_state(b),muthMatrix)
	return muthState

def p_output_simple(eta,b,n,mu):
	muthState = 1.0
	for i in range(0,int(b)):
		muthState = muthState * (1-(1-(((n-i)/n)*eta))**mu)
	return muthState

def p_output_toosimple(eta,b,n,mu):
	muthState = ((math.factorial(n)/math.factorial(n-b))*((1/n)**b))*((1-((1-eta)**mu))**b)
	return muthState
#print np.dot(initial_state(float(b)),np.linalg.matrix_power(transition_matrix(float(eta),float(b),float(n)),2))

def mark_wtwo_btwo(eta,b,n,mu):
	muthState = (1 + (1-eta)**mu - 2*((1-(eta/2))**mu))
	return muthState

def slap_the_matrix(matrix):
	#slappyTotal = matrix[0][2]/(matrix[0][1]*matrix[1][2]) + matrix[1][2]
	slappyTotal = ((1-matrix[0][0])*(1-matrix[1][1]))/matrix[0][2]
	return slappyTotal



X,Y = [],[]
Ytwo = []
Yold = []
Yoldtwo = []
Ymark = []
Yslap = []
error = []
corrected = []
for mu in range(0,mu_max+1,step):
#for mu in inputFluxes:
	eta = float(eta)
	mu = float(mu)
	n = float(n)
	b = float(b)
	X.append(mu)
	prob = p_output(eta,b,n,mu)[-1]
	Y.append(prob)
	print 'mu:',mu,'probability:',Y[-1],'\r'
	B = (1-eta)**mu
	Ytwo.append(p_output_simple(eta,b,n,mu))
	#Yold.append((1.0 - (1-eta)**mu)**b)
	Yold.append(p_output_toosimple(eta,b,n,mu))
	Yoldtwo.append((1.0 - 2.718281828459**(-mu*eta))**b)
	print 'error',slap_the_matrix(np.linalg.matrix_power(transition_matrix(eta,b,n),int(mu)))
	Ymark.append(mark_wtwo_btwo(eta,b,n,mu))
	error.append(slap_the_matrix(np.linalg.matrix_power(transition_matrix(eta,b,n),int(mu))))
	corrected.append(Ytwo[-1]/error[-1])
	print 'ca:',mu,'probability:',Ytwo[-1]

#plt.loglog(inputFluxes,experimentalCounts,'ro')
plt.xlabel(r'Photons per pulse ($\mu$)')
plt.ylabel(r'Trigger probability $P_{output}$')
plt.title(r'$\eta$: '+str(eta)+', $n$: '+str(int(n))+', $b$: '+str(int(b)))
plt.grid(b=True)
plt.loglog(X,Y,'b-',label = r'Matrix')
plt.loglog(X,Ytwo,'r-',label = r'$\prod_{i=0}^{b-1}\left(1-\left(1-\frac{n-i}{n}\eta\right)^{\mu}\right)$')
#plt.loglog(X,Yoldtwo,'y--', label = r'$\left(1-e^{-\mu\eta}\right)^b$')
#plt.loglog(X,Yold,'g--',label = r'$\frac{n!}{(n-b)!} \left(\frac{1}{n}\right)^b \left(1-(1-\eta)^{\mu}\right)^b$')
plt.loglog(X,Ymark,'g--',label = r'$w=2$,$b=2$ (email)')
#plt.loglog(X,error,'r--',label = r'so-called error')
#plt.loglog(X,corrected,'g-', label = r'corrected')
plt.legend(loc='lower right', shadow=False)
#plt.loglog(X,Yoldtwo,'b--')
plt.show()