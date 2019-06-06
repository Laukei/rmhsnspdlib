import scipy.fftpack as fft
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt

from pulseshapegenerator import PulseGenerator, create_dataset

class _BaseFFT:
	'''
	Base class for converting signal to frequency domain
	'''
	def __init__(self,dataset,timestep,*args,**kwargs):
		self.dataset = dataset
		self.timestep = timestep
		self._fft()


	def _fft(self):
		self._fdata = fft.rfft(self.dataset)
		self._freqs = fft.rfftfreq(self._fdata.size,d=self.timestep)
		# plt.plot(freqs,fdata)
		# plt.show()
		return self._freqs,self._fdata


	def _ifft(self):
		self._idata = fft.irfft(self._fdata)
		return self._idata


	def get_data(self):
		return self._ifft()


class BandPassFilter(_BaseFFT):
	'''
	band pass filter blocking all frequencies below low_frequency and above high_frequency
	frequencies supplied in Hz
	'''
	def __init__(self,dataset,timestep,*args,**kwargs):
		super().__init__(dataset,timestep)
		self.low_frequency = kwargs.get('low_frequency',0)
		self.high_frequency = kwargs.get('high_frequency',float('inf'))
		self._filter()


	def _filter(self):
		self._fdata = [f if (self._freqs[i] >= self.low_frequency and self._freqs[i] <= self.high_frequency) else 0 for i,f in enumerate(self._fdata)]


class Amplifier(_BaseFFT):
	'''
	Amplifier
	'''
	def __init__(self,dataset,timestep,filename,*args,**kwargs):
		super().__init__(dataset,timestep)
		self.filename = filename
		self._get_amplifier_data()
		self._amplify()


	def _get_amplifier_data(self):
		self.amp_data = np.genfromtxt(self.filename,skip_header=18,skip_footer=1,delimiter=',')
		self.amp_data = np.transpose(self.amp_data)
		self._update_lookup()


	def _amplify(self):
		self._fdata = self._fdata * (10 ** (self.lookup_gain(self._freqs)/10.0))


	def _update_lookup(self):
		self.lookup_gain = interp1d(self.amp_data[0],self.amp_data[1],bounds_error=False,fill_value='extrapolate')




if __name__ == "__main__":
	time_factor = 10 #datapoints per ns
	p = PulseGenerator(t_factor=time_factor)
	dataset = np.transpose(create_dataset(p,70000,125))[0]

	# b = BandPassFilter(dataset,1e-9/time_factor,low_frequency=10e6,high_frequency=580e6)
	# data = b.get_data()

	# plt.plot(range(len(b.dataset)),b.dataset)
	# plt.plot(range(len(data)),data)
	# plt.show()
	
	x_values = np.arange(len(dataset))/time_factor
	plt.plot(x_values,dataset,label='Simulated pulse')

	for filename in ['A1.csv','A2.csv','A3.csv','C1.csv']:
		a = Amplifier(dataset,1e-9/time_factor,'fieldfox\\'+filename)
		data = a.get_data()
		plt.plot(x_values,data,label=filename)
	plt.legend()
	plt.grid()
	plt.xlabel('Time (ns)')
	plt.ylabel('Relative voltage normalized to peak input')
	plt.show()